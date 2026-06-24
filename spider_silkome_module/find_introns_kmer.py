"""
find_introns_kmer.py — k-mer repeat-profile intron finder for full-length
spidroin gene sequences (sense strand, ATG start, stop end).

Why this method (and not the obvious ones):
  * Reading-frame / stop-codon scanning fails: the repeat region is codon
    biased and has no stop codon in any of the three frames, so frame breaks
    are invisible.
  * Bare GT..AG search is far too degenerate to localize an intron.

Core idea — repeats vs. non-repeats, frame-independent:
  1. Count every k-mer (k=14) within the sequence and assign each base the
     MINIMUM copy number among the k-mers covering it.  copy == 1 -> single
     copy (non-repeat); copy >= 2 -> repeat-region exon.
  2. Take contiguous single-copy runs of >= min_run nt as candidate regions
     (short scattered single-copy stretches are k-mer edge noise around point
     mutations and are ignored).  The 5'-most and 3'-most runs are the NTD and
     CTD domains, NOT introns, so they are dropped; only runs inside the repeat
     body are candidate introns.
  3. Flank-repeat filter: a real intron is flanked by the high-copy repeat
     body, so require the median k-mer copy of both flanks (flank_window nt on
     each side) to be >= min_flank_cov.  This removes single-copy stretches
     sitting at the edge of degenerate / low-copy regions.
  4. Boundaries: within +/- splice_search nt of each candidate end, search for
     GT (donor) .. AG (acceptor).  Among all pairs require the intron length
     ≡ frame_remainder (mod frame_mod) (domain rule: length ≡ 1 (mod 3)) and
     pick the pair whose length is closest to the candidate span.

Thresholds were calibrated on all 1660 full-length sequences under
data/processed/typing_results and are exposed as parameters (defaults = the
calibrated values); nothing is hard-coded.

Outputs (per species, under output_dir/<species>/):
  kmer_introns.gff      # intron rows on the spidroin_full_length sequences
                        # (1-based relative coords), aligned with the
                        # protein_confirmation `fullprot_miniprot.gff` schema so
                        # that agents/evidence/junction_validation.py can reuse it.
  kmer_coverage.tsv     # per-base coverage, one '# <spidroin_id>' block per
                        # sequence then columns position(1-based), base,
                        # kmer_coverage, and rna_coverage (4th column, only when a
                        # matching RNA-seq bigWig exists; disable --no-emit-coverage).
  spliced_proteins.faa  # one protein per sequence: the sense CDS with every
                        # predicted intron spliced out, translated in frame 0.
                        # An internal '*' (premature_stop=yes in the header) means
                        # the intron call is almost certainly wrong; disable with
                        # --no-emit-proteins.
  figures/<id>.png      # one coverage plot per sequence (x = position): k-mer
                        # coverage as a log-scaled filled area with introns / NTD /
                        # CTD shaded, plus the RNA-seq coverage as a green line on a
                        # secondary axis (disable with --no-emit-plots).
Cross-species `intron_summary.tsv` and `protein_summary.tsv` are also written
under output_dir; the output_dir itself carries a run timestamp suffix.

RNA-seq bigWigs are matched per species by name, preferring 3rd-gen ONT long
reads (<bw_dir>/ONT/<species>.bw) and falling back to 2nd-gen BGI short reads
(<bw_dir>/BGI/<species>.bw); the order is configurable via --bw-priority.  The
genomic locus is taken from each FASTA header ('chrom:start-end(strand)');
'-' strand coverage is reversed to align with the sense sequence.
"""

from collections import Counter, deque
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime
import math
import os
from pathlib import Path
import re
from statistics import median

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from loguru import logger
import polars as pl
import pyBigWig
from tqdm import tqdm
import typer

from spider_silkome_module.config import INTERIM_DATA_DIR, PROCESSED_DATA_DIR, PROJ_ROOT

app = typer.Typer()

FASTA_FILE = "spidroin_full_length.fasta"
GFF_SOURCE = "kmer_intron_finder"
COVERAGE_FILE = "kmer_coverage.tsv"
PROTEIN_FILE = "spliced_proteins.faa"
FIGURE_DIRNAME = "figures"
# Base bigWig cache; per-platform subdirs (cache/bw/ONT, cache/bw/BGI) hold one
# <species>.bw each.  Coverage prefers 3rd-gen (ONT) long reads and falls back to
# 2nd-gen (BGI) short reads for species without an ONT bigWig (see --bw-priority).
DEFAULT_BW_DIR = PROJ_ROOT / "cache" / "bw"
DEFAULT_BW_PRIORITY = "ONT,BGI"

# ── calibrated defaults (see module docstring) ───────────────────────────────
DEFAULT_K = 12
DEFAULT_MIN_RUN = 70
DEFAULT_MIN_FLANK_COV = 3
DEFAULT_FLANK_WINDOW = 30
DEFAULT_SPLICE_SEARCH = 30
DEFAULT_FRAME_MOD = 3
DEFAULT_FRAME_REMAINDER = 1
DEFAULT_MIN_INTRON_LEN = 20


# ── data structure ───────────────────────────────────────────────────────────


@dataclass
class Intron:
    """One predicted intron on a (sense-strand) full-length spidroin sequence."""

    start: int              # 1-based inclusive, donor G position
    end: int                # 1-based inclusive, acceptor G position
    length: int
    donor: str              # donor dinucleotide, e.g. "GT"
    acceptor: str           # acceptor dinucleotide, e.g. "AG"
    candidate_start: int    # single-copy run that seeded this intron
    candidate_end: int
    left_flank_cov: float   # median k-mer copy number of the left flank
    right_flank_cov: float  # median k-mer copy number of the right flank
    span_diff: int          # |intron_length - candidate_span|


# ── k-mer profile ─────────────────────────────────────────────────────────────


def kmer_min_copy_profile(seq: str, k: int = DEFAULT_K) -> list[int]:
    """
    Per-base minimum copy number among all k-mers covering that base.

    Base j is covered by the k-mers starting at i in [j-k+1, j] (clipped to the
    valid range), so the profile is a sliding-window minimum over the per-k-mer
    copy counts, computed in O(n) with a monotonic deque.
    """
    n = len(seq)
    if n < k:
        return [1] * n
    counts = Counter(seq[i:i + k] for i in range(n - k + 1))
    kmer_copy = [counts[seq[i:i + k]] for i in range(n - k + 1)]  # length m = n-k+1
    m = n - k + 1

    cov = [1] * n
    dq: deque[int] = deque()  # k-mer indices, increasing kmer_copy
    for j in range(n):
        if j < m:  # k-mer index j enters the window exactly once
            while dq and kmer_copy[dq[-1]] >= kmer_copy[j]:
                dq.pop()
            dq.append(j)
        lo = max(0, j - k + 1)  # oldest k-mer index still covering base j
        while dq[0] < lo:
            dq.popleft()
        cov[j] = kmer_copy[dq[0]]
    return cov


def single_copy_runs(profile: list[int], min_run: int = DEFAULT_MIN_RUN) -> list[tuple[int, int]]:
    """Return 1-based inclusive (start, end) runs where copy == 1 and length >= min_run."""
    runs: list[tuple[int, int]] = []
    i, n = 0, len(profile)
    while i < n:
        if profile[i] == 1:
            j = i
            while j < n and profile[j] == 1:
                j += 1
            if j - i >= min_run:
                runs.append((i + 1, j))
            i = j
        else:
            i += 1
    return runs


def flank_copy(
    profile: list[int], cand_start: int, cand_end: int, window: int = DEFAULT_FLANK_WINDOW
) -> tuple[float, float]:
    """Median k-mer copy of the `window` bases immediately flanking a candidate run."""
    left = profile[max(0, cand_start - 1 - window):cand_start - 1]
    right = profile[cand_end:cand_end + window]
    left_cov = median(left) if left else 1.0
    right_cov = median(right) if right else 1.0
    return left_cov, right_cov


# ── boundary refinement ───────────────────────────────────────────────────────


def refine_boundaries(
    seq: str,
    cand_start: int,
    cand_end: int,
    *,
    splice_search: int = DEFAULT_SPLICE_SEARCH,
    frame_mod: int = DEFAULT_FRAME_MOD,
    frame_remainder: int = DEFAULT_FRAME_REMAINDER,
    min_intron_len: int = DEFAULT_MIN_INTRON_LEN,
) -> tuple[int, int, int, int] | None:
    """
    Search GT (donor) near cand_start and AG (acceptor) near cand_end within
    +/- splice_search nt.  Keep pairs whose intron length ≡ frame_remainder
    (mod frame_mod) and return the one closest to the candidate span as
    (donor_start, acceptor_end, length, span_diff), all 1-based; None if none.
    """
    n = len(seq)
    span = cand_end - cand_start + 1
    donors = [
        d for d in range(max(1, cand_start - splice_search), min(n - 1, cand_start + splice_search) + 1)
        if seq[d - 1:d + 1] == "GT"
    ]
    acceptors = [
        a for a in range(max(2, cand_end - splice_search), min(n, cand_end + splice_search) + 1)
        if seq[a - 2:a] == "AG"
    ]
    best: tuple[int, int, int, int] | None = None  # (span_diff, donor, acceptor, length)
    for d in donors:
        for a in acceptors:
            length = a - d + 1
            if length < min_intron_len:
                continue
            if frame_mod and length % frame_mod != frame_remainder:
                continue
            diff = abs(length - span)
            if best is None or diff < best[0]:
                best = (diff, d, a, length)
    if best is None:
        return None
    diff, donor, acceptor, length = best
    return donor, acceptor, length, diff


def _resolve_overlaps(introns: list[Intron]) -> list[Intron]:
    """Drop overlapping introns, keeping the one with the smallest span_diff."""
    kept: list[Intron] = []
    for it in sorted(introns, key=lambda x: x.start):
        if kept and it.start <= kept[-1].end:
            if it.span_diff < kept[-1].span_diff:
                kept[-1] = it
        else:
            kept.append(it)
    return kept


# ── main API ──────────────────────────────────────────────────────────────────


def introns_from_profile(
    seq: str,
    profile: list[int],
    *,
    min_run: int = DEFAULT_MIN_RUN,
    min_flank_cov: int = DEFAULT_MIN_FLANK_COV,
    flank_window: int = DEFAULT_FLANK_WINDOW,
    splice_search: int = DEFAULT_SPLICE_SEARCH,
    frame_mod: int = DEFAULT_FRAME_MOD,
    frame_remainder: int = DEFAULT_FRAME_REMAINDER,
    min_intron_len: int = DEFAULT_MIN_INTRON_LEN,
) -> list[Intron]:
    """Collect introns from a precomputed k-mer copy profile (callers reuse the profile)."""
    runs = single_copy_runs(profile, min_run)
    if len(runs) < 2:  # need a 5' NTD and a 3' CTD run to delimit the repeat body
        return []

    introns: list[Intron] = []
    for cand_start, cand_end in runs[1:-1]:  # drop NTD (first) and CTD (last)
        left_cov, right_cov = flank_copy(profile, cand_start, cand_end, flank_window)
        if min(left_cov, right_cov) < min_flank_cov:
            continue
        refined = refine_boundaries(
            seq, cand_start, cand_end,
            splice_search=splice_search, frame_mod=frame_mod,
            frame_remainder=frame_remainder, min_intron_len=min_intron_len,
        )
        if refined is None:
            continue
        donor, acceptor, length, span_diff = refined
        introns.append(Intron(
            start=donor,
            end=acceptor,
            length=length,
            donor=seq[donor - 1:donor + 1],
            acceptor=seq[acceptor - 2:acceptor],
            candidate_start=cand_start,
            candidate_end=cand_end,
            left_flank_cov=left_cov,
            right_flank_cov=right_cov,
            span_diff=span_diff,
        ))
    return _resolve_overlaps(introns)


def find_introns(
    seq: str,
    *,
    k: int = DEFAULT_K,
    min_run: int = DEFAULT_MIN_RUN,
    min_flank_cov: int = DEFAULT_MIN_FLANK_COV,
    flank_window: int = DEFAULT_FLANK_WINDOW,
    splice_search: int = DEFAULT_SPLICE_SEARCH,
    frame_mod: int = DEFAULT_FRAME_MOD,
    frame_remainder: int = DEFAULT_FRAME_REMAINDER,
    min_intron_len: int = DEFAULT_MIN_INTRON_LEN,
) -> list[Intron]:
    """Find introns in a sense-strand full-length spidroin sequence (see module docstring)."""
    seq = seq.upper()
    if len(seq) < k:
        return []
    profile = kmer_min_copy_profile(seq, k)
    return introns_from_profile(
        seq, profile,
        min_run=min_run, min_flank_cov=min_flank_cov, flank_window=flank_window,
        splice_search=splice_search, frame_mod=frame_mod,
        frame_remainder=frame_remainder, min_intron_len=min_intron_len,
    )


# ── protein translation ──────────────────────────────────────────────────────

_TYPE_RE = re.compile(r"type=(\S+)")


def parse_type(description: str) -> str | None:
    """Parse the spidroin type from a 'type=<...>' field in a FASTA header."""
    m = _TYPE_RE.search(description)
    return m.group(1) if m else None


def splice_out_introns(seq: str, introns: list[Intron]) -> str:
    """Drop every predicted intron (1-based inclusive start..end) from a sense sequence."""
    if not introns:
        return seq
    keep = [True] * len(seq)
    for it in introns:
        for pos in range(it.start - 1, it.end):  # 1-based inclusive -> 0-based slice
            keep[pos] = False
    return "".join(base for base, k in zip(seq, keep) if k)


def translate_cds(cds: str) -> tuple[str, bool]:
    """
    Translate a sense-strand CDS in frame 0 (the sequences start at the ATG).

    The CDS is trimmed to whole codons and translated without early stopping; the
    single terminal stop is dropped.  Returns (protein, has_premature_stop) where
    has_premature_stop is True when an internal '*' remains — for a spliced
    full-length spidroin that flags a mis-called intron rather than a real gene.
    """
    usable = cds[: len(cds) - (len(cds) % 3)]
    if len(usable) < 3:
        return "", False
    protein = str(Seq(usable).translate(to_stop=False)).rstrip("*")
    return protein, "*" in protein


# ── GFF output ────────────────────────────────────────────────────────────────


def write_introns_gff(
    introns_by_sid: dict[str, list[Intron]], output_path: Path, source: str = GFF_SOURCE
) -> None:
    """Write intron rows (1-based relative coords) aligned with fullprot_miniprot.gff."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as out:
        out.write("##gff-version 3\n")
        out.write(f"# Source: {source} (k-mer repeat-profile intron finder)\n")
        for sid in sorted(introns_by_sid):
            for i, it in enumerate(introns_by_sid[sid], 1):
                attrs = (
                    f"ID={sid}.intron{i};length={it.length};"
                    f"donor={it.donor};acceptor={it.acceptor};splice={it.donor}..{it.acceptor};"
                    f"candidate={it.candidate_start}-{it.candidate_end};"
                    f"left_flank_cov={it.left_flank_cov:.0f};"
                    f"right_flank_cov={it.right_flank_cov:.0f};span_diff={it.span_diff}"
                )
                # Full-length sequences are already sense-strand -> introns on '+'.
                out.write(
                    f"{sid}\t{source}\tintron\t{it.start}\t{it.end}\t.\t+\t.\t{attrs}\n"
                )


# ── RNA-seq bigWig coverage ───────────────────────────────────────────────────

_LOCUS_RE = re.compile(r"(\S+):(\d+)-(\d+)\(([+-])\)")


def parse_locus(description: str) -> tuple[str, int, int, str] | None:
    """Parse 'chrom:start-end(strand)' from a spidroin_full_length FASTA header."""
    m = _LOCUS_RE.search(description)
    if m is None:
        return None
    return m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)


def rna_coverage(bw, description: str, n: int) -> list[float] | None:
    """
    Per-base RNA-seq coverage aligned to a sense-strand sequence of length n.

    Reads the genomic interval from the open bigWig, fills missing values with
    0, and reverses it for '-' strand loci so the array maps to relative
    positions 1..n.  Returns None if the locus cannot be mapped (header
    unparseable, chromosome absent, or length mismatch).
    """
    locus = parse_locus(description)
    if locus is None:
        return None
    chrom, start, end, strand = locus
    if chrom not in bw.chroms() or (end - start + 1) != n:
        return None
    try:
        raw = bw.values(chrom, start - 1, end)
    except (RuntimeError, OverflowError):
        return None
    vals = [0.0 if (v is None or math.isnan(v)) else float(v) for v in raw]
    if strand == "-":
        vals.reverse()
    return vals


# ── per-base coverage output ──────────────────────────────────────────────────


def coverage_block(
    seq_id: str, seq: str, profile: list[int], rna: list[float] | None = None
) -> str:
    """
    Render one sequence's per-base coverage as a tab-separated text block.

    Columns: position (1-based), base, k-mer coverage (= min k-mer copy number),
    and — when an RNA-seq bigWig is provided — RNA-seq coverage (4th column;
    'NA' where the base could not be mapped).  A leading '# <seq_id>' header
    line separates sequences within a species file.
    """
    lines = [f"# {seq_id}"]
    if rna is None:
        lines.extend(f"{i}\t{base}\t{cov}" for i, (base, cov) in enumerate(zip(seq, profile), 1))
    else:
        for i, (base, cov, r) in enumerate(zip(seq, profile, rna), 1):
            rcell = "NA" if (r is None or math.isnan(r)) else f"{r:g}"
            lines.append(f"{i}\t{base}\t{cov}\t{rcell}")
    return "\n".join(lines) + "\n"


def plot_coverage(
    seq_id: str,
    profile: list[int],
    introns: list[Intron],
    output_path: Path,
    *,
    ntd_ctd: tuple[tuple[int, int], tuple[int, int]] | None = None,
    rna: list[float] | None = None,
) -> None:
    """
    Plot per-base coverage for one sequence: x = position (bp), y = coverage.

    The left (log) y-axis shows k-mer coverage as a filled area, so single-copy
    "valleys" (NTD/CTD and introns) and the high-copy repeat body are both
    visible.  Predicted introns are shaded red; the terminal NTD/CTD single-copy
    domains are shaded grey.  When an RNA-seq bigWig is provided, its coverage is
    drawn as a green line on a secondary (linear) y-axis — RNA-seq drops to zero
    across genuine introns, independently corroborating the prediction.
    """
    import matplotlib

    matplotlib.use("Agg")  # headless backend, safe inside worker processes
    import matplotlib.pyplot as plt

    n = len(profile)
    fig, ax = plt.subplots(figsize=(12, 3.2))
    if n:
        ax.fill_between(range(1, n + 1), profile, 0.8, color="#4C72B0", lw=0, alpha=0.55)
    ax.set_yscale("log")
    ax.set_ylim(0.8, (max(profile) * 1.4) if profile else 10)
    ax.set_xlim(1, max(n, 1))
    ax.set_xlabel("Position (bp)")
    ax.set_ylabel("k-mer coverage (min copy)", color="#1f3b73")
    ax.set_title(f"{seq_id}  (length {n} bp)")
    ax.spines[["top"]].set_visible(False)

    top = ax.get_ylim()[1]
    if ntd_ctd is not None:
        for (start, end), label in zip(ntd_ctd, ("NTD", "CTD")):
            ax.axvspan(start, end, color="#9E9E9E", alpha=0.18)
            ax.text((start + end) / 2, top, label, ha="center", va="top",
                    fontsize=8, color="#555555")
    for i, it in enumerate(introns, 1):
        ax.axvspan(it.start, it.end, color="#D55E00", alpha=0.30)
        ax.text((it.start + it.end) / 2, 1.0, f"intron{i}\n{it.donor}..{it.acceptor}",
                ha="center", va="bottom", fontsize=7, color="#A23B00")

    if rna is not None and any(not math.isnan(r) for r in rna):
        ax2 = ax.twinx()
        ax2.plot(range(1, n + 1), rna, color="#2CA02C", lw=0.8, alpha=0.9)
        ax2.set_ylabel("RNA-seq coverage", color="#1b7a1b")
        ax2.set_ylim(bottom=0)
        ax2.spines[["top"]].set_visible(False)

    fig.savefig(output_path, dpi=120, bbox_inches="tight")
    plt.close(fig)


# ── per-species worker (top-level for ProcessPoolExecutor pickling) ───────────


def process_species_file(
    fasta_path: Path,
    output_gff: Path,
    params: dict,
    force: bool,
    emit_coverage: bool = True,
    emit_plots: bool = True,
    emit_proteins: bool = True,
    bw_dirs: list[Path] | None = None,
) -> dict:
    """Find introns for every sequence in one species fasta; write GFF (+ coverage + proteins + plots)."""
    species = fasta_path.parent.name
    out_dir = output_gff.parent
    cov_path = out_dir / COVERAGE_FILE
    prot_path = out_dir / PROTEIN_FILE
    fig_dir = out_dir / FIGURE_DIRNAME
    expected = [output_gff]
    if emit_coverage:
        expected.append(cov_path)
    if emit_proteins:
        expected.append(prot_path)
    if emit_plots:
        expected.append(fig_dir)
    if all(p.exists() for p in expected) and not force:
        return {"species": species, "n_seqs": 0, "n_introns": 0, "n_rna": 0,
                "rna_platform": None, "skipped": True, "rows": [], "protein_rows": []}

    k = params["k"]
    min_run = params["min_run"]
    intron_params = {key: val for key, val in params.items() if key != "k"}
    out_dir.mkdir(parents=True, exist_ok=True)
    if emit_plots:
        fig_dir.mkdir(parents=True, exist_ok=True)

    # Open the species RNA-seq bigWig, preferring the first platform dir that has
    # one (ONT 3rd-gen) and falling back to later ones (BGI 2nd-gen).
    bw = None
    rna_platform = None
    for d in bw_dirs or []:
        cand = d / f"{species}.bw"
        if not cand.exists():
            continue
        try:
            bw = pyBigWig.open(str(cand))
            rna_platform = d.name
        except RuntimeError:
            bw = None
        if bw is not None:
            break

    introns_by_sid: dict[str, list[Intron]] = {}
    rows: list[dict] = []
    protein_records: list[SeqRecord] = []
    protein_rows: list[dict] = []
    n_seqs = n_rna = 0
    cov_fh = open(cov_path, "w") if emit_coverage else None
    try:
        if cov_fh is not None:
            cols = "position(1-based), base, kmer_coverage" + (
                ", rna_coverage" if bw is not None else "")
            src = f"; rna source: {rna_platform}" if rna_platform else ""
            cov_fh.write(f"# per-base coverage; columns: {cols}{src}\n")
        for rec in SeqIO.parse(str(fasta_path), "fasta"):
            n_seqs += 1
            seq = str(rec.seq).upper()
            # Profile is computed once and reused for coverage, intron calling and plotting.
            profile = kmer_min_copy_profile(seq, k) if len(seq) >= k else [1] * len(seq)
            # RNA-seq coverage aligned to the sense sequence (None if no bigWig);
            # NaN placeholder when the species has a bigWig but this locus is unmappable.
            rna = None
            if bw is not None:
                rna = rna_coverage(bw, rec.description, len(seq))
                if rna is None:
                    rna = [math.nan] * len(seq)
                else:
                    n_rna += 1
            if cov_fh is not None:
                cov_fh.write(coverage_block(rec.id, seq, profile, rna))
            introns = introns_from_profile(seq, profile, **intron_params)
            if emit_plots:
                runs = single_copy_runs(profile, min_run)
                ntd_ctd = (runs[0], runs[-1]) if len(runs) >= 2 else None
                plot_coverage(rec.id, profile, introns, fig_dir / f"{rec.id}.png",
                              ntd_ctd=ntd_ctd, rna=rna)
            if introns:
                introns_by_sid[rec.id] = introns
                for i, it in enumerate(introns, 1):
                    rows.append({
                        "species": species,
                        "spidroin_id": rec.id,
                        "intron_index": i,
                        "start": it.start,
                        "end": it.end,
                        "length": it.length,
                        "donor": it.donor,
                        "acceptor": it.acceptor,
                        "candidate_start": it.candidate_start,
                        "candidate_end": it.candidate_end,
                        "left_flank_cov": round(it.left_flank_cov, 1),
                        "right_flank_cov": round(it.right_flank_cov, 1),
                        "span_diff": it.span_diff,
                    })
            if emit_proteins:
                # Splice out every predicted intron and translate the sense CDS;
                # an internal stop flags a mis-called intron (premature_stop=yes).
                protein, premature = translate_cds(splice_out_introns(seq, introns))
                locus_field = rec.description.split(None, 1)
                desc = (
                    f"{locus_field[1] if len(locus_field) > 1 else ''} "
                    f"introns={len(introns)} "
                    f"premature_stop={'yes' if premature else 'no'} "
                    f"length={len(protein)}"
                ).strip()
                protein_records.append(SeqRecord(Seq(protein), id=rec.id, description=desc))
                protein_rows.append({
                    "species": species,
                    "spidroin_id": rec.id,
                    "spidroin_type": parse_type(rec.description) or "",
                    "n_introns": len(introns),
                    "protein_length": len(protein),
                    "has_premature_stop": premature,
                })
    finally:
        if cov_fh is not None:
            cov_fh.close()
        if bw is not None:
            bw.close()

    write_introns_gff(introns_by_sid, output_gff)
    if emit_proteins:
        SeqIO.write(protein_records, str(prot_path), "fasta")
    return {
        "species": species,
        "n_seqs": n_seqs,
        "n_introns": len(rows),
        "n_rna": n_rna,
        "rna_platform": rna_platform,
        "skipped": False,
        "rows": rows,
        "protein_rows": protein_rows,
    }


# ── CLI ────────────────────────────────────────────────────────────────────────


@app.command()
def main(
    typing_dir: Path = PROCESSED_DATA_DIR / "typing_results",
    output_dir: Path = INTERIM_DATA_DIR / f"intron_kmer_{datetime.now():%Y%m%d_%H%M%S}",
    k: int = DEFAULT_K,
    min_run: int = DEFAULT_MIN_RUN,
    min_flank_cov: int = DEFAULT_MIN_FLANK_COV,
    flank_window: int = DEFAULT_FLANK_WINDOW,
    splice_search: int = DEFAULT_SPLICE_SEARCH,
    frame_mod: int = DEFAULT_FRAME_MOD,
    frame_remainder: int = DEFAULT_FRAME_REMAINDER,
    min_intron_len: int = DEFAULT_MIN_INTRON_LEN,
    fasta_name: str = FASTA_FILE,
    species: str | None = None,
    workers: int = 75,
    force: bool = False,
    emit_coverage: bool = True,
    emit_plots: bool = True,
    emit_proteins: bool = True,
    bw_dir: Path = DEFAULT_BW_DIR,
    bw_priority: str = DEFAULT_BW_PRIORITY,
):
    """Find k-mer-based introns across all species under typing_dir (parallel, resumable)."""
    # Pin BLAS/OMP threads so per-worker libraries don't oversubscribe cores.
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ.setdefault(var, "1")

    # Resolve RNA-seq bigWig platform dirs in priority order (ONT first, BGI fallback).
    bw_dirs = [bw_dir / name for name in bw_priority.split(",") if (bw_dir / name).is_dir()]
    if not bw_dirs:
        logger.warning(f"no bigWig platform dirs under {bw_dir} ({bw_priority}); RNA-seq coverage disabled")
        bw_dirs = None
    else:
        logger.info(f"RNA-seq bigWig priority: {[d.name for d in bw_dirs]}")

    species_dirs = sorted(d for d in typing_dir.iterdir() if d.is_dir())
    if species is not None:
        species_dirs = [d for d in species_dirs if d.name == species]
    tasks = [
        (d / fasta_name, output_dir / d.name / "kmer_introns.gff")
        for d in species_dirs
        if (d / fasta_name).exists()
    ]
    if not tasks:
        logger.error(f"No '{fasta_name}' found under {typing_dir} matching filter")
        raise typer.Exit(1)

    params = dict(
        k=k, min_run=min_run, min_flank_cov=min_flank_cov, flank_window=flank_window,
        splice_search=splice_search, frame_mod=frame_mod, frame_remainder=frame_remainder,
        min_intron_len=min_intron_len,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Scanning {len(tasks)} species | min_run={min_run} min_flank_cov={min_flank_cov}")

    all_rows: list[dict] = []
    all_protein_rows: list[dict] = []
    n_introns = n_seqs = n_skipped = n_rna = 0
    platform_counts: Counter[str] = Counter()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = {
            ex.submit(process_species_file, fa, gff, params, force,
                      emit_coverage, emit_plots, emit_proteins, bw_dirs): fa
            for fa, gff in tasks
        }
        for fut in tqdm(as_completed(futures), total=len(futures), desc="Finding introns"):
            res = fut.result()
            if res["skipped"]:
                n_skipped += 1
                continue
            all_rows.extend(res["rows"])
            all_protein_rows.extend(res["protein_rows"])
            n_introns += res["n_introns"]
            n_seqs += res["n_seqs"]
            n_rna += res["n_rna"]
            if res["rna_platform"]:
                platform_counts[res["rna_platform"]] += 1

    if all_rows:
        pl.DataFrame(all_rows).write_csv(output_dir / "intron_summary.tsv", separator="\t")
    if all_protein_rows:
        pl.DataFrame(all_protein_rows).write_csv(output_dir / "protein_summary.tsv", separator="\t")
    n_premature = sum(1 for r in all_protein_rows if r["has_premature_stop"])
    platforms = ", ".join(f"{name}:{platform_counts[name]}" for name in sorted(platform_counts)) or "none"
    logger.success(
        f"Done: {n_introns} introns over {n_seqs} sequences "
        f"({n_rna} with RNA-seq [{platforms}], {n_skipped} species skipped); "
        f"{len(all_protein_rows)} proteins, {n_premature} with premature stop -> {output_dir}"
    )


if __name__ == "__main__":
    app()
