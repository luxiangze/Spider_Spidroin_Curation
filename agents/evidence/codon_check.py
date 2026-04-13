"""
codon_check.py — verify and correct start/stop codon boundaries for spidroin candidates.

For each locus with an NTD and/or CTD hit, extracts a genomic window flanking the
domain boundary, converts to mRNA orientation, and locates the nearest ATG (start)
or stop codon.  When a codon is found the locus coordinates are corrected to include
the actual codon.

Algorithm follows extract_spidroin_sequences() in
references/2025_Schoneberg_data/analyse_spidroins.py:
  1. Start codon — extract [ntd_start − extension, ntd_end] ('+') or
     [ntd_start, ntd_end + extension] ('−', reverse-complemented) and find the
     first ATG reading 5′→3′.
  2. Stop codon  — extract [ctd_end − overlap, ctd_end + extension] ('+') or
     [ctd_start − extension, ctd_start + overlap] ('−', reverse-complemented), reverse
     the mRNA sequence, and find the first reversed-stop pattern
     (AAT / GAT / AGT = reversed TAA / TAG / TGA), which corresponds to the
     rightmost / last stop codon in the mRNA window.
  3. Update locus.start / locus.end and the NTD/CTD sub-coordinates accordingly.

Genome FASTA path conventions (read-only):
  Old annotation: data/raw/01.ref_gff/<species_id>/*.genome.fa.gz
  New annotation: data/raw/spider_anno2/<NNN>spd/genome.softmasked.fa.gz

Coordinates throughout are 1-based closed intervals (GFF3 convention).
pysam fetch uses 0-based half-open intervals; conversions are documented inline.
"""

from itertools import groupby
from pathlib import Path
import re

from loguru import logger
import pysam

# ── constants ────────────────────────────────────────────────────────────────

# Extension added beyond the domain boundary when searching for codons (bp).
_EXTENSION_BP = 150

# Bases included on the "inside" of the CTD boundary so that a stop codon
# embedded at the very tail of the HMM alignment is not missed.
_BOUNDARY_OVERLAP = 50

# Reversed stop-codon triplets.  Reversing a sequence and searching for the
# first of these patterns locates the *last* (rightmost in mRNA 5′→3′) stop
# codon in the window — equivalent to the inv_ctd trick in analyse_spidroins.py.
_STOP_REVERSED_RE = re.compile(r"AAT|GAT|AGT")  # reversed TAA | TAG | TGA

_REVCOMP = str.maketrans("ACGTacgt", "TGCAtgca")


# ── helpers ───────────────────────────────────────────────────────────────────

def _revcomp(seq: str) -> str:
    """Return the reverse complement of a DNA string."""
    return seq.translate(_REVCOMP)[::-1]


def _find_genome_fasta(
    fasta_root: Path | None,
    species_id: str,
    fasta_root_new: Path | None = None,
) -> Path | None:
    """
    Locate a bgzip-compressed FASTA file for species_id.
    Requires a sibling .fai index (created by `samtools faidx`).

    Search order:
      1. Old annotation layout: {fasta_root}/{species_id}/*.genome.fa.gz
      2. New annotation layout: {fasta_root_new}/{NNN}spd/genome.softmasked.fa.gz
         where NNN is the numeric prefix of species_id (e.g. "001" from
         "001.Allagelena_difficilis").

    Returns None if no suitable file is found.
    """
    if fasta_root is not None:
        sp_dir = fasta_root / species_id
        if sp_dir.is_dir():
            for pattern in ("*.genome.fa.gz", "*.genome.fna.gz", "*.fa.gz", "*.fna.gz"):
                for candidate in sp_dir.glob(pattern):
                    if Path(str(candidate) + ".fai").exists():
                        return candidate

    if fasta_root_new is not None:
        nnn = species_id.split(".")[0]
        new_sp_dir = fasta_root_new / f"{nnn}spd"
        if new_sp_dir.is_dir():
            for pattern in ("genome.softmasked.fa.gz", "genome.fa.gz", "*.fa.gz", "*.fna.gz"):
                for candidate in new_sp_dir.glob(pattern):
                    if Path(str(candidate) + ".fai").exists():
                        return candidate

    logger.debug(
        f"No bgzip+fai-indexed genome FASTA found for {species_id}. "
        "Run `samtools faidx` on the bgzip-compressed FASTA to generate the index."
    )
    return None


# ── codon search ─────────────────────────────────────────────────────────────

def _find_start_codon(
    fasta: pysam.FastaFile,
    chrom: str,
    ntd_start: int,
    ntd_end: int,
    strand: str,
    extension_bp: int = _EXTENSION_BP,
) -> tuple[bool | None, int | None]:
    """
    Find the ATG nearest to the NTD 5′ boundary (inner-to-outer search).

    Search window (1-based, inclusive):
      '+' strand: [ntd_start − extension_bp, ntd_start]
                  → rfind("ATG") from right to left: returns the RIGHTMOST
                    ATG, i.e. the one closest to the nHMMER boundary.
      '−' strand: [ntd_end, ntd_end + extension_bp]
                  → reverse-complement → rfind("ATG") from right to left
                    in the revcomp sequence, which returns the ATG closest
                    to ntd_end from the 5′-upstream direction.

    Short-circuit: if the boundary position itself already carries ATG
    (or its revcomp on '−'), that coordinate is returned immediately
    without further searching.

    Returns:
        (True,  new_coord) — ATG found.
            '+': 1-based genomic position of the A in ATG.
            '−': 1-based genomic position of the *highest* base of the ATG
                 triplet (= 5′ mRNA boundary, to be stored as ntd_end).
        (False, None)  — no ATG within the window.
        (None,  None)  — chrom absent from FASTA or retrieval error.
    """
    if chrom not in fasta.references:
        return None, None
    chrom_len = fasta.get_reference_length(chrom)
    try:
        if strand == "+":
            # Short-circuit: check if ntd_start itself is ATG.
            if ntd_start + 2 <= chrom_len:
                if fasta.fetch(chrom, ntd_start - 1, ntd_start + 2).upper() == "ATG":
                    return True, ntd_start
            # Search upstream of ntd_start.  rfind returns the RIGHTMOST ATG
            # = the one closest to the nHMMER boundary (inner-to-outer order).
            # ntd_start is already shifted −3 from the raw nHMMER end (via the CTD
            # symmetric logic). Extend 5 bp into the domain so an ATG that starts at
            # or 2 bp inside the original nHMMER boundary is also found.
            fetch_start = max(0, ntd_start - 1 - extension_bp)
            fetch_end   = min(chrom_len, ntd_end)  # same window as before
            seq = fasta.fetch(chrom, fetch_start, fetch_end).upper()
            inner_idx = ntd_start - 1 - fetch_start  # 0-based index of ntd_start in seq
            idx = seq[:inner_idx + 5].rfind("ATG")
            if idx == -1:
                return False, None
            return True, fetch_start + idx + 1
        else:
            # '−' strand: 5′ mRNA end is at the HIGHER genomic coordinate (ntd_end).
            # Short-circuit: check if ntd_end already points to ATG (mRNA direction).
            if ntd_end >= 3:
                if _revcomp(fasta.fetch(chrom, ntd_end - 3, ntd_end).upper()) == "ATG":
                    return True, ntd_end
            # Fetch the full window (NTD domain + upstream extension) to provide
            # sequence context for codons that straddle ntd_end.
            # After revcomp: seq[i] ↔ 1-based (fetch_end − i).
            # seq[0] = outermost (highest genomic = ntd_end + extension_bp).
            # seq[inner_idx] ≈ ntd_end (inner boundary).
            # rfind in seq[0:inner_idx] returns the ATG with the LARGEST index
            # j < inner_idx, i.e. the ATG closest to ntd_end from upstream.
            fetch_start = max(0, ntd_start - 1)   # same window as before
            fetch_end   = min(chrom_len, ntd_end + extension_bp)
            seq = _revcomp(fasta.fetch(chrom, fetch_start, fetch_end).upper())
            # inner_idx: 0-based revcomp position corresponding to ntd_end.
            # Positions 0..inner_idx−1 in seq are upstream of ntd_end (outer region).
            inner_idx = fetch_end - ntd_end        # = extension_bp (if no clipping)
            # ntd_end is already shifted +3 into the domain from the raw nHMMER end.
            # Extend 5 bp into the domain so an ATG at or 2 bp inside the original
            # nHMMER boundary is also found.
            idx = seq[:inner_idx + 5].rfind("ATG")
            if idx == -1:
                return False, None
            # ATG first base (A) in mRNA ↔ 1-based genomic (fetch_end − idx)
            return True, fetch_end - idx
    except Exception:
        return None, None


def _find_stop_codon(
    fasta: pysam.FastaFile,
    chrom: str,
    ctd_start: int,
    ctd_end: int,
    strand: str,
    extension_bp: int = _EXTENSION_BP,
    boundary_overlap: int = _BOUNDARY_OVERLAP,
) -> tuple[bool | None, int | None]:
    """
    Find the stop codon nearest to the CTD 3′ boundary (inner-to-outer search).

    Search window (1-based, inclusive):
      '+' strand: [ctd_end − boundary_overlap, ctd_end + extension_bp]
                  → forward scan from the inner boundary (ctd_end) outward:
                    returns the FIRST (leftmost) stop codon encountered.
      '−' strand: [ctd_start − extension_bp, ctd_start + boundary_overlap]
                  → reverse-complement → forward scan in revcomp from the
                    inner boundary outward: returns the FIRST stop codon
                    encountered (= lowest genomic coord = closest to ctd_start).

    The forward scan (inner → outer) ensures the stop codon closest to the
    nHMMER boundary is selected, rather than one farther downstream.

    Short-circuit: if the boundary position itself carries a stop codon,
    that coordinate is returned immediately without further searching.

    Returns:
        (True,  new_coord) — stop codon found.
            '+': 1-based genomic position of the LAST base of the stop codon
                 (= highest coord; to be stored as ctd_end).
            '−': 1-based genomic position of the FIRST base of the stop codon
                 (= lowest coord; to be stored as ctd_start).
        (False, None)  — no stop codon within the window.
        (None,  None)  — chrom absent from FASTA or retrieval error.
    """
    if chrom not in fasta.references:
        return None, None
    chrom_len = fasta.get_reference_length(chrom)
    _STOP_CODONS = {"TAA", "TAG", "TGA"}
    try:
        if strand == "+":
            # Short-circuit: check if ctd_end is already the last base of a stop codon.
            if ctd_end >= 3:
                codon = fasta.fetch(chrom, ctd_end - 3, ctd_end).upper()
                if codon in _STOP_CODONS:
                    return True, ctd_end
            fetch_start = max(0, ctd_end - 1 - boundary_overlap)
            fetch_end   = min(chrom_len, ctd_end + extension_bp)
            seq = fasta.fetch(chrom, fetch_start, fetch_end).upper()
            # inner_idx: 0-based position of ctd_end in seq (= boundary_overlap if no clipping).
            # ctd_end is already extended by +3 beyond the raw nHMMER end.  To make sure a
            # stop codon that starts right at the original nHMMER end (3 bp inside ctd_end)
            # is not missed, start scanning 5 positions before inner_idx (3 for the +3
            # extension offset + 2 extra for codons straddling the original boundary).
            inner_idx = ctd_end - 1 - fetch_start  # = boundary_overlap (if no clipping)
            for i in range(max(0, inner_idx - 5), len(seq) - 2):
                if seq[i : i + 3] in _STOP_CODONS:
                    # 1-based last base of stop codon = fetch_start + i + 3
                    return True, fetch_start + i + 3
            return False, None
        else:
            # '−' strand: 3′ mRNA end is at the LOWER genomic coordinate (ctd_start).
            # Short-circuit: check if ctd_start is the first base of a stop codon (mRNA direction).
            if ctd_start + 2 <= chrom_len:
                codon = _revcomp(fasta.fetch(chrom, ctd_start - 1, ctd_start + 2).upper())
                if codon in _STOP_CODONS:
                    return True, ctd_start
            # After revcomp of [fetch_start, fetch_end): seq[i] ↔ 1-based (fetch_end − i).
            # seq[0] = innermost (genomic ctd_start + boundary_overlap, inside CTD boundary).
            # seq[inner_idx] = ctd_start (inner boundary, after −3 extension).
            # seq[-1] = outermost (genomic ctd_start − extension_bp).
            # ctd_start is already extended by −3 beyond the raw nHMMER start.  Scan from
            # inner_idx − 5 outward so stop codons at the original nHMMER boundary are found.
            fetch_start = max(0, ctd_start - 1 - extension_bp)
            fetch_end   = min(chrom_len, ctd_start + boundary_overlap)
            seq = _revcomp(fasta.fetch(chrom, fetch_start, fetch_end).upper())
            inner_idx = fetch_end - ctd_start  # = boundary_overlap (if no clipping)
            for i in range(max(0, inner_idx - 5), len(seq) - 2):
                if seq[i : i + 3] in _STOP_CODONS:
                    # 1-based lowest genomic coord of stop codon = fetch_end − i − 2
                    return True, fetch_end - i - 2
            return False, None
    except Exception:
        return None, None


# ── main entry point ──────────────────────────────────────────────────────────

def annotate_loci(
    loci: list[dict],
    fasta_root: Path | None,
    fasta_root_new: Path | None = None,
) -> list[dict]:
    """
    Locate start/stop codons near NTD/CTD boundaries and update locus coordinates.

    For each locus:
      • If ntd_start/ntd_end are set, calls _find_start_codon and — when a codon
        is found — updates ntd_start ('+') or ntd_end ('−') and the locus boundary.
      • If ctd_start/ctd_end are set, calls _find_stop_codon and — when a codon
        is found — updates ctd_end ('+') or ctd_start ('−') and the locus boundary.

    Sets per-locus keys:
      ntd_start_codon (bool | None) — True if ATG was found near the NTD boundary
      ctd_stop_codon  (bool | None) — True if a stop codon was found near the CTD boundary

    Appends to locus["notes"]:
      "boundary:start+stop"  both codons confirmed
      "boundary:no_start"    ATG not found within the search window
      "boundary:no_stop"     stop codon not found within the search window
    """
    for locus in loci:
        locus.setdefault("ntd_start_codon", None)
        locus.setdefault("ctd_stop_codon", None)

    # Sort by species so groupby can batch FASTA opens; stable sort preserves
    # intra-species locus order.
    loci.sort(key=lambda loc: loc["species"])

    for species_id, species_iter in groupby(loci, key=lambda loc: loc["species"]):
        species_loci = list(species_iter)

        fasta_path = _find_genome_fasta(fasta_root, species_id, fasta_root_new)
        if fasta_path is None:
            continue

        try:
            fasta = pysam.FastaFile(str(fasta_path))
        except Exception as exc:
            logger.warning(f"Failed to open FASTA {fasta_path}: {exc}")
            continue

        with fasta:
            logger.debug(
                f"Opened bgzip FASTA for {species_id}: {fasta_path} "
                f"({len(fasta.references)} sequences)"
            )
            for locus in species_loci:
                seqid  = locus["seqid"]
                strand = locus.get("strand", "+")

                # ── start codon ───────────────────────────────────────────
                ntd_start_codon: bool | None = None
                if locus.get("ntd_start") is not None:
                    found, new_coord = _find_start_codon(
                        fasta, seqid,
                        locus["ntd_start"], locus["ntd_end"],
                        strand,
                    )
                    ntd_start_codon = found
                    if found and new_coord is not None:
                        if strand == "+":
                            locus["ntd_start"] = new_coord
                            locus["start"] = min(locus["start"], new_coord)
                        else:
                            locus["ntd_end"] = new_coord
                            locus["end"] = max(locus["end"], new_coord)

                # ── stop codon ────────────────────────────────────────────
                ctd_stop_codon: bool | None = None
                if locus.get("ctd_start") is not None:
                    found, new_coord = _find_stop_codon(
                        fasta, seqid,
                        locus["ctd_start"], locus["ctd_end"],
                        strand,
                    )
                    ctd_stop_codon = found
                    if found and new_coord is not None:
                        if strand == "+":
                            locus["ctd_end"] = new_coord
                            locus["end"] = max(locus["end"], new_coord)
                        else:
                            locus["ctd_start"] = new_coord
                            locus["start"] = min(locus["start"], new_coord)

                locus["ntd_start_codon"] = ntd_start_codon
                locus["ctd_stop_codon"]  = ctd_stop_codon

                # ── boundary notes ────────────────────────────────────────
                existing = locus.get("notes", "") or ""
                boundary_notes: list[str] = []

                if ntd_start_codon is False:
                    boundary_notes.append("boundary:no_start")
                if ctd_stop_codon is False:
                    boundary_notes.append("boundary:no_stop")
                if ntd_start_codon is True and ctd_stop_codon is True and not boundary_notes:
                    boundary_notes.append("boundary:start+stop")

                if boundary_notes:
                    sep = "; " if existing else ""
                    locus["notes"] = existing + sep + "; ".join(boundary_notes)

    return loci
