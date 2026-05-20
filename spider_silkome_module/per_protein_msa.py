"""
per_protein_msa.py — produce per-protein domain-level multiple sequence
alignments for confirmed spidroin predictions.

For each predicted spidroin, the NTD and CTD regions are sliced from the
predicted protein (using hmmscan coordinates from protein_validation), then
aligned with same-type and same-family-group reference panels using MAFFT
L-INS-i, then trimmed with trimAl.

Outputs (per spidroin):
  msa/<species>/<spidroin_id>/
    ntd_strict.aln    NTD vs same-type reference panel
    ntd_loose.aln     NTD vs same-family-group reference panel
    ctd_strict.aln
    ctd_loose.aln
    overview.aln      full predicted protein vs same-type panel (MAFFT --auto)

Reference panels are read from data/interim/spidroin_proteins/<TYPE>_<NTD|CTD>.fa.
"""

from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from loguru import logger
import polars as pl
from tqdm import tqdm
import typer

from spider_silkome_module.config import INTERIM_DATA_DIR, PROCESSED_DATA_DIR
from spider_silkome_module.protein_validation import FAMILY_GROUPS
from spider_silkome_module.utils.run_cmd import run_cmd

app = typer.Typer()

DEFAULT_TASK_NAME = "protein_confirmation"
DEFAULT_REF_DIR = INTERIM_DATA_DIR / "spidroin_proteins"

# A panel needs at least this many references for a useful MSA.
_MIN_PANEL_SIZE = 3
# Cap reference panel sizes so MAFFT runs in a reasonable time.
_MAX_PANEL_SIZE = 80


# ── reference panel handling ─────────────────────────────────────────────────


def _load_reference_panel(ref_dir: Path, type_name: str, kind: str) -> list[SeqRecord]:
    """Load <type>_{NTD|CTD}.fa from ref_dir; return [] if missing."""
    fa = ref_dir / f"{type_name}_{kind}.fa"
    if not fa.exists():
        return []
    return list(SeqIO.parse(fa, "fasta"))


def build_loose_panel(
    ref_dir: Path,
    target_type: str,
    kind: str,
) -> list[SeqRecord]:
    """Concatenate all types belonging to the same FAMILY_GROUPS bucket."""
    target_group = FAMILY_GROUPS.get(target_type)
    if target_group is None:
        return _load_reference_panel(ref_dir, target_type, kind)
    members = [t for t, g in FAMILY_GROUPS.items() if g == target_group]
    seen: set[str] = set()
    panel: list[SeqRecord] = []
    for t in members:
        for rec in _load_reference_panel(ref_dir, t, kind):
            if rec.id in seen:
                continue
            seen.add(rec.id)
            panel.append(rec)
    return panel


def _cap_panel(panel: list[SeqRecord], n: int = _MAX_PANEL_SIZE) -> list[SeqRecord]:
    return panel if len(panel) <= n else panel[:n]


# ── slicing predicted proteins ───────────────────────────────────────────────


def _load_predicted(pred_fa: Path) -> dict[str, str]:
    return {rec.id: str(rec.seq).upper().rstrip("*") for rec in SeqIO.parse(pred_fa, "fasta")}


def slice_domain(seq: str, start: int | None, end: int | None,
                 fallback: str = "ntd") -> str:
    """
    Slice [start, end] 1-based inclusive.  When coords are missing, fallback to
    the first 100 aa (NTD) or the last 100 aa (CTD).
    """
    if seq == "":
        return ""
    if start is not None and end is not None and 1 <= start <= end <= len(seq):
        return seq[start - 1:end]
    if fallback == "ntd":
        return seq[:min(120, len(seq))]
    return seq[-min(120, len(seq)):]


# ── MSA runner ───────────────────────────────────────────────────────────────


def run_mafft_linsi(input_fa: Path, output_aln: Path, threads: int,
                    force: bool = False) -> None:
    cmd = (f"mafft --quiet --localpair --maxiterate 100 --thread {threads} "
           f"{input_fa} > {output_aln}")
    run_cmd(cmd, [output_aln], force=force)


def run_mafft_auto(input_fa: Path, output_aln: Path, threads: int,
                   force: bool = False) -> None:
    cmd = f"mafft --quiet --retree 2 --thread {threads} {input_fa} > {output_aln}"
    run_cmd(cmd, [output_aln], force=force)


def run_trimal(input_aln: Path, output_aln: Path, force: bool = False) -> None:
    cmd = f"trimal -in {input_aln} -out {output_aln} -gappyout"
    run_cmd(cmd, [output_aln], force=force)


# ── per-spidroin orchestration ───────────────────────────────────────────────


def write_panel_fasta(records: list[SeqRecord], output_fa: Path) -> None:
    output_fa.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, output_fa, "fasta")


def make_msa(
    seq_id: str,
    seq_str: str,
    panel: list[SeqRecord],
    workdir: Path,
    tag: str,
    threads: int,
    use_linsi: bool,
    force: bool,
) -> Path | None:
    """
    Combine query + panel into a fasta and run MAFFT.  Returns the trimmed
    .aln path or None if the panel was too small.
    """
    if len(panel) < _MIN_PANEL_SIZE:
        return None
    panel = _cap_panel(panel)

    combined = workdir / f"{tag}.input.fa"
    raw_aln = workdir / f"{tag}.raw.aln"
    final_aln = workdir / f"{tag}.aln"
    if final_aln.exists() and not force:
        return final_aln

    query = SeqRecord(Seq(seq_str), id=seq_id, description="QUERY")
    write_panel_fasta([query] + panel, combined)

    if use_linsi:
        run_mafft_linsi(combined, raw_aln, threads, force=True)
    else:
        run_mafft_auto(combined, raw_aln, threads, force=True)
    run_trimal(raw_aln, final_aln, force=True)
    # Clean up intermediate artifacts; keep .aln only.
    combined.unlink(missing_ok=True)
    raw_aln.unlink(missing_ok=True)
    return final_aln


def process_spidroin(
    row: dict,
    pred_seq: str,
    species_msa_dir: Path,
    ref_dir: Path,
    threads: int,
    force: bool,
) -> dict:
    """Generate the five .aln files for one spidroin entry; return a status dict."""
    sid = row["spidroin_id"]
    spt = row.get("spidroin_type") or ""
    ntd_panel_strict = _load_reference_panel(ref_dir, spt, "NTD")
    ctd_panel_strict = _load_reference_panel(ref_dir, spt, "CTD")
    ntd_panel_loose = build_loose_panel(ref_dir, spt, "NTD")
    ctd_panel_loose = build_loose_panel(ref_dir, spt, "CTD")

    workdir = species_msa_dir / sid
    workdir.mkdir(parents=True, exist_ok=True)

    ntd_seq = slice_domain(pred_seq, row.get("hmm_ntd_from"),
                           row.get("hmm_ntd_to"), fallback="ntd")
    ctd_seq = slice_domain(pred_seq, row.get("hmm_ctd_from"),
                           row.get("hmm_ctd_to"), fallback="ctd")

    status = {"spidroin_id": sid}
    status["ntd_strict"] = bool(make_msa(sid, ntd_seq, ntd_panel_strict,
                                          workdir, "ntd_strict",
                                          threads, False, force))
    status["ntd_loose"] = bool(make_msa(sid, ntd_seq, ntd_panel_loose,
                                         workdir, "ntd_loose",
                                         threads, False, force))
    status["ctd_strict"] = bool(make_msa(sid, ctd_seq, ctd_panel_strict,
                                          workdir, "ctd_strict",
                                          threads, False, force))
    status["ctd_loose"] = bool(make_msa(sid, ctd_seq, ctd_panel_loose,
                                         workdir, "ctd_loose",
                                         threads, False, force))
    # Overview MSA: full-length predicted protein vs same-type panel.
    overview_panel = _load_reference_panel(ref_dir, spt, "NTD") + \
                     _load_reference_panel(ref_dir, spt, "CTD")
    if len(overview_panel) >= _MIN_PANEL_SIZE and pred_seq:
        make_msa(sid, pred_seq, overview_panel, workdir,
                 "overview", threads, False, force)
        status["overview"] = (workdir / "overview.aln").exists()
    else:
        status["overview"] = False
    return status


def process_species(
    species_name: str,
    confirmation_species_dir: Path,
    msa_species_dir: Path,
    ref_dir: Path,
    threads: int,
    force: bool,
    only_validated: bool = True,
) -> int:
    """Iterate spidroins of a species; return number of processed entries."""
    confirmation_tsv = confirmation_species_dir / "protein_confirmation.tsv"
    pred_fa = confirmation_species_dir / "predicted_proteins.fa"
    if not confirmation_tsv.exists() or not pred_fa.exists():
        logger.warning(f"[{species_name}] missing inputs, skipping")
        return 0

    df = pl.read_csv(confirmation_tsv, separator="\t", infer_schema_length=1000)
    if only_validated and "validation_status" in df.columns:
        df = df.filter(pl.col("validation_status").is_in(["validated", "partial"]))
    if df.is_empty():
        logger.info(f"[{species_name}] no rows pass filter, skipping")
        return 0

    seqs = _load_predicted(pred_fa)
    msa_species_dir.mkdir(parents=True, exist_ok=True)

    n = 0
    for row in df.to_dicts():
        sid = row["spidroin_id"]
        seq = seqs.get(sid, "")
        if not seq:
            continue
        process_spidroin(row, seq, msa_species_dir, ref_dir, threads, force)
        n += 1
    return n


# ── CLI ──────────────────────────────────────────────────────────────────────


@app.command()
def main(
    confirmation_dir: Path = PROCESSED_DATA_DIR / DEFAULT_TASK_NAME / "confirmation",
    msa_dir: Path = PROCESSED_DATA_DIR / DEFAULT_TASK_NAME / "msa",
    ref_dir: Path = DEFAULT_REF_DIR,
    threads: int = 8,
    species: str | None = None,
    only_validated: bool = True,
    force: bool = False,
):
    """
    Generate per-spidroin domain-level MSAs (NTD/CTD strict + loose, overview).
    """
    species_dirs = sorted(d for d in confirmation_dir.iterdir() if d.is_dir())
    if species is not None:
        species_dirs = [d for d in species_dirs if d.name == species]
    if not species_dirs:
        logger.error(f"No species under {confirmation_dir} matching filter")
        raise typer.Exit(1)

    total = 0
    for sp_dir in tqdm(species_dirs, desc="Generating MSAs"):
        n = process_species(
            sp_dir.name, sp_dir,
            msa_dir / sp_dir.name,
            ref_dir, threads, force, only_validated,
        )
        logger.info(f"[{sp_dir.name}] {n} spidroins MSA-aligned")
        total += n
    logger.success(f"MSA generation complete for {total} spidroins")


if __name__ == "__main__":
    app()
