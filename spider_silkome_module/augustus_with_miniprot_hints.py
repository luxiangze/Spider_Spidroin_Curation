"""
augustus_with_miniprot_hints.py — refine spidroin gene structure with Augustus
using miniprot CDS rows as strong P (protein) hints.

Input per species:
  miniprot_raw_dir/<species>.gff           (full-length protein miniprot output)
  typing_dir/<species>/spidroin_full_length.fasta

Output per species (under output_dir/<species>/):
  augustus_output.gff
  augustus_output.aa

Hints are written to:
  hints_dir/<species>/hints_from_miniprot.gff

Hint translation rules:
  - first CDS block start  → 'start' hint
  - last CDS block end     → 'stop' hint
  - each CDS block         → 'CDSpart' hint
  - between consecutive CDS blocks → 'intron' hint (if gap > 0)
"""

from pathlib import Path

from loguru import logger
from tqdm import tqdm
import typer

from spider_silkome_module.config import EXTERNAL_DATA_DIR, INTERIM_DATA_DIR, PROCESSED_DATA_DIR
from spider_silkome_module.protein_confirmation import (
    parse_miniprot_gff,
    select_best_per_spidroin,
)
from spider_silkome_module.utils.run_cmd import run_cmd

app = typer.Typer()

DEFAULT_TASK_NAME = "protein_confirmation"
FASTA_FILE = "spidroin_full_length.fasta"
HINTS_FILE = "hints_from_miniprot.gff"
GFF_OUTPUT = "augustus_output.gff"
AA_OUTPUT = "augustus_output.aa"


def write_hints(species_gff: Path, hints_path: Path) -> int:
    """
    Convert miniprot best-hit CDS rows into Augustus CDS/intron/start/stop
    hints with source 'P' (protein).  Returns number of hint lines written.
    """
    records_per_id = parse_miniprot_gff(species_gff)
    best = select_best_per_spidroin(records_per_id)
    hints_path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(hints_path, "w") as out:
        for sid, r in sorted(best.items()):
            if not r.cds_blocks:
                continue
            first_s, _, _ = r.cds_blocks[0]
            _, last_e, _ = r.cds_blocks[-1]
            # start codon hint
            out.write(
                f"{sid}\tminiprot\tstart\t{first_s}\t{first_s + 2}\t"
                f"{r.score}\t{r.strand}\t.\tsrc=P;grp={sid};pri=4\n"
            )
            n += 1
            # stop codon hint
            out.write(
                f"{sid}\tminiprot\tstop\t{last_e - 2}\t{last_e}\t"
                f"{r.score}\t{r.strand}\t.\tsrc=P;grp={sid};pri=4\n"
            )
            n += 1
            # CDS-part + intron hints
            for s, e, _ in r.cds_blocks:
                out.write(
                    f"{sid}\tminiprot\tCDSpart\t{s}\t{e}\t"
                    f"{r.score}\t{r.strand}\t.\tsrc=P;grp={sid};pri=4\n"
                )
                n += 1
            for i in range(len(r.cds_blocks) - 1):
                istart = r.cds_blocks[i][1] + 1
                iend = r.cds_blocks[i + 1][0] - 1
                if iend < istart:
                    continue
                out.write(
                    f"{sid}\tminiprot\tintron\t{istart}\t{iend}\t"
                    f"{r.score}\t{r.strand}\t.\tsrc=P;grp={sid};pri=4\n"
                )
                n += 1
    return n


def run_augustus(
    fasta: Path,
    hints: Path,
    extrinsic_cfg: Path,
    augustus_species: str,
    gff_out: Path,
    force: bool = False,
) -> None:
    """Invoke augustus with hints; write GFF3."""
    hints_arg = (f"--hintsfile={hints} --extrinsicCfgFile={extrinsic_cfg}"
                 if hints.exists() and hints.stat().st_size > 0 else "")
    cmd = (
        f"augustus --strand=both --singlestrand=true "
        f"--alternatives-from-evidence=true --gff3=on --uniqueGeneId=true "
        f"--genemodel=partial --UTR=off --species={augustus_species} "
        f"{hints_arg} {fasta} > {gff_out}"
    )
    run_cmd(cmd, [gff_out], force=force)


def extract_protein(fasta: Path, gff_out: Path, aa_out: Path,
                    force: bool = False) -> None:
    if not gff_out.exists():
        logger.warning(f"No {gff_out.name} found, skipping protein extraction")
        return
    cmd = f"getAnnoFasta.pl --seqfile={fasta} {gff_out}"
    run_cmd(cmd, [aa_out], force=force)


@app.command()
def main(
    miniprot_raw_dir: Path = INTERIM_DATA_DIR / DEFAULT_TASK_NAME / "miniprot_raw",
    typing_dir: Path = PROCESSED_DATA_DIR / "typing_results",
    hints_dir: Path = INTERIM_DATA_DIR / DEFAULT_TASK_NAME / "augustus_hints",
    output_dir: Path = PROCESSED_DATA_DIR / DEFAULT_TASK_NAME / "augustus",
    extrinsic_cfg: Path = EXTERNAL_DATA_DIR / "extrinsic.cfg",
    augustus_species: str = "parasteatoda",
    species: str | None = None,
    force: bool = False,
):
    """Run Augustus refinement using miniprot CDS as hints, per species."""
    gff_files = sorted(miniprot_raw_dir.glob("*.gff"))
    if species is not None:
        gff_files = [g for g in gff_files if g.stem == species]
    if not gff_files:
        logger.error(f"No miniprot GFF files in {miniprot_raw_dir} matching filter")
        raise typer.Exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)
    hints_dir.mkdir(parents=True, exist_ok=True)

    for gff in tqdm(gff_files, desc="Augustus refinement"):
        species_name = gff.stem
        fasta = typing_dir / species_name / FASTA_FILE
        if not fasta.exists():
            logger.warning(f"[{species_name}] no fasta, skipping")
            continue
        sp_hints = hints_dir / species_name / HINTS_FILE
        n = write_hints(gff, sp_hints)
        logger.info(f"[{species_name}] {n} hint lines written")

        sp_out = output_dir / species_name
        sp_out.mkdir(parents=True, exist_ok=True)
        gff_out = sp_out / GFF_OUTPUT
        aa_out = sp_out / AA_OUTPUT
        run_augustus(fasta, sp_hints, extrinsic_cfg, augustus_species,
                     gff_out, force=force)
        extract_protein(fasta, gff_out, aa_out, force=force)

    logger.success("Augustus refinement complete.")


if __name__ == "__main__":
    app()
