from pathlib import Path

from loguru import logger
from tqdm import tqdm
import typer

from spider_silkome_module.config import EXTERNAL_DATA_DIR, PROCESSED_DATA_DIR
from spider_silkome_module.utils.run_cmd import run_cmd

app = typer.Typer()

FASTA_FILE = "spidroin_full_length.fasta"
HINTS_FILE = "hints.gff"
GFF_OUTPUT = "augustus_output.gff"
AA_OUTPUT = "augustus_output.aa"


def run_augustus_for_species(
    species_dir: Path,
    extrinsic_cfg: Path,
    augustus_species: str,
    force: bool = False,
) -> bool:
    """Run Augustus gene prediction for a single species directory."""
    fasta = species_dir / FASTA_FILE
    hints = species_dir / HINTS_FILE
    gff_out = species_dir / GFF_OUTPUT

    if not fasta.exists():
        logger.warning(f"No {FASTA_FILE} found in {species_dir}, skipping")
        return False

    hints_arg = f"--hintsfile={hints} --extrinsicCfgFile={extrinsic_cfg}" if hints.exists() else ""

    cmd = (
        f"augustus --strand=forward --singlestrand=true "
        f"--alternatives-from-evidence=true --gff3=on --uniqueGeneId=true "
        f"--genemodel=exactlyone --UTR=off --species={augustus_species} "
        f"{hints_arg} {fasta} > {gff_out}"
    )
    run_cmd(cmd, [gff_out], force=force)
    return True


def extract_protein(species_dir: Path, force: bool = False) -> None:
    """Extract protein sequences from Augustus GFF output."""
    fasta = species_dir / FASTA_FILE
    gff_out = species_dir / GFF_OUTPUT
    aa_out = species_dir / AA_OUTPUT

    if not gff_out.exists():
        logger.warning(f"No {GFF_OUTPUT} found in {species_dir}, skipping protein extraction")
        return

    cmd = f"getAnnoFasta.pl --seqfile={fasta} {gff_out}"
    run_cmd(cmd, [aa_out], force=force)


def collect_species_dirs(input_path: Path) -> list[Path]:
    """
    Collect species directories to process.
    Supports both a single species directory and a parent directory
    containing multiple species subdirectories.
    """
    if (input_path / FASTA_FILE).exists() or (input_path / HINTS_FILE).exists():
        return [input_path]
    return sorted(d for d in input_path.iterdir() if d.is_dir())


@app.command()
def main(
    input_path: Path = PROCESSED_DATA_DIR / "typing_results",
    extrinsic_cfg: Path = EXTERNAL_DATA_DIR / "extrinsic.cfg",
    augustus_species: str = "parasteatoda",
    force: bool = False,
):
    """
    Run Augustus gene prediction on spidroin sequences for each species.

    input_path can be the typing_results root directory (all species)
    or a single species directory (for testing).
    """
    species_dirs = collect_species_dirs(input_path)
    logger.info(f"Found {len(species_dirs)} species to process")

    for species_dir in tqdm(species_dirs, desc="Running Augustus"):
        logger.info(f"Processing: {species_dir.name}")
        ok = run_augustus_for_species(species_dir, extrinsic_cfg, augustus_species, force)
        if ok:
            extract_protein(species_dir, force)

    logger.success("Augustus gene prediction complete.")


if __name__ == "__main__":
    app()
