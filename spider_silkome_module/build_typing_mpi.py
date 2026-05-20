"""
build_typing_mpi.py — build miniprot mpi indexes from typing_results.

For each species directory under typing_results/, locates spidroin_full_length.fasta
and builds a miniprot .mpi index named after the species directory.

Output:
    {output_path}/{species_dir_name}.mpi
"""

from pathlib import Path

from loguru import logger
from tqdm import tqdm
import typer

from spider_silkome_module.config import INTERIM_DATA_DIR, PROCESSED_DATA_DIR
from spider_silkome_module.utils.run_cmd import run_cmd

app = typer.Typer()

FASTA_FILE = "spidroin_full_length.fasta"


def find_spidroin_fastas(typing_dir: Path) -> list[tuple[str, Path]]:
    """
    Find spidroin_full_length.fasta in each species subdirectory of typing_dir.
    Returns a list of (species_dir_name, fasta_path) tuples.
    Skips species directories that do not contain a fasta file or whose fasta is empty.
    """
    results: list[tuple[str, Path]] = []
    for sp_dir in sorted(typing_dir.iterdir()):
        if not sp_dir.is_dir():
            continue
        fa = sp_dir / FASTA_FILE
        if not fa.exists() or fa.stat().st_size == 0:
            logger.warning(f"No usable {FASTA_FILE} in {sp_dir.name}, skipping")
            continue
        results.append((sp_dir.name, fa))
    return results


@app.command()
def main(
    typing_dir: Path = PROCESSED_DATA_DIR / "typing_results",
    output_path: Path = INTERIM_DATA_DIR / "protein_confirmation" / "mpi",
    threads: int = 70,
    force: bool = False,
):
    """
    Build miniprot .mpi index for each species' spidroin_full_length.fasta.
    """
    logger.info(f"Typing dir: {typing_dir}")
    logger.info(f"Output path: {output_path}")

    species_fastas = find_spidroin_fastas(typing_dir)
    if not species_fastas:
        logger.error(f"No fasta files found under {typing_dir}")
        raise typer.Exit(1)

    logger.info(f"Found {len(species_fastas)} species to index")
    output_path.mkdir(parents=True, exist_ok=True)

    for species_name, fa in tqdm(species_fastas, desc="Building MPI indexes"):
        mpi_output = output_path / f"{species_name}.mpi"
        cmd = f"miniprot -t {threads} -d {mpi_output} {fa}"
        run_cmd(cmd, [mpi_output], force=force)

    logger.success(f"Completed MPI index building for {len(species_fastas)} species")


if __name__ == "__main__":
    app()
