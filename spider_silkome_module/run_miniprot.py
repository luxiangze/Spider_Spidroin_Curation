from pathlib import Path

from loguru import logger
from tqdm import tqdm
import typer

from spider_silkome_module.config import INTERIM_DATA_DIR, RAW_DATA_DIR
from spider_silkome_module.utils.run_cmd import run_cmd

app = typer.Typer()


def find_mpi_files(input_path: Path) -> list[tuple[str, Path]]:
    """
    Find all .mpi files in the input path.
    Returns a list of (species_name, mpi_path) tuples.
    """
    mpi_files = []

    if input_path.is_file() and input_path.suffix == ".mpi":
        species_name = input_path.stem
        mpi_files.append((species_name, input_path))
    elif input_path.is_dir():
        for mpi in sorted(input_path.glob("*.mpi")):
            species_name = mpi.stem
            mpi_files.append((species_name, mpi))

    return mpi_files


def run_miniprot(
    mpi_path: Path,
    protein_fasta: Path,
    output_gff: Path,
    threads: int = 70,
    outc: float = 0.8,
    force: bool = False,
) -> None:
    """
    Run miniprot alignment for a single genome.
    """
    cmd = (
        f"miniprot -S --gff-delim . --outc {outc} -t {threads} -I --gff-only "
        f"{mpi_path} {protein_fasta} > {output_gff}"
    )
    run_cmd(cmd, [output_gff], force=force)


@app.command()
def main(
    mpi_input_path: Path = RAW_DATA_DIR / "spider_genome",
    protein_fasta: Path = INTERIM_DATA_DIR / "miniprot_mapping_20260127" / "cdhit" / "cdhit_shortest_seq.fa",
    output_path: Path = INTERIM_DATA_DIR / "miniprot_mapping_20260127" / "miniprot",
    threads: int = 70,
    outc: float = 0.8,
    force: bool = False,
):
    """
    Run miniprot alignment for one or multiple genomes.

    If mpi_input_path is a directory, it will process all subdirectories containing .mpi files.
    If mpi_input_path is a single .mpi file, it will process only that file.
    """
    logger.info(f"MPI input path: {mpi_input_path}")
    logger.info(f"Protein FASTA: {protein_fasta}")

    mpi_files = find_mpi_files(mpi_input_path)
    if not mpi_files:
        logger.error(f"No .mpi files found in {mpi_input_path}")
        raise typer.Exit(1)

    logger.info(f"Found {len(mpi_files)} genome(s) to process")

    output_path.mkdir(parents=True, exist_ok=True)

    for species_name, mpi_path in tqdm(mpi_files, desc="Processing genomes"):
        logger.info(f"Processing {species_name}: {mpi_path}")
        output_gff = output_path / f"{species_name}.gff"

        run_miniprot(mpi_path, protein_fasta, output_gff, threads, outc, force)
        logger.info(f"Output: {output_gff}")

    logger.success(f"Completed miniprot alignment for {len(mpi_files)} genome(s)")


if __name__ == "__main__":
    app()
