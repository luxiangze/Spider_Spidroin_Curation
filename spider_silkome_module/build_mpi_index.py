from pathlib import Path

from loguru import logger
from tqdm import tqdm
import typer

from spider_silkome_module.config import RAW_DATA_DIR
from spider_silkome_module.utils.run_cmd import run_cmd

app = typer.Typer()


def extract_species_name(folder_name: str) -> str:
    """
    Extract species name from folder name by removing the numeric prefix.
    e.g. '001.Allagelena_difficilis' -> 'Allagelena_difficilis'
    """
    if "." in folder_name:
        return folder_name.split(".", 1)[1]
    return folder_name


def find_genome_files(input_path: Path) -> list[tuple[str, Path]]:
    """
    Find all genome FASTA files in subdirectories of input_path.
    Returns a list of (species_name, genome_fasta_path) tuples.
    """
    genome_files = []
    for subdir in sorted(input_path.iterdir()):
        if not subdir.is_dir():
            continue

        fasta_list = list(subdir.glob("*.fa")) + list(subdir.glob("*.fna"))
        if fasta_list:
            genome_files.append((subdir.name, fasta_list[0]))
        else:
            logger.warning(f"No genome FASTA found in {subdir}")

    return genome_files


@app.command()
def main(
    input_path: Path = RAW_DATA_DIR / "01.ref_gff",
    output_path: Path = RAW_DATA_DIR / "spider_genome",
    threads: int = 70,
    force: bool = False,
):
    """
    Batch build miniprot .mpi index files from genome FASTA files.

    Each output .mpi file is named by species name (e.g. Allagelena_difficilis.mpi).
    """
    logger.info(f"Input path: {input_path}")
    logger.info(f"Output path: {output_path}")

    genome_files = find_genome_files(input_path)
    if not genome_files:
        logger.error(f"No genome FASTA files found in {input_path}")
        raise typer.Exit(1)

    logger.info(f"Found {len(genome_files)} genome(s) to index")
    output_path.mkdir(parents=True, exist_ok=True)

    for species_name, genome_fasta in tqdm(genome_files, desc="Building MPI indexes"):
        mpi_output = output_path / f"{species_name}.mpi"
        logger.info(f"Indexing {species_name}: {genome_fasta} -> {mpi_output}")

        cmd = f"miniprot -t {threads} -d {mpi_output} {genome_fasta}"
        run_cmd(cmd, [mpi_output], force=force)

    logger.success(f"Completed MPI index building for {len(genome_files)} genome(s)")


if __name__ == "__main__":
    app()
