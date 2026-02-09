from pathlib import Path

from loguru import logger
from tqdm import tqdm
import typer

from spider_silkome_module.config import INTERIM_DATA_DIR, RAW_DATA_DIR, REFERENCES_DIR
from spider_silkome_module.utils.run_cmd import run_cmd

app = typer.Typer()


def find_genome_files(input_path: Path) -> list[tuple[str, Path]]:
    """
    Find all genome fasta files in the input path.
    Returns a list of (species_name, genome_path) tuples.

    Supports both:
    - Single file input
    - Directory with subdirectories (one species per subdirectory)
    """
    genome_files = []

    if input_path.is_file() and input_path.suffix in [".fa", ".fasta", ".fna"]:
        species_name = input_path.stem
        genome_files.append((species_name, input_path))
    elif input_path.is_dir():
        for subdir in input_path.iterdir():
            if subdir.is_dir():
                folder_name = subdir.name
                if "." in folder_name:
                    species_name = folder_name.split(".", 1)[1]
                else:
                    species_name = folder_name

                fa_list = list(subdir.glob("*.fa")) + list(subdir.glob("*.fasta")) + list(subdir.glob("*.fna"))
                if fa_list:
                    genome_files.append((species_name, fa_list[0]))
                else:
                    logger.warning(f"No genome file found in {subdir}")

    return genome_files


def press_hmm_models(hmm_dir: Path, force: bool = False) -> None:
    """
    Press HMM models using hmmpress.
    """
    hmm_files = list(hmm_dir.glob("*TD.hmm"))
    logger.info(f"Found {len(hmm_files)} HMM models to press")

    for hmm_file in tqdm(hmm_files, desc="Pressing HMM models"):
        h3m_file = hmm_file.with_suffix(".h3m")

        cmd = f"hmmpress {hmm_file}"
        run_cmd(cmd, [h3m_file], force=force)


def run_nhmmer(
    genome_path: Path,
    hmm_file: Path,
    output_dir: Path,
    threads: int = 70,
    force: bool = False,
) -> None:
    """
    Run nhmmer search for a single genome and HMM model.
    """
    model_name = hmm_file.stem.split(".")[0]
    tbl_file = output_dir / f"{model_name}.tbl"
    out_file = output_dir / f"{model_name}.out"

    cmd = f"nhmmer --cpu {threads} --tblout {tbl_file} {hmm_file} {genome_path} > {out_file}"
    run_cmd(cmd, [out_file], force=force)


@app.command()
def main(
    genome_path: Path = RAW_DATA_DIR / "spider_genome",
    hmm_dir: Path = REFERENCES_DIR / "2025_Schoneberg_data" / "hmmer_nucl_profile_trimmed",
    output_path: Path = INTERIM_DATA_DIR / "automated_spidroin_annotation" / "nhmmer_search",
    threads: int = 70,
    force: bool = False,
    skip_press: bool = False,
):
    """
    Run nhmmer search for one or multiple genomes.

    If genome_path is a directory, it will process all subdirectories containing genome files.
    If genome_path is a single fasta file, it will process only that file.
    """
    logger.info(f"Genome path: {genome_path}")
    logger.info(f"HMM directory: {hmm_dir}")
    logger.info(f"Output path: {output_path}")

    if not skip_press:
        logger.info("Pressing HMM models...")
        press_hmm_models(hmm_dir, force=force)

    genome_files = find_genome_files(genome_path)
    if not genome_files:
        logger.error(f"No genome files found in {genome_path}")
        raise typer.Exit(1)

    logger.info(f"Found {len(genome_files)} genome(s) to process")

    hmm_files = list(hmm_dir.glob("*TD.hmm"))
    logger.info(f"Found {len(hmm_files)} HMM models")

    output_path.mkdir(parents=True, exist_ok=True)

    for species_name, genome_file in tqdm(genome_files, desc="Processing genomes"):
        logger.info(f"Processing {species_name}: {genome_file}")

        species_output_dir = output_path / species_name
        species_output_dir.mkdir(parents=True, exist_ok=True)

        for hmm_file in hmm_files:
            run_nhmmer(genome_file, hmm_file, species_output_dir, threads, force)

        logger.info(f"Completed: {species_name}")

    logger.success(f"Completed nhmmer search for {len(genome_files)} genome(s)")


if __name__ == "__main__":
    app()
