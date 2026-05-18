from pathlib import Path

from loguru import logger
from tqdm import tqdm
import typer

from spider_silkome_module.config import INTERIM_DATA_DIR, RAW_DATA_DIR, REFERENCES_DIR
from spider_silkome_module.utils.run_cmd import run_cmd

app = typer.Typer()


FASTA_SUFFIXES = (".fa", ".fasta", ".fna")
FASTA_GZ_SUFFIXES = tuple(f"{s}.gz" for s in FASTA_SUFFIXES)
# Files in raw genome directories that are NOT the whole-genome assembly.
NON_GENOME_NAMES = {"pep.fa", "mito.fa", "cds.fa", "transcripts.fa"}


def _is_fasta(path: Path) -> bool:
    name = path.name.lower()
    return name.endswith(FASTA_SUFFIXES) or name.endswith(FASTA_GZ_SUFFIXES)


def _pick_genome_in_dir(subdir: Path) -> Path | None:
    """
    Pick the most likely whole-genome FASTA from a species subdirectory.

    Priority:
    1. Files whose name starts with "genome" (handles `genome.softmasked.fa.gz`, etc.)
    2. Any FASTA file not in NON_GENOME_NAMES (excludes pep.fa, mito.fa, ...)
    Compressed (.gz) files are accepted.
    """
    candidates = [p for p in subdir.iterdir() if p.is_file() and _is_fasta(p)]
    if not candidates:
        return None

    genome_named = [p for p in candidates if p.name.lower().startswith("genome")]
    if genome_named:
        # Prefer the softmasked one if present
        softmasked = [p for p in genome_named if "softmasked" in p.name.lower()]
        return softmasked[0] if softmasked else genome_named[0]

    filtered = [p for p in candidates if p.name not in NON_GENOME_NAMES]
    return filtered[0] if filtered else None


def find_genome_files(input_path: Path) -> list[tuple[str, Path]]:
    """
    Find all genome fasta files in the input path.
    Returns a list of (species_name, genome_path) tuples.

    Supports:
    - Single FASTA file input (plain or gzipped)
    - Directory with subdirectories (one species per subdirectory). Subdirectory
      name may be prefixed with an index like `001.Species_name`; the index is stripped.
    """
    genome_files: list[tuple[str, Path]] = []

    if input_path.is_file() and _is_fasta(input_path):
        species_name = input_path.name
        for suffix in FASTA_GZ_SUFFIXES + FASTA_SUFFIXES:
            if species_name.lower().endswith(suffix):
                species_name = species_name[: -len(suffix)]
                break
        genome_files.append((species_name, input_path))
        return genome_files

    if input_path.is_dir():
        for subdir in sorted(input_path.iterdir()):
            if not subdir.is_dir():
                continue
            folder_name = subdir.name
            # Strip numeric prefix like "001." -> "Species_name"
            if "." in folder_name and folder_name.split(".", 1)[0].isdigit():
                species_name = folder_name.split(".", 1)[1]
            else:
                species_name = folder_name

            genome = _pick_genome_in_dir(subdir)
            if genome is not None:
                genome_files.append((species_name, genome))
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
        # hmmpress appends .h3{m,i,f,p} to the original filename (e.g. X.hmm.h3m)
        h3m_file = Path(f"{hmm_file}.h3m")

        cmd = f"hmmpress -f {hmm_file}"
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

    # nhmmer cannot read gzipped FASTA directly; stream via pigz (parallel) when needed.
    if genome_path.suffix == ".gz":
        cmd = (
            f"pigz -dc -p {threads} {genome_path} | "
            f"nhmmer --cpu {threads} --tblout {tbl_file} {hmm_file} - > {out_file}"
        )
    else:
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
