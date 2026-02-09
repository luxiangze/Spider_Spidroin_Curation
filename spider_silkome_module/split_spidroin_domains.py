"""Split spidroin CTD and NTD sequences from FASTA file into separate files."""

from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from loguru import logger
import typer

from spider_silkome_module.config import EXTERNAL_DATA_DIR, INTERIM_DATA_DIR

app = typer.Typer()


def parse_fasta(fasta_path: Path) -> dict[str, list[SeqRecord]]:
    """Parse FASTA file and group sequences by spidroin type and domain.

    Args:
        fasta_path: Path to the input FASTA file.

    Returns:
        Dictionary mapping (spidroin_type_domain) to list of SeqRecord.
    """
    sequences: dict[str, list[SeqRecord]] = defaultdict(list)

    for record in SeqIO.parse(fasta_path, "fasta"):
        parts = record.id.split("|")
        if len(parts) < 8:
            logger.warning(f"Skipping malformed header: {record.id}")
            continue

        spidroin_type = parts[5]
        domain = parts[7]
        key = f"{spidroin_type}_{domain}"
        sequences[key].append(record)

    return sequences


def write_fasta_files(
    sequences: dict[str, list[SeqRecord]],
    output_dir: Path,
) -> None:
    """Write sequences to separate FASTA files.

    Args:
        sequences: Dictionary mapping key to list of SeqRecord.
        output_dir: Directory to write output files.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    for key, records in sequences.items():
        output_file = output_dir / f"{key}.fa"
        SeqIO.write(records, output_file, "fasta")
        logger.info(f"Wrote {len(records)} sequences to {output_file.name}")


@app.command()
def main(
    input_path: Path = EXTERNAL_DATA_DIR
    / "spider-silkome-database.v1.prot.fixed.fasta",
    output_dir: Path = INTERIM_DATA_DIR / "spidroin_proteins",
) -> None:
    """Split spidroin sequences by type and domain (CTD/NTD).

    Args:
        input_path: Path to input FASTA file.
        output_dir: Directory to write output files.
    """
    logger.info(f"Reading sequences from {input_path}")

    if not input_path.exists():
        logger.error(f"Input file not found: {input_path}")
        raise typer.Exit(code=1)

    sequences = parse_fasta(input_path)
    logger.info(f"Found {len(sequences)} spidroin type-domain combinations")

    write_fasta_files(sequences, output_dir)
    logger.success(f"Split complete. Output files in {output_dir}")


if __name__ == "__main__":
    app()
