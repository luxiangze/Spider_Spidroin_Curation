from pathlib import Path

from Bio import SeqIO
from loguru import logger
import typer

from spider_silkome_module.config import EXTERNAL_DATA_DIR

app = typer.Typer()


def parse_header(header: str) -> str:
    """
    Parse the original header and return the new format.

    Original format: >1|7047|Oecobiidae|Oecobius|navus|MiSp|MiSp|CTD
    New format: >1_MiSp_CTD
    """
    parts = header.split("|")
    if len(parts) < 8:
        logger.warning(f"Unexpected header format: {header}")
        return header

    seq_id = parts[0]
    protein_type = parts[6]
    domain_type = parts[7]

    return f"{seq_id}_{protein_type}_{domain_type}"


@app.command()
def main(
    input_path: Path = EXTERNAL_DATA_DIR / "spider-silkome-database.v1.prot.fixed.fasta",
    output_path: Path = EXTERNAL_DATA_DIR / "spider-silkome-database.v1.prot.fixed.renamed.fasta",
):
    """
    Rename FASTA headers from original format to simplified format.

    Original: >1|7047|Oecobiidae|Oecobius|navus|MiSp|MiSp|CTD
    New: >1_MiSp_CTD
    """
    logger.info(f"Reading FASTA from {input_path}")

    records = []
    for record in SeqIO.parse(input_path, "fasta"):
        new_id = parse_header(record.id)
        record.id = new_id
        record.description = ""
        records.append(record)

    logger.info(f"Writing {len(records)} records to {output_path}")
    SeqIO.write(records, output_path, "fasta")

    logger.success(f"Renamed {len(records)} sequences successfully.")


if __name__ == "__main__":
    app()
