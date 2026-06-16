from __future__ import annotations

from pathlib import Path

from loguru import logger
from tqdm import tqdm
import typer

app = typer.Typer(help="Split miniprot --gff mixed output into GFF3, PAF, and unclassified files.")


GFF_FEATURES = {
    "gene",
    "mRNA",
    "CDS",
    "exon",
    "intron",
    "start_codon",
    "stop_codon",
    "five_prime_UTR",
    "three_prime_UTR",
}


def is_int_string(value: str) -> bool:
    try:
        int(value)
        return True
    except ValueError:
        return False


def is_gff_line(fields: list[str]) -> bool:
    return len(fields) >= 9 and fields[2] in GFF_FEATURES


def is_paf_line(fields: list[str]) -> bool:
    if len(fields) < 12 or fields[4] not in {"+", "-"}:
        return False
    int_cols = [1, 2, 3, 6, 7, 8, 9, 10, 11]
    return all(is_int_string(fields[i]) for i in int_cols)


def split_miniprot_mixed_output(
    mixed_path: Path,
    gff_path: Path,
    paf_path: Path,
    unclassified_path: Path,
) -> dict[str, int]:
    counts = {"gff": 0, "paf": 0, "unclassified": 0}
    gff_path.parent.mkdir(parents=True, exist_ok=True)
    paf_path.parent.mkdir(parents=True, exist_ok=True)
    unclassified_path.parent.mkdir(parents=True, exist_ok=True)

    with mixed_path.open() as inp, gff_path.open("w") as gff, paf_path.open("w") as paf, unclassified_path.open("w") as unk:
        for line in tqdm(inp, desc=f"Split {mixed_path.name}"):
            if not line.strip():
                continue
            if line.startswith("#"):
                gff.write(line)
                counts["gff"] += 1
                continue
            fields = line.rstrip("\n").split("\t")
            if is_gff_line(fields):
                gff.write(line)
                counts["gff"] += 1
            elif is_paf_line(fields):
                paf.write(line)
                counts["paf"] += 1
            else:
                unk.write(line)
                counts["unclassified"] += 1
    return counts


@app.command()
def main(
    mixed_path: Path,
    gff_path: Path,
    paf_path: Path,
    unclassified_path: Path,
) -> None:
    logger.info(f"Splitting miniprot output: {mixed_path}")
    counts = split_miniprot_mixed_output(mixed_path, gff_path, paf_path, unclassified_path)
    logger.success(f"Split counts: {counts}")


if __name__ == "__main__":
    app()

