from __future__ import annotations

from pathlib import Path

from loguru import logger
from tqdm import tqdm
import typer

from spider_silkome_module.parse_miniprot_window_gff import (
    MiniprotCdsBlock,
    MiniprotWindowGene,
    build_models,
    write_jsonl,
    write_summary,
)

app = typer.Typer(help="Parse miniprot window GFF into mRNA dataclass models and summaries.")


@app.command()
def main(
    miniprot_gff: Path,
    window_fasta: Path,
    query_fasta: Path,
    output_jsonl: Path,
    output_summary: Path,
) -> None:
    logger.info(f"Parsing miniprot models: {miniprot_gff}")
    models = build_models(miniprot_gff, window_fasta, query_fasta)
    for _ in tqdm(range(1), total=1, desc="Write model outputs"):
        models.sort(key=lambda model: (model.window_id, model.genomic_start, model.genomic_end, model.mrna_id))
        write_jsonl(models, output_jsonl)
        write_summary(models, output_summary)
    logger.success(f"Parsed {len(models)} miniprot mRNA models")
    logger.info(f"JSONL: {output_jsonl}")
    logger.info(f"Summary: {output_summary}")


__all__ = [
    "MiniprotCdsBlock",
    "MiniprotWindowGene",
    "build_models",
    "write_jsonl",
    "write_summary",
]


if __name__ == "__main__":
    app()

