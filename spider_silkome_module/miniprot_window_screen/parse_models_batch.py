from __future__ import annotations

from pathlib import Path
from typing import Any

from loguru import logger
from tqdm import tqdm
import typer

from spider_silkome_module.miniprot_window_screen.batch_common import (
    fail_if_any_failed,
    read_tsv,
    summary_path,
    write_tsv,
)
from spider_silkome_module.parse_miniprot_window_gff import (
    build_models,
    write_jsonl,
    write_summary,
)

app = typer.Typer(help="Batch parse miniprot mixed output into mRNA model JSONL/summary files.")


@app.command()
def main(
    species_manifest: Path = typer.Option(..., "--species-manifest"),
    interim_root: Path = typer.Option(..., "--interim-root"),
    processed_root: Path = typer.Option(..., "--processed-root"),
    query_proteins: Path = typer.Option(..., "--query-proteins"),
    force: bool = False,
) -> None:
    del processed_root
    rows = read_tsv(species_manifest)
    summary_rows: list[dict[str, Any]] = []
    logger.info(f"Parsing miniprot models for {len(rows)} species")

    for row in tqdm(rows, desc="Parse models"):
        species = row["species"]
        mixed = Path(row["miniprot_mixed"])
        window_fasta = Path(row["window_fasta"])
        jsonl = Path(row["models_jsonl"])
        summary = Path(row["models_summary"])
        try:
            if not mixed.exists() or not window_fasta.exists():
                status = "skipped"
                reason = "missing miniprot output or window FASTA"
                n_models = 0
            elif not force and jsonl.exists() and summary.exists():
                status = "skipped"
                reason = "outputs exist"
                n_models = ""
            else:
                models = build_models(mixed, window_fasta, query_proteins)
                models.sort(key=lambda model: (model.window_id, model.genomic_start, model.genomic_end, model.mrna_id))
                write_jsonl(models, jsonl)
                write_summary(models, summary)
                status = "ok"
                reason = ""
                n_models = len(models)
            summary_rows.append({
                "species": species,
                "status": status,
                "reason": reason,
                "input": str(mixed),
                "models_jsonl": str(jsonl),
                "models_summary": str(summary),
                "n_models": n_models,
            })
        except Exception as exc:
            logger.exception(f"Model parsing failed for {species}")
            summary_rows.append({
                "species": species,
                "status": "failed",
                "reason": repr(exc),
                "input": str(mixed),
                "models_jsonl": str(jsonl),
                "models_summary": str(summary),
                "n_models": "",
            })

    out = summary_path(interim_root, "parse_models")
    write_tsv(out, summary_rows)
    logger.success(f"Parse models summary: {out}")
    fail_if_any_failed(summary_rows)


if __name__ == "__main__":
    app()
