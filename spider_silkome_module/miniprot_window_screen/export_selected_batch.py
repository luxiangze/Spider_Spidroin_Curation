from __future__ import annotations

from pathlib import Path
from typing import Any

from loguru import logger
from tqdm import tqdm
import typer

from spider_silkome_module.miniprot_window_screen.batch_common import (
    count_fasta_records,
    count_tsv_records,
    fail_if_any_failed,
    read_tsv,
    summary_path,
    write_tsv,
)
from spider_silkome_module.miniprot_window_screen.export_selected import (
    annotate_selected_stop_fields,
    load_models,
    load_selected_rows,
    write_selected_fasta,
    write_selected_gff,
    write_selected_rows,
)

app = typer.Typer(help="Batch export selected miniprot models to GFF3 and protein FASTA.")


@app.command()
def main(
    species_manifest: Path = typer.Option(..., "--species-manifest"),
    interim_root: Path = typer.Option(..., "--interim-root"),
    processed_root: Path = typer.Option(..., "--processed-root"),
    force: bool = False,
) -> None:
    del processed_root
    rows = read_tsv(species_manifest)
    summary_rows: list[dict[str, Any]] = []
    logger.info(f"Exporting selected models for {len(rows)} species")

    for row in tqdm(rows, desc="Export selected"):
        species = row["species"]
        selected_tsv = Path(row["selected_tsv"])
        models_jsonl = Path(row["models_jsonl"])
        miniprot_mixed = Path(row["miniprot_mixed"])
        output_gff = Path(row["selected_gff"])
        output_faa = Path(row["selected_faa"])
        try:
            if not selected_tsv.exists() or not models_jsonl.exists() or not miniprot_mixed.exists():
                status = "skipped"
                reason = "missing selected TSV, models, or miniprot output"
                n_gff = 0
                n_fasta = 0
            elif not force and output_gff.exists() and output_faa.exists():
                status = "skipped"
                reason = "outputs exist"
                n_gff = ""
                n_fasta = count_fasta_records(output_faa)
            else:
                selected_rows = load_selected_rows(selected_tsv)
                models_by_id = load_models(models_jsonl)
                selected_rows = annotate_selected_stop_fields(selected_rows, models_by_id)
                write_selected_rows(selected_tsv, selected_rows)
                n_gff = write_selected_gff(selected_rows, miniprot_mixed, output_gff)
                n_fasta = write_selected_fasta(selected_rows, models_by_id, output_faa)
                status = "ok"
                reason = ""
            summary_rows.append({
                "species": species,
                "status": status,
                "reason": reason,
                "selected_tsv": str(selected_tsv),
                "selected_gff": str(output_gff),
                "selected_faa": str(output_faa),
                "n_selected": count_tsv_records(selected_tsv),
                "n_gff_records": n_gff,
                "n_fasta_records": n_fasta,
            })
        except Exception as exc:
            logger.exception(f"Export failed for {species}")
            summary_rows.append({
                "species": species,
                "status": "failed",
                "reason": repr(exc),
                "selected_tsv": str(selected_tsv),
                "selected_gff": str(output_gff),
                "selected_faa": str(output_faa),
                "n_selected": "",
                "n_gff_records": "",
                "n_fasta_records": "",
            })

    out = summary_path(interim_root, "export_selected")
    write_tsv(out, summary_rows)
    logger.success(f"Export summary: {out}")
    fail_if_any_failed(summary_rows)


if __name__ == "__main__":
    app()
