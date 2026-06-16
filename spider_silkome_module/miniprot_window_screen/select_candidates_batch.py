from __future__ import annotations

from pathlib import Path
from typing import Any

from loguru import logger
from tqdm import tqdm
import typer

from spider_silkome_module.miniprot_window_screen.batch_common import (
    count_tsv_records,
    existing_path,
    fail_if_any_failed,
    read_tsv,
    summary_path,
    write_tsv,
)
from spider_silkome_module.miniprot_window_screen.common import load_windows_bed
from spider_silkome_module.miniprot_window_screen.select_candidates import (
    gff_rows,
    load_domain_hits,
    load_models,
    load_typing_loci,
    make_domain_pairs,
    select_candidates,
    write_rows,
)

app = typer.Typer(help="Batch select likely spidroin MPIDs from parsed miniprot models.")


@app.command()
def main(
    species_manifest: Path = typer.Option(..., "--species-manifest"),
    interim_root: Path = typer.Option(..., "--interim-root"),
    processed_root: Path = typer.Option(..., "--processed-root"),
    strong_positive: float = typer.Option(0.60, "--strong-positive"),
    strong_query_coverage: float = typer.Option(0.80, "--strong-query-coverage"),
    rescue_positive: float = typer.Option(0.55, "--rescue-positive"),
    rescue_query_coverage: float = typer.Option(0.70, "--rescue-query-coverage"),
    min_aa: int = typer.Option(1000, "--min-aa"),
    force: bool = False,
) -> None:
    del processed_root
    rows = read_tsv(species_manifest)
    summary_rows: list[dict[str, Any]] = []
    logger.info(f"Selecting candidates for {len(rows)} species")

    for row in tqdm(rows, desc="Select candidates"):
        species = row["species"]
        models_jsonl = Path(row["models_jsonl"])
        windows_bed = Path(row["merged_bed"])
        output_tsv = Path(row["selected_tsv"])
        try:
            if not models_jsonl.exists() or not windows_bed.exists():
                status = "skipped"
                reason = "missing models or windows"
                selection_mode = row.get("selection_mode", "fallback")
                typing_loci_count = 0
                domain_pairs_count = 0
                n_selected = 0
            elif not force and output_tsv.exists():
                status = "skipped"
                reason = "output exists"
                selection_mode = row.get("selection_mode", "fallback")
                typing_loci_count = row.get("typing_locus_count", "")
                domain_pairs_count = ""
                n_selected = count_tsv_records(output_tsv)
            else:
                models = load_models(models_jsonl)
                windows = load_windows_bed(windows_bed)
                nhmmer_gff = existing_path(row.get("nhmmer_gff"))
                domain_hits = load_domain_hits(nhmmer_gff, windows) if nhmmer_gff else []
                typing_loci = load_typing_loci(
                    existing_path(row.get("typing_tsv")),
                    existing_path(row.get("typing_gff")),
                    windows,
                )
                domain_pairs = [] if typing_loci else make_domain_pairs(domain_hits)
                miniprot_evidence = existing_path(row.get("miniprot_evidence_gff"))
                evidence_rows = gff_rows(miniprot_evidence) if miniprot_evidence else []
                selected = select_candidates(
                    species,
                    models,
                    windows,
                    domain_hits,
                    domain_pairs,
                    evidence_rows,
                    typing_loci=typing_loci,
                    strong_positive=strong_positive,
                    strong_query_coverage=strong_query_coverage,
                    rescue_positive=rescue_positive,
                    rescue_query_coverage=rescue_query_coverage,
                    min_aa=min_aa,
                )
                write_rows(selected, output_tsv)
                status = "ok"
                reason = ""
                selection_mode = "typing_guided" if typing_loci else "fallback"
                typing_loci_count = len(typing_loci)
                domain_pairs_count = len(domain_pairs)
                n_selected = len(selected)
            summary_rows.append({
                "species": species,
                "status": status,
                "reason": reason,
                "selection_mode": selection_mode,
                "typing_loci": typing_loci_count,
                "domain_pairs": domain_pairs_count,
                "output": str(output_tsv),
                "n_selected": n_selected,
            })
        except Exception as exc:
            logger.exception(f"Candidate selection failed for {species}")
            summary_rows.append({
                "species": species,
                "status": "failed",
                "reason": repr(exc),
                "selection_mode": row.get("selection_mode", ""),
                "typing_loci": "",
                "domain_pairs": "",
                "output": str(output_tsv),
                "n_selected": "",
            })

    out = summary_path(interim_root, "select_candidates")
    write_tsv(out, summary_rows)
    logger.success(f"Selection summary: {out}")
    fail_if_any_failed(summary_rows)


if __name__ == "__main__":
    app()
