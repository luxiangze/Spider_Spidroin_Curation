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
from spider_silkome_module.miniprot_window_screen.build_tracks import main as build_tracks_main

app = typer.Typer(help="Batch build pyGenomeTracks GTF/INI/command files.")


@app.command()
def main(
    species_manifest: Path = typer.Option(..., "--species-manifest"),
    interim_root: Path = typer.Option(..., "--interim-root"),
    processed_root: Path = typer.Option(..., "--processed-root"),
    region_flank: int = typer.Option(1000, "--region-flank", min=0),
    clip_region_to_window: bool = typer.Option(True, "--clip-region-to-window/--no-clip-region-to-window"),
    track_label_fraction: float = typer.Option(0.35, "--track-label-fraction", min=0.0, max=1.0),
    per_mrna_track_height: float = typer.Option(0.65, "--per-mrna-track-height", min=0.1),
    plot_font_size: int = typer.Option(7, "--plot-font-size", min=1),
    plot_width: int = typer.Option(42, "--plot-width", min=1),
    force: bool = False,
) -> None:
    del processed_root
    rows = read_tsv(species_manifest)
    summary_rows: list[dict[str, Any]] = []
    logger.info(f"Building tracks for {len(rows)} species")

    for row in tqdm(rows, desc="Build tracks"):
        species = row["species"]
        models_jsonl = Path(row["models_jsonl"])
        miniprot_evidence = Path(row["miniprot_evidence_gff"])
        windows_bed = Path(row["merged_bed"])
        output_dir = Path(row["interim_dir"])
        commands_tsv = Path(row["pygt_commands"])
        try:
            if not models_jsonl.exists() or not windows_bed.exists() or not miniprot_evidence.exists():
                status = "skipped"
                reason = "missing models, windows, or miniprot evidence"
            elif not force and commands_tsv.exists():
                status = "skipped"
                reason = "output exists"
            else:
                build_tracks_main(
                    models_jsonl,
                    miniprot_evidence,
                    windows_bed,
                    output_dir,
                    nhmmer_gff=existing_path(row.get("nhmmer_gff")),
                    min_positive=None,
                    region_flank=region_flank,
                    clip_region_to_window=clip_region_to_window,
                    track_label_fraction=track_label_fraction,
                    per_mrna_track_height=per_mrna_track_height,
                    plot_font_size=plot_font_size,
                    plot_width=plot_width,
                    selected_mpid_tsv=existing_path(row.get("selected_tsv")),
                )
                status = "ok"
                reason = ""
            summary_rows.append({
                "species": species,
                "status": status,
                "reason": reason,
                "commands_tsv": str(commands_tsv),
                "n_plot_commands": count_tsv_records(commands_tsv),
            })
        except Exception as exc:
            logger.exception(f"Build tracks failed for {species}")
            summary_rows.append({
                "species": species,
                "status": "failed",
                "reason": repr(exc),
                "commands_tsv": str(commands_tsv),
                "n_plot_commands": "",
            })

    out = summary_path(interim_root, "build_tracks")
    write_tsv(out, summary_rows)
    logger.success(f"Build tracks summary: {out}")
    fail_if_any_failed(summary_rows)


if __name__ == "__main__":
    app()
