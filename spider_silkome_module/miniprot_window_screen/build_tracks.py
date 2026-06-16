from __future__ import annotations

from pathlib import Path

from loguru import logger
from tqdm import tqdm
import typer

from spider_silkome_module.build_pygenometracks_window_tracks import (
    collect_plot_region,
    group_models_by_window,
    load_models,
    load_selected_mpids,
    load_windows,
    parse_gff_rows,
    write_command_table,
    write_miniprot_evidence_gtf,
    write_nhmmer_gtf,
    write_ordered_window_miniprot_gtfs,
    write_selected_miniprot_gtf,
    write_track_ini,
    write_window_miniprot_gtf,
)

app = typer.Typer(help="Build pyGenomeTracks inputs for miniprot window screening.")


@app.command()
def main(
    models_jsonl: Path,
    miniprot_evidence_gff: Path,
    windows_bed: Path,
    output_dir: Path,
    nhmmer_gff: Path | None = None,
    min_positive: float | None = None,
    region_flank: int = typer.Option(1000, "--region-flank", min=0),
    clip_region_to_window: bool = typer.Option(
        True,
        "--clip-region-to-window/--no-clip-region-to-window",
    ),
    track_label_fraction: float = typer.Option(0.35, "--track-label-fraction", min=0.0, max=1.0),
    per_mrna_track_height: float = typer.Option(0.65, "--per-mrna-track-height", min=0.1),
    plot_font_size: int = typer.Option(7, "--plot-font-size", min=1),
    plot_width: int = typer.Option(42, "--plot-width", min=1),
    selected_mpid_tsv: Path | None = typer.Option(None, "--selected-mpid-tsv"),
) -> None:
    logger.info(f"Building pyGenomeTracks configuration in {output_dir}")
    pygt_dir = output_dir / "pygenometracks"
    gtf_dir = pygt_dir / "gtf"
    tracks_dir = pygt_dir / "tracks"
    plots_dir = pygt_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    models = load_models(models_jsonl, min_positive=min_positive)
    selected_lookup = load_selected_mpids(selected_mpid_tsv)
    models_by_window = group_models_by_window(models)
    windows = load_windows(windows_bed)
    nhmmer_rows = parse_gff_rows(nhmmer_gff) if nhmmer_gff and nhmmer_gff.exists() else []
    miniprot_evidence_rows = parse_gff_rows(miniprot_evidence_gff)

    window_miniprot_gtf = gtf_dir / "window_miniprot_models.gtf"
    ordered_window_miniprot_dir = gtf_dir / "window_miniprot_ordered"
    nhmmer_gtf = gtf_dir / "nhmmer_evidence.gtf"
    miniprot_evidence_gtf = gtf_dir / "miniprot_evidence_models.gtf"
    selected_miniprot_gtf = gtf_dir / "selected_miniprot_models.gtf"
    commands_tsv = pygt_dir / "pygenometracks_commands.tsv"

    n_window = write_window_miniprot_gtf(models, window_miniprot_gtf)
    n_selected = write_selected_miniprot_gtf(models, selected_lookup, selected_miniprot_gtf)
    ordered_tracks_by_window = write_ordered_window_miniprot_gtfs(
        models_by_window,
        windows,
        ordered_window_miniprot_dir,
        selected_lookup,
    )
    n_nhmmer = write_nhmmer_gtf(nhmmer_rows, windows, nhmmer_gtf)
    n_evidence = write_miniprot_evidence_gtf(miniprot_evidence_rows, windows, miniprot_evidence_gtf)

    plot_regions = {
        window["window_id"]: collect_plot_region(
            window,
            models_by_window,
            nhmmer_rows,
            miniprot_evidence_rows,
            region_flank,
            clip_region_to_window,
        )
        for window in windows
    }

    for window in tqdm(windows, desc="Write pyGenomeTracks INI"):
        write_track_ini(
            tracks_dir / f"{window['safe_name']}.ini",
            ordered_tracks_by_window.get(window["window_id"], []),
            selected_miniprot_gtf,
            nhmmer_gtf,
            miniprot_evidence_gtf,
            per_mrna_track_height,
        )
    write_command_table(
        windows,
        plot_regions,
        tracks_dir,
        plots_dir,
        commands_tsv,
        track_label_fraction,
        plot_font_size,
        plot_width,
        per_mrna_track_height,
    )

    logger.success(
        f"tracks={commands_tsv}; window_models={n_window}; selected_models={n_selected}; "
        f"nhmmer_hits={n_nhmmer}; evidence_models={n_evidence}"
    )


if __name__ == "__main__":
    app()

