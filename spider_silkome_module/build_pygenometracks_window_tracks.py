from __future__ import annotations

import csv
import json
from pathlib import Path
import re
import shlex
import shutil
from typing import Any

import typer

app = typer.Typer()


Feature = dict[str, Any]
Window = dict[str, Any]


PLOT_FEATURES = {"mRNA", "CDS"}


def parse_attrs(attr_str: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in attr_str.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value
        elif " " in item:
            key, value = item.split(" ", 1)
            attrs[key] = value.strip().strip('"')
    return attrs


def parse_float(value: str | None) -> float | None:
    if value in (None, "", "."):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def gtf_escape(value: Any) -> str:
    text = "" if value is None else str(value)
    return text.replace("\\", "\\\\").replace('"', '\\"')


def format_gtf_attrs(attrs: dict[str, Any]) -> str:
    return " ".join(f'{key} "{gtf_escape(value)}";' for key, value in attrs.items())


def format_metric(value: Any, digits: int = 3) -> str:
    if value is None:
        return "NA"
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return str(value)


def get_protein_sequence_length(model: dict[str, Any]) -> int | None:
    value = model.get("protein_sequence_length")
    if value is not None:
        try:
            return int(value)
        except (TypeError, ValueError):
            pass
    protein_sequence = model.get("protein_sequence")
    if isinstance(protein_sequence, str):
        return len(protein_sequence)
    return None


def write_gtf_line(
    handle,
    seqid: str,
    source: str,
    feature: str,
    start: int,
    end: int,
    score: float | int | None,
    strand: str,
    frame: str,
    attrs: dict[str, Any],
) -> None:
    score_text = "." if score is None else str(score)
    handle.write(
        "\t".join([
            seqid,
            source,
            feature,
            str(start),
            str(end),
            score_text,
            strand,
            frame,
            format_gtf_attrs(attrs),
        ])
        + "\n"
    )


def load_models(models_jsonl: Path, min_positive: float | None = None) -> list[dict[str, Any]]:
    models: list[dict[str, Any]] = []
    with models_jsonl.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            model = json.loads(line)
            if min_positive is not None:
                positive = model.get("positive")
                if positive is None or positive < min_positive:
                    continue
            models.append(model)
    return models


def load_windows(windows_bed: Path) -> list[Window]:
    windows: list[Window] = []
    with windows_bed.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            start0 = int(fields[1])
            end = int(fields[2])
            start1 = start0 + 1
            window_id = f"{chrom}_{start1}-{end}:."
            windows.append({
                "chrom": chrom,
                "start0": start0,
                "start": start1,
                "end": end,
                "window_id": window_id,
                "safe_name": safe_name(f"{chrom}_{start1}-{end}"),
                "region": f"{chrom}:{start1}-{end}",
            })
    return windows


def row_overlaps_window(row: Feature, window: Window) -> bool:
    return (
        row["seqid"] == window["chrom"]
        and row["end"] >= window["start"]
        and row["start"] <= window["end"]
    )


def parse_gff_rows(gff: Path) -> list[Feature]:
    rows: list[Feature] = []
    with gff.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attrs = parse_attrs(fields[8])
            rows.append({
                "seqid": fields[0],
                "source": fields[1],
                "feature": fields[2],
                "start": int(fields[3]),
                "end": int(fields[4]),
                "score": parse_float(fields[5]),
                "strand": fields[6],
                "phase": fields[7],
                "attrs": attrs,
            })
    return rows


def model_sort_key(model: dict[str, Any]) -> tuple[Any, ...]:
    positive = model.get("positive")
    protein_sequence_length = get_protein_sequence_length(model)
    score = model.get("score")
    query_coverage = model.get("query_coverage")
    return (
        -(positive if positive is not None else float("-inf")),
        -(protein_sequence_length if protein_sequence_length is not None else float("-inf")),
        model["genomic_start"],
        -(score if score is not None else float("-inf")),
        -(query_coverage if query_coverage is not None else float("-inf")),
        model["mrna_id"],
    )


def group_models_by_window(models: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    grouped: dict[str, list[dict[str, Any]]] = {}
    for model in models:
        grouped.setdefault(model["window_id"], []).append(model)
    for window_id, window_models in grouped.items():
        grouped[window_id] = sorted(window_models, key=model_sort_key)
    return grouped


def model_common_attrs(model: dict[str, Any]) -> dict[str, Any]:
    gene_id = f"window_gene_{model['mrna_id']}"
    transcript_id = f"window_tx_{model['mrna_id']}"
    return {
        "gene_id": gene_id,
        "transcript_id": transcript_id,
        "gene_name": f"{model.get('target_protein') or model['mrna_id']}|{model['mrna_id']}",
        "transcript_name": f"{model.get('target_protein') or model['mrna_id']}|{model['mrna_id']}",
        "source_window": model["window_id"],
        "target_protein": model.get("target_protein"),
        "identity": model.get("identity"),
        "positive": model.get("positive"),
        "query_coverage": model.get("query_coverage"),
        "exon_count": model.get("exon_count"),
    }


def write_model_gtf(model: dict[str, Any], handle) -> None:
    common_attrs = model_common_attrs(model)
    write_gtf_line(
        handle,
        model["genomic_seqid"],
        "window_miniprot",
        "gene",
        model["genomic_start"],
        model["genomic_end"],
        model.get("score"),
        model["strand"],
        ".",
        common_attrs,
    )
    write_gtf_line(
        handle,
        model["genomic_seqid"],
        "window_miniprot",
        "transcript",
        model["genomic_start"],
        model["genomic_end"],
        model.get("score"),
        model["strand"],
        ".",
        common_attrs,
    )
    for i, block in enumerate(model.get("cds_blocks", []), 1):
        block_attrs = dict(common_attrs)
        block_attrs["exon_number"] = i
        write_gtf_line(
            handle,
            block["genomic_seqid"],
            "window_miniprot",
            "exon",
            block["genomic_start"],
            block["genomic_end"],
            block.get("score"),
            block["strand"],
            ".",
            block_attrs,
        )
        write_gtf_line(
            handle,
            block["genomic_seqid"],
            "window_miniprot",
            "CDS",
            block["genomic_start"],
            block["genomic_end"],
            block.get("score"),
            block["strand"],
            block.get("phase") or ".",
            block_attrs,
        )


def write_window_miniprot_gtf(models: list[dict[str, Any]], output_gtf: Path) -> int:
    output_gtf.parent.mkdir(parents=True, exist_ok=True)
    sorted_models = sorted(
        models,
        key=lambda m: (m["genomic_seqid"], m["genomic_start"], m["genomic_end"], m["mrna_id"]),
    )
    with output_gtf.open("w") as out:
        for model in sorted_models:
            write_model_gtf(model, out)
    return len(sorted_models)


def track_title(model: dict[str, Any]) -> str:
    return (
        f"P={format_metric(model.get('positive'))} "
        f"aa={get_protein_sequence_length(model) or 'NA'} "
        f"{model['mrna_id']}"
    )


def write_ordered_window_miniprot_gtfs(
    models_by_window: dict[str, list[dict[str, Any]]],
    windows: list[Window],
    output_root: Path,
) -> dict[str, list[dict[str, Any]]]:
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    tracks_by_window: dict[str, list[dict[str, Any]]] = {}
    for window in windows:
        window_models = models_by_window.get(window["window_id"], [])
        window_dir = output_root / window["safe_name"]
        window_dir.mkdir(parents=True, exist_ok=True)
        entries: list[dict[str, Any]] = []
        for rank, model in enumerate(window_models, 1):
            positive = format_metric(model.get("positive"), digits=4)
            gtf_path = window_dir / (
                f"{rank:03d}_{model['genomic_start']}_{positive}_{safe_name(model['mrna_id'])}.gtf"
            )
            with gtf_path.open("w") as out:
                write_model_gtf(model, out)
            entries.append({
                "rank": rank,
                "model": model,
                "gtf": gtf_path,
                "title": track_title(model),
            })
        tracks_by_window[window["window_id"]] = entries
    return tracks_by_window


def write_nhmmer_gtf(rows: list[Feature], windows: list[Window], output_gtf: Path) -> int:
    output_gtf.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with output_gtf.open("w") as out:
        for row in rows:
            if not any(row_overlaps_window(row, window) for window in windows):
                continue
            attrs = row["attrs"]
            hit_id = attrs.get("ID") or f"nhmmer_{row['seqid']}_{row['start']}_{row['end']}_{n + 1}"
            gene_id = f"nhmmer_gene_{hit_id}"
            transcript_id = f"nhmmer_tx_{hit_id}"
            label = attrs.get("Name") or row["feature"]
            common_attrs = {
                "gene_id": gene_id,
                "transcript_id": transcript_id,
                "gene_name": label,
                "transcript_name": label,
                "hmm_feature": row["feature"],
                "evalue": attrs.get("E-value"),
            }
            write_gtf_line(
                out,
                row["seqid"],
                "nhmmer",
                "gene",
                row["start"],
                row["end"],
                row["score"],
                row["strand"],
                ".",
                common_attrs,
            )
            write_gtf_line(
                out,
                row["seqid"],
                "nhmmer",
                "transcript",
                row["start"],
                row["end"],
                row["score"],
                row["strand"],
                ".",
                common_attrs,
            )
            write_gtf_line(
                out,
                row["seqid"],
                "nhmmer",
                "exon",
                row["start"],
                row["end"],
                row["score"],
                row["strand"],
                ".",
                common_attrs,
            )
            n += 1
    return n


def write_miniprot_evidence_gtf(rows: list[Feature], windows: list[Window], output_gtf: Path) -> int:
    output_gtf.parent.mkdir(parents=True, exist_ok=True)
    mrnas: dict[str, Feature] = {}
    cds_by_parent: dict[str, list[Feature]] = {}

    for row in rows:
        if row["feature"] == "mRNA":
            mrna_id = row["attrs"].get("ID")
            if mrna_id:
                mrnas[mrna_id] = row
        elif row["feature"] == "CDS":
            parent = row["attrs"].get("Parent")
            if parent:
                cds_by_parent.setdefault(parent, []).append(row)

    n = 0
    with output_gtf.open("w") as out:
        for mrna_id, row in sorted(
            mrnas.items(),
            key=lambda item: (item[1]["seqid"], item[1]["start"], item[1]["end"], item[0]),
        ):
            if not any(row_overlaps_window(row, window) for window in windows):
                continue
            target = row["attrs"].get("Target", "").split()
            target_protein = target[0] if target else mrna_id
            common_attrs = {
                "gene_id": f"evidence_gene_{mrna_id}",
                "transcript_id": f"evidence_tx_{mrna_id}",
                "gene_name": target_protein,
                "transcript_name": target_protein,
                "identity": row["attrs"].get("Identity"),
                "positive": row["attrs"].get("Positive"),
                "rank": row["attrs"].get("Rank"),
            }
            write_gtf_line(
                out,
                row["seqid"],
                "miniprot_evidence",
                "gene",
                row["start"],
                row["end"],
                row["score"],
                row["strand"],
                ".",
                common_attrs,
            )
            write_gtf_line(
                out,
                row["seqid"],
                "miniprot_evidence",
                "transcript",
                row["start"],
                row["end"],
                row["score"],
                row["strand"],
                ".",
                common_attrs,
            )
            for i, cds in enumerate(sorted(cds_by_parent.get(mrna_id, []), key=lambda item: item["start"]), 1):
                block_attrs = dict(common_attrs)
                block_attrs["exon_number"] = i
                write_gtf_line(
                    out,
                    cds["seqid"],
                    "miniprot_evidence",
                    "exon",
                    cds["start"],
                    cds["end"],
                    cds["score"],
                    cds["strand"],
                    ".",
                    block_attrs,
                )
                write_gtf_line(
                    out,
                    cds["seqid"],
                    "miniprot_evidence",
                    "CDS",
                    cds["start"],
                    cds["end"],
                    cds["score"],
                    cds["strand"],
                    cds["phase"],
                    block_attrs,
                )
            n += 1
    return n


def collect_plot_region(
    window: Window,
    models_by_window: dict[str, list[dict[str, Any]]],
    nhmmer_rows: list[Feature],
    miniprot_evidence_rows: list[Feature],
    region_flank: int,
    clip_region_to_window: bool,
) -> dict[str, Any]:
    starts: list[int] = []
    ends: list[int] = []

    for model in models_by_window.get(window["window_id"], []):
        starts.append(model["genomic_start"])
        ends.append(model["genomic_end"])
        for block in model.get("cds_blocks", []):
            starts.append(block["genomic_start"])
            ends.append(block["genomic_end"])

    for row in nhmmer_rows:
        if row_overlaps_window(row, window):
            starts.append(row["start"])
            ends.append(row["end"])

    for row in miniprot_evidence_rows:
        if row["feature"] in PLOT_FEATURES and row_overlaps_window(row, window):
            starts.append(row["start"])
            ends.append(row["end"])

    if not starts:
        plot_start = window["start"]
        plot_end = window["end"]
    else:
        plot_start = min(starts) - region_flank
        plot_end = max(ends) + region_flank

    if clip_region_to_window:
        plot_start = max(plot_start, window["start"])
        plot_end = min(plot_end, window["end"])

    plot_start = max(1, plot_start)
    if plot_end < plot_start:
        plot_start = window["start"]
        plot_end = window["end"]

    return {
        "plot_start": plot_start,
        "plot_end": plot_end,
        "plot_region": f"{window['chrom']}:{plot_start}-{plot_end}",
        "full_region": window["region"],
        "region_flank": region_flank,
    }


def write_track_ini(
    ini_path: Path,
    ordered_window_tracks: list[dict[str, Any]],
    nhmmer_gtf: Path,
    miniprot_evidence_gtf: Path,
    per_mrna_track_height: float,
) -> None:
    ini_path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []

    for entry in ordered_window_tracks:
        model = entry["model"]
        section_name = f"window_miniprot_{entry['rank']:03d}_{safe_name(model['mrna_id'])}"
        lines.extend([
            f"[{section_name}]",
            f"file = {entry['gtf'].resolve()}",
            f"title = {entry['title']}",
            f"height = {per_mrna_track_height}",
            "merge_transcripts = false",
            "prefered_name = transcript_name",
            "labels = false",
            "max_labels = 1",
            "style = UCSC",
            "display = collapsed",
            "color = #1f78b4",
            "border_color = #0b4f71",
            "",
        ])

    lines.extend([
        "[nhmmer_evidence]",
        f"file = {nhmmer_gtf.resolve()}",
        "title = NHMMER NTD/CTD",
        "height = 2.5",
        "merge_transcripts = false",
        "prefered_name = gene_name",
        "labels = true",
        "max_labels = 80",
        "style = UCSC",
        "display = stacked",
        "gene_rows = 8",
        "color = #e31a1c",
        "border_color = #8f1111",
        "",
        "[miniprot_evidence_models]",
        f"file = {miniprot_evidence_gtf.resolve()}",
        "title = original miniprot evidence",
        "height = 6",
        "merge_transcripts = false",
        "prefered_name = transcript_name",
        "labels = true",
        "max_labels = 100",
        "style = UCSC",
        "display = stacked",
        "gene_rows = 16",
        "color = #33a02c",
        "border_color = #1b6b19",
        "",
        "[x-axis]",
        "where = bottom",
        "",
    ])
    ini_path.write_text("\n".join(lines))


def write_command_table(
    windows: list[Window],
    plot_regions: dict[str, dict[str, Any]],
    tracks_dir: Path,
    plots_dir: Path,
    output_tsv: Path,
    track_label_fraction: float,
    plot_font_size: int,
    plot_width: int,
    per_mrna_track_height: float,
) -> None:
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str | int | float]] = []
    for window in windows:
        ini = tracks_dir / f"{window['safe_name']}.ini"
        png = plots_dir / f"{window['safe_name']}.png"
        region_info = plot_regions[window["window_id"]]
        command = (
            "pyGenomeTracks "
            f"--tracks {shlex.quote(str(ini.resolve()))} "
            f"--region {shlex.quote(region_info['plot_region'])} "
            f"--trackLabelFraction {track_label_fraction} "
            f"--trackLabelHAlign left --fontSize {plot_font_size} "
            f"--width {plot_width} --dpi 130 "
            f"--outFileName {shlex.quote(str(png.resolve()))}"
        )
        rows.append({
            "window_id": window["window_id"],
            "full_region": region_info["full_region"],
            "plot_region": region_info["plot_region"],
            "plot_start": region_info["plot_start"],
            "plot_end": region_info["plot_end"],
            "region_flank": region_info["region_flank"],
            "track_label_fraction": track_label_fraction,
            "font_size": plot_font_size,
            "plot_width": plot_width,
            "per_mrna_track_height": per_mrna_track_height,
            "ini": str(ini),
            "output_png": str(png),
            "command": command,
        })

    fieldnames = [
        "window_id",
        "full_region",
        "plot_region",
        "plot_start",
        "plot_end",
        "region_flank",
        "track_label_fraction",
        "font_size",
        "plot_width",
        "per_mrna_track_height",
        "ini",
        "output_png",
        "command",
    ]
    with output_tsv.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


@app.command()
def main(
    models_jsonl: Path,
    nhmmer_gff: Path,
    miniprot_evidence_gff: Path,
    windows_bed: Path,
    output_dir: Path,
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
) -> None:
    pygt_dir = output_dir / "pygenometracks"
    gtf_dir = pygt_dir / "gtf"
    tracks_dir = pygt_dir / "tracks"
    plots_dir = pygt_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    models = load_models(models_jsonl, min_positive=min_positive)
    models_by_window = group_models_by_window(models)
    windows = load_windows(windows_bed)
    nhmmer_rows = parse_gff_rows(nhmmer_gff)
    miniprot_evidence_rows = parse_gff_rows(miniprot_evidence_gff)

    window_miniprot_gtf = gtf_dir / "window_miniprot_models.gtf"
    ordered_window_miniprot_dir = gtf_dir / "window_miniprot_ordered"
    nhmmer_gtf = gtf_dir / "nhmmer_evidence.gtf"
    miniprot_evidence_gtf = gtf_dir / "miniprot_evidence_models.gtf"
    commands_tsv = pygt_dir / "pygenometracks_commands.tsv"

    n_window = write_window_miniprot_gtf(models, window_miniprot_gtf)
    ordered_tracks_by_window = write_ordered_window_miniprot_gtfs(
        models_by_window,
        windows,
        ordered_window_miniprot_dir,
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

    for window in windows:
        write_track_ini(
            tracks_dir / f"{window['safe_name']}.ini",
            ordered_tracks_by_window.get(window["window_id"], []),
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

    typer.echo(f"window miniprot models: {n_window}")
    typer.echo(f"nhmmer evidence hits: {n_nhmmer}")
    typer.echo(f"miniprot evidence models: {n_evidence}")
    typer.echo(f"windows: {len(windows)}")
    typer.echo(f"region flank: {region_flank}")
    typer.echo(f"clip region to window: {clip_region_to_window}")
    typer.echo(f"track label fraction: {track_label_fraction}")
    typer.echo(f"plot font size: {plot_font_size}")
    typer.echo(f"plot width: {plot_width}")
    typer.echo(f"per-mRNA track height: {per_mrna_track_height}")
    typer.echo(f"GTF dir: {gtf_dir}")
    typer.echo(f"commands: {commands_tsv}")


if __name__ == "__main__":
    app()
