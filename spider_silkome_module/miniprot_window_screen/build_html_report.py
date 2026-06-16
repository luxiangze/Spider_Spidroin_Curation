from __future__ import annotations

import base64
from collections import Counter, defaultdict
import html
from pathlib import Path
import statistics
import struct
from typing import Any

from loguru import logger
import typer

from spider_silkome_module.config import INTERIM_DATA_DIR, PROCESSED_DATA_DIR, REPORTS_DIR
from spider_silkome_module.miniprot_window_screen.batch_common import read_tsv

app = typer.Typer(help="Build a self-contained HTML report for miniprot window screening.")


def esc(value: Any) -> str:
    return html.escape("" if value is None else str(value), quote=True)


def as_float(value: Any, default: float = 0.0) -> float:
    try:
        if value in ("", None, "."):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def as_int(value: Any, default: int = 0) -> int:
    try:
        if value in ("", None, "."):
            return default
        return int(float(value))
    except (TypeError, ValueError):
        return default


def fmt_int(value: Any) -> str:
    return f"{as_int(value):,}"


def fmt_float(value: Any, digits: int = 3) -> str:
    return f"{as_float(value):.{digits}f}"


def stat_summary(values: list[float]) -> dict[str, float]:
    if not values:
        return {"min": 0.0, "median": 0.0, "max": 0.0, "mean": 0.0}
    return {
        "min": min(values),
        "median": statistics.median(values),
        "max": max(values),
        "mean": statistics.mean(values),
    }


def read_png_dimensions(path: Path) -> tuple[int, int]:
    with path.open("rb") as handle:
        header = handle.read(24)
    if len(header) >= 24 and header[:8] == b"\x89PNG\r\n\x1a\n":
        width, height = struct.unpack(">II", header[16:24])
        return width, height
    return 1600, 700


def svg_bar_chart(
    data: list[tuple[str, float]],
    title: str,
    x_label: str = "",
    width: int = 920,
    row_height: int = 24,
    color: str = "#2f6f8f",
) -> str:
    if not data:
        return empty_svg(title, "No data")
    margin_left = 230
    margin_right = 80
    margin_top = 46
    margin_bottom = 44
    height = max(180, margin_top + margin_bottom + row_height * len(data))
    plot_width = width - margin_left - margin_right
    max_value = max(value for _, value in data) or 1
    parts = [
        f'<svg class="chart" viewBox="0 0 {width} {height}" role="img" '
        f'aria-label="{esc(title)}">',
        f'<text class="svg-title" x="{width / 2:.1f}" y="24" text-anchor="middle">{esc(title)}</text>',
        f'<line class="axis" x1="{margin_left}" y1="{height - margin_bottom}" '
        f'x2="{width - margin_right}" y2="{height - margin_bottom}"/>',
    ]
    for idx, (label, value) in enumerate(data):
        y = margin_top + idx * row_height
        bar_width = 0 if max_value == 0 else value / max_value * plot_width
        parts.extend([
            f'<text class="svg-label" x="{margin_left - 10}" y="{y + 15}" '
            f'text-anchor="end">{esc(label)}</text>',
            f'<rect x="{margin_left}" y="{y + 3}" width="{bar_width:.1f}" '
            f'height="{row_height - 8}" rx="3" fill="{color}"/>',
            f'<text class="svg-value" x="{margin_left + bar_width + 6:.1f}" '
            f'y="{y + 15}">{value:g}</text>',
        ])
    if x_label:
        parts.append(
            f'<text class="svg-note" x="{margin_left + plot_width / 2:.1f}" '
            f'y="{height - 12}" text-anchor="middle">{esc(x_label)}</text>'
        )
    parts.append("</svg>")
    return "\n".join(parts)


def svg_histogram(
    values: list[float],
    title: str,
    x_label: str,
    bins: int = 20,
    width: int = 920,
    height: int = 290,
    color: str = "#668f3f",
) -> str:
    if not values:
        return empty_svg(title, "No data")
    min_v = min(values)
    max_v = max(values)
    if min_v == max_v:
        counts = [len(values)]
    else:
        step = (max_v - min_v) / bins
        counts = [0 for _ in range(bins)]
        for value in values:
            idx = min(bins - 1, int((value - min_v) / step))
            counts[idx] += 1
    margin_left = 64
    margin_right = 24
    margin_top = 46
    margin_bottom = 56
    plot_width = width - margin_left - margin_right
    plot_height = height - margin_top - margin_bottom
    max_count = max(counts) or 1
    bar_gap = 3
    bar_width = (plot_width - bar_gap * (len(counts) - 1)) / len(counts)
    parts = [
        f'<svg class="chart" viewBox="0 0 {width} {height}" role="img" '
        f'aria-label="{esc(title)}">',
        f'<text class="svg-title" x="{width / 2:.1f}" y="24" text-anchor="middle">{esc(title)}</text>',
        f'<line class="axis" x1="{margin_left}" y1="{height - margin_bottom}" '
        f'x2="{width - margin_right}" y2="{height - margin_bottom}"/>',
        f'<line class="axis" x1="{margin_left}" y1="{margin_top}" x2="{margin_left}" '
        f'y2="{height - margin_bottom}"/>',
    ]
    for idx, count in enumerate(counts):
        x = margin_left + idx * (bar_width + bar_gap)
        bar_h = count / max_count * plot_height
        y = height - margin_bottom - bar_h
        parts.append(
            f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_width:.1f}" height="{bar_h:.1f}" '
            f'rx="2" fill="{color}"/>'
        )
    parts.extend([
        f'<text class="svg-note" x="{margin_left}" y="{height - 30}" text-anchor="middle">'
        f'{min_v:.2g}</text>',
        f'<text class="svg-note" x="{width - margin_right}" y="{height - 30}" text-anchor="end">'
        f'{max_v:.2g}</text>',
        f'<text class="svg-note" x="{width / 2:.1f}" y="{height - 10}" text-anchor="middle">'
        f'{esc(x_label)}</text>',
        f'<text class="svg-note" x="18" y="{margin_top + 12}" transform="rotate(-90 18 {margin_top + 12})">'
        "count</text>",
        "</svg>",
    ])
    return "\n".join(parts)


def svg_donut(data: list[tuple[str, float]], title: str, width: int = 520, height: int = 300) -> str:
    if not data:
        return empty_svg(title, "No data", width=width, height=height)
    colors = ["#2f6f8f", "#c7633d", "#668f3f", "#9a6fb0", "#d6a43a", "#4f7db8"]
    total = sum(value for _, value in data) or 1
    cx = 150
    cy = 158
    radius = 92
    stroke = 34
    circumference = 2 * 3.141592653589793 * radius
    offset = 0.0
    parts = [
        f'<svg class="chart" viewBox="0 0 {width} {height}" role="img" '
        f'aria-label="{esc(title)}">',
        f'<text class="svg-title" x="{width / 2:.1f}" y="24" text-anchor="middle">{esc(title)}</text>',
    ]
    for idx, (label, value) in enumerate(data):
        dash = value / total * circumference
        parts.append(
            f'<circle cx="{cx}" cy="{cy}" r="{radius}" fill="none" stroke="{colors[idx % len(colors)]}" '
            f'stroke-width="{stroke}" stroke-dasharray="{dash:.2f} {circumference - dash:.2f}" '
            f'stroke-dashoffset="{-offset:.2f}" transform="rotate(-90 {cx} {cy})"/>'
        )
        offset += dash
    parts.append(f'<text class="donut-total" x="{cx}" y="{cy - 4}" text-anchor="middle">{int(total):,}</text>')
    parts.append(f'<text class="svg-note" x="{cx}" y="{cy + 20}" text-anchor="middle">selected MPID</text>')
    for idx, (label, value) in enumerate(data):
        y = 94 + idx * 34
        parts.extend([
            f'<rect x="300" y="{y - 13}" width="16" height="16" rx="3" fill="{colors[idx % len(colors)]}"/>',
            f'<text class="svg-label" x="324" y="{y}">{esc(label)}: {int(value):,} '
            f'({value / total * 100:.1f}%)</text>',
        ])
    parts.append("</svg>")
    return "\n".join(parts)


def empty_svg(title: str, message: str, width: int = 920, height: int = 180) -> str:
    return "\n".join([
        f'<svg class="chart" viewBox="0 0 {width} {height}" role="img" aria-label="{esc(title)}">',
        f'<text class="svg-title" x="{width / 2:.1f}" y="28" text-anchor="middle">{esc(title)}</text>',
        f'<text class="svg-note" x="{width / 2:.1f}" y="{height / 2:.1f}" text-anchor="middle">{esc(message)}</text>',
        "</svg>",
    ])


def workflow_svg() -> str:
    steps = [
        ("Genome + references", "基因组与全长蛛丝蛋白参考序列"),
        ("Miniprot + HMMER", "初筛同源比对与 NTD/CTD 结构域证据"),
        ("Typing-guided locus", "用 typing 结果约束候选边界和基因座"),
        ("Candidate windows", "候选窗口提取并重比对全长蛋白"),
        ("Selected MPID", "筛选高可信 MPID、GFF3 与 FASTA"),
    ]
    width = 1180
    height = 260
    box_w = 190
    box_h = 92
    gap = 44
    x0 = 28
    y = 82
    parts = [
        f'<svg class="workflow" viewBox="0 0 {width} {height}" role="img" aria-label="workflow">',
        '<text class="workflow-title" x="590" y="34" text-anchor="middle">蛛丝蛋白候选序列鉴定流程</text>',
    ]
    for idx, (title, subtitle) in enumerate(steps):
        x = x0 + idx * (box_w + gap)
        color = ["#2f6f8f", "#4f7db8", "#668f3f", "#d6a43a", "#c7633d"][idx]
        parts.extend([
            f'<rect x="{x}" y="{y}" width="{box_w}" height="{box_h}" rx="12" fill="{color}" opacity="0.95"/>',
            f'<text class="workflow-step" x="{x + box_w / 2}" y="{y + 34}" text-anchor="middle">{esc(title)}</text>',
            f'<text class="workflow-sub" x="{x + box_w / 2}" y="{y + 61}" text-anchor="middle">{esc(subtitle)}</text>',
        ])
        if idx < len(steps) - 1:
            ax = x + box_w + 8
            ay = y + box_h / 2
            parts.extend([
                f'<line x1="{ax}" y1="{ay}" x2="{ax + gap - 18}" y2="{ay}" stroke="#3b4752" stroke-width="3"/>',
                f'<polygon points="{ax + gap - 18},{ay - 7} {ax + gap - 18},{ay + 7} '
                f'{ax + gap - 6},{ay}" fill="#3b4752"/>',
            ])
    parts.append("</svg>")
    return "\n".join(parts)


def table_html(rows: list[dict[str, Any]], columns: list[tuple[str, str]], limit: int | None = None) -> str:
    display_rows = rows[:limit] if limit is not None else rows
    head = "".join(f"<th>{esc(label)}</th>" for _, label in columns)
    body_parts = []
    for row in display_rows:
        body_parts.append("<tr>" + "".join(f"<td>{esc(row.get(key, ''))}</td>" for key, _ in columns) + "</tr>")
    if not body_parts:
        body_parts.append(f'<tr><td colspan="{len(columns)}" class="muted">No records</td></tr>')
    return f"<table><thead><tr>{head}</tr></thead><tbody>{''.join(body_parts)}</tbody></table>"


def summarize_step(path: Path) -> dict[str, Any]:
    rows = read_tsv(path)
    status_counter = Counter(row.get("status", "") or "unknown" for row in rows)
    return {
        "step": path.stem.replace("_summary", ""),
        "rows": len(rows),
        "ok": status_counter.get("ok", 0),
        "failed": status_counter.get("failed", 0),
        "skipped": status_counter.get("skipped", 0),
        "other": len(rows) - status_counter.get("ok", 0) - status_counter.get("failed", 0) - status_counter.get("skipped", 0),
    }


def selected_plot_key(path: Path, selected_window_counts: dict[tuple[str, str], int]) -> tuple[int, str, str]:
    species = path.parents[2].name
    window = f"{path.stem}:."
    return (-selected_window_counts.get((species, window), 0), species, path.name)


def image_svg(path: Path, title: str, caption: str) -> str:
    width, height = read_png_dimensions(path)
    data = base64.b64encode(path.read_bytes()).decode("ascii")
    view_height = height + 72
    return "\n".join([
        f'<figure class="track-figure"><svg class="track-svg" viewBox="0 0 {width} {view_height}" '
        f'role="img" aria-label="{esc(title)}">',
        f'<text class="track-title" x="18" y="28">{esc(title)}</text>',
        f'<text class="track-caption" x="18" y="52">{esc(caption)}</text>',
        f'<image x="0" y="72" width="{width}" height="{height}" '
        f'href="data:image/png;base64,{data}"/>',
        "</svg></figure>",
    ])


def locus_svg_from_gff(gff_path: Path, title: str) -> str:
    mrnas: list[dict[str, Any]] = []
    cds_by_parent: dict[str, list[tuple[int, int]]] = defaultdict(list)
    if not gff_path.exists():
        return empty_svg(title, "GFF3 not found")
    with gff_path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            feature = fields[2]
            start = as_int(fields[3])
            end = as_int(fields[4])
            attrs = parse_gff_attrs(fields[8])
            if feature == "mRNA":
                mrnas.append({
                    "id": attrs.get("ID", ""),
                    "seqid": fields[0],
                    "start": start,
                    "end": end,
                    "strand": fields[6],
                    "positive": attrs.get("positive", attrs.get("Positive", "")),
                    "type": attrs.get("typing_spidroin_type", ""),
                    "aa": attrs.get("protein_sequence_length", ""),
                })
            elif feature == "CDS" and attrs.get("Parent"):
                cds_by_parent[attrs["Parent"]].append((start, end))
    if not mrnas:
        return empty_svg(title, "No selected models")
    mrnas = sorted(mrnas, key=lambda row: (-as_float(row["positive"]), row["seqid"], row["start"]))[:8]
    min_start = min(row["start"] for row in mrnas)
    max_end = max(row["end"] for row in mrnas)
    for row in mrnas:
        for start, end in cds_by_parent.get(row["id"], []):
            min_start = min(min_start, start)
            max_end = max(max_end, end)
    span = max(1, max_end - min_start + 1)
    width = 980
    margin_left = 170
    margin_right = 40
    row_height = 40
    height = 70 + len(mrnas) * row_height
    plot_width = width - margin_left - margin_right
    parts = [
        f'<svg class="chart" viewBox="0 0 {width} {height}" role="img" aria-label="{esc(title)}">',
        f'<text class="svg-title" x="{width / 2:.1f}" y="24" text-anchor="middle">{esc(title)}</text>',
        f'<text class="svg-note" x="{margin_left}" y="48">{esc(mrnas[0]["seqid"])}:{min_start:,}-{max_end:,}</text>',
    ]
    for idx, row in enumerate(mrnas):
        y = 72 + idx * row_height
        tx = margin_left + (row["start"] - min_start) / span * plot_width
        tw = max(2, (row["end"] - row["start"] + 1) / span * plot_width)
        label = f'{row["id"]} {row["type"]} P={fmt_float(row["positive"])} aa={row["aa"]}'
        parts.extend([
            f'<text class="svg-label" x="{margin_left - 10}" y="{y + 5}" text-anchor="end">{esc(label)}</text>',
            f'<line x1="{tx:.1f}" y1="{y}" x2="{tx + tw:.1f}" y2="{y}" stroke="#3b4752" stroke-width="2"/>',
        ])
        for start, end in cds_by_parent.get(row["id"], []):
            x = margin_left + (start - min_start) / span * plot_width
            w = max(2, (end - start + 1) / span * plot_width)
            parts.append(f'<rect x="{x:.1f}" y="{y - 8}" width="{w:.1f}" height="16" rx="2" fill="#2f6f8f"/>')
    parts.append("</svg>")
    return "\n".join(parts)


def parse_gff_attrs(attr_text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in attr_text.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value
    return attrs


def css() -> str:
    return """
:root {
  --ink: #1f2933;
  --muted: #61707c;
  --line: #d9e1e7;
  --bg: #f6f8fa;
  --panel: #ffffff;
  --blue: #2f6f8f;
  --green: #668f3f;
  --orange: #c7633d;
}
* { box-sizing: border-box; }
body {
  margin: 0;
  background: var(--bg);
  color: var(--ink);
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", "Noto Sans CJK SC", "Microsoft YaHei", Arial, sans-serif;
  line-height: 1.55;
}
header {
  background: #1f2933;
  color: white;
  padding: 34px 48px 30px;
}
header h1 { margin: 0 0 8px; font-size: 30px; font-weight: 760; }
header p { margin: 0; color: #d7e0e7; }
main { max-width: 1240px; margin: 0 auto; padding: 28px 28px 56px; }
section {
  background: var(--panel);
  border: 1px solid var(--line);
  border-radius: 10px;
  padding: 22px;
  margin-bottom: 22px;
  box-shadow: 0 1px 2px rgba(31, 41, 51, 0.04);
}
h2 { margin: 0 0 14px; font-size: 22px; }
h3 { margin: 18px 0 10px; font-size: 17px; }
.muted { color: var(--muted); }
.grid { display: grid; gap: 14px; }
.grid.metrics { grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); }
.metric {
  border: 1px solid var(--line);
  border-radius: 8px;
  padding: 14px;
  background: #fbfcfd;
}
.metric .value { font-size: 28px; font-weight: 780; color: var(--blue); }
.metric .label { color: var(--muted); font-size: 13px; margin-top: 3px; }
.charts { display: grid; grid-template-columns: repeat(auto-fit, minmax(460px, 1fr)); gap: 18px; align-items: start; }
.chart, .workflow, .track-svg { width: 100%; height: auto; display: block; }
.svg-title, .workflow-title { font-size: 18px; font-weight: 760; fill: var(--ink); }
.svg-label { font-size: 12px; fill: var(--ink); }
.svg-value { font-size: 12px; fill: var(--muted); }
.svg-note, .track-caption { font-size: 12px; fill: var(--muted); }
.donut-total { font-size: 28px; font-weight: 800; fill: var(--ink); }
.axis { stroke: #aab7c1; stroke-width: 1; }
.workflow-step { fill: #fff; font-size: 16px; font-weight: 740; }
.workflow-sub { fill: #eef5f8; font-size: 11px; }
table {
  width: 100%;
  border-collapse: collapse;
  font-size: 13px;
  margin-top: 10px;
}
th, td {
  border-bottom: 1px solid var(--line);
  padding: 7px 8px;
  text-align: left;
  vertical-align: top;
}
th {
  background: #eef3f6;
  font-weight: 700;
  position: sticky;
  top: 0;
}
.table-wrap { overflow-x: auto; border: 1px solid var(--line); border-radius: 8px; }
.track-grid { display: grid; grid-template-columns: 1fr; gap: 18px; }
.track-figure {
  margin: 0;
  padding: 12px;
  border: 1px solid var(--line);
  border-radius: 8px;
  background: #fbfcfd;
}
.track-title { font-size: 20px; font-weight: 760; fill: var(--ink); }
code { background: #eef3f6; padding: 2px 5px; border-radius: 4px; }
footer { color: var(--muted); font-size: 12px; margin-top: 24px; }
"""


def build_report(
    processed_root: Path,
    interim_root: Path,
    output: Path,
    max_track_plots: int,
    top_n_species: int,
    top_n_examples: int,
) -> None:
    selected_path = processed_root / "all_species_selected_spidroin_mpid.tsv"
    proteins_path = processed_root / "all_species_selected_spidroin_proteins.faa"
    manifest_path = processed_root / "species_run_manifest.tsv"
    selected_rows = read_tsv(selected_path)
    manifest_rows = read_tsv(manifest_path)
    task_name = processed_root.name

    if not selected_rows:
        raise FileNotFoundError(f"No selected rows found: {selected_path}")
    if not manifest_rows:
        raise FileNotFoundError(f"No manifest rows found: {manifest_path}")

    species_with_selected = sum(1 for row in manifest_rows if as_int(row.get("n_selected")) > 0)
    ok_species = sum(1 for row in manifest_rows if row.get("status") == "ok")
    status_counts = Counter(row.get("selection_status", "") for row in selected_rows)
    mode_counts = Counter(row.get("selection_mode", "") for row in selected_rows)
    internal_stop_candidates = sum(1 for row in selected_rows if as_int(row.get("internal_stop_count")) > 0)
    terminal_stop_only_candidates = sum(
        1 for row in selected_rows if str(row.get("terminal_stop_only", "")).lower() == "true"
    )
    type_counts = Counter(row.get("typing_spidroin_type", "") or "unknown" for row in selected_rows)
    species_counts = Counter(row.get("species", "") for row in selected_rows)

    positive_values = [as_float(row.get("positive")) for row in selected_rows if row.get("positive")]
    coverage_values = [as_float(row.get("query_coverage")) for row in selected_rows if row.get("query_coverage")]
    aa_values = [as_float(row.get("protein_sequence_length")) for row in selected_rows if row.get("protein_sequence_length")]
    exon_values = [as_float(row.get("exon_count")) for row in selected_rows if row.get("exon_count")]

    selected_window_counts: dict[tuple[str, str], int] = Counter(
        (row.get("species", ""), row.get("window_id", "")) for row in selected_rows
    )
    plot_paths = sorted(
        interim_root.glob("*/pygenometracks/plots/*.png"),
        key=lambda path: selected_plot_key(path, selected_window_counts),
    )
    embedded_plots = plot_paths[:max_track_plots]

    step_summaries = []
    for step in [
        "prepare_windows_summary.tsv",
        "run_miniprot_summary.tsv",
        "split_miniprot_summary.tsv",
        "parse_models_summary.tsv",
        "select_candidates_summary.tsv",
        "export_selected_summary.tsv",
        "build_tracks_summary.tsv",
        "pygenometracks_summary.tsv",
    ]:
        path = interim_root / "manifests" / step
        if path.exists():
            step_summaries.append(summarize_step(path))
        else:
            step_summaries.append({
                "step": step.replace("_summary.tsv", ""),
                "rows": 0,
                "ok": 0,
                "failed": 0,
                "skipped": 0,
                "other": "missing",
            })

    representative_rows = sorted(
        selected_rows,
        key=lambda row: (
            -as_float(row.get("positive")),
            -as_float(row.get("query_coverage")),
            -as_float(row.get("protein_sequence_length")),
            row.get("species", ""),
            row.get("mpid", ""),
        ),
    )[:top_n_examples]

    manifest_display = sorted(
        manifest_rows,
        key=lambda row: (-as_int(row.get("n_selected")), row.get("species", "")),
    )[:top_n_species]

    selected_display = sorted(
        selected_rows,
        key=lambda row: (row.get("species", ""), row.get("typing_spidroin_type", ""), row.get("mpid", "")),
    )[:top_n_examples]

    species_chart = svg_bar_chart(
        [(species, count) for species, count in species_counts.most_common(top_n_species)],
        f"Selected MPID 数量最多的 {top_n_species} 个物种",
        "selected MPID count",
        color="#2f6f8f",
    )
    type_chart = svg_bar_chart(
        [(spidroin_type, count) for spidroin_type, count in type_counts.most_common(16)],
        "蛛丝蛋白类型分布",
        "selected MPID count",
        color="#668f3f",
    )
    selection_chart = svg_donut(
        [(label or "unknown", value) for label, value in status_counts.most_common()],
        "筛选状态比例",
    )

    locus_gff = processed_root / "001.Allagelena_difficilis" / "selected_spidroin_models.gff3"
    track_svgs = []
    for path in embedded_plots:
        species = path.parents[2].name
        window = f"{path.stem}:."
        n_selected = selected_window_counts.get((species, window), 0)
        title = f"{species} | {window}"
        caption = f"pyGenomeTracks; selected MPID in this window: {n_selected}; source: {path}"
        track_svgs.append(image_svg(path, title, caption))

    track_message = (
        f"已发现 {len(plot_paths)} 张 pyGenomeTracks PNG，本报告嵌入 {len(embedded_plots)} 张。"
        if plot_paths
        else "当前还没有发现 pyGenomeTracks PNG；可在绘图步骤完成后重新生成本报告。"
    )

    pos_stats = stat_summary(positive_values)
    cov_stats = stat_summary(coverage_values)
    aa_stats = stat_summary(aa_values)
    exon_stats = stat_summary(exon_values)

    html_text = f"""<!doctype html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Miniprot Window Spidroin Screen Report - {esc(task_name)}</title>
<style>{css()}</style>
</head>
<body>
<header>
  <h1>蛛丝蛋白候选序列鉴定报告</h1>
  <p>Task: <code>{esc(task_name)}</code> | Processed: <code>{esc(processed_root)}</code></p>
</header>
<main>
  <section>
    <h2>执行摘要</h2>
    <div class="grid metrics">
      <div class="metric"><div class="value">{len(manifest_rows):,}</div><div class="label">纳入物种</div></div>
      <div class="metric"><div class="value">{ok_species:,}</div><div class="label">运行成功物种</div></div>
      <div class="metric"><div class="value">{species_with_selected:,}</div><div class="label">有 selected MPID 的物种</div></div>
      <div class="metric"><div class="value">{len(selected_rows):,}</div><div class="label">selected MPID 总数</div></div>
      <div class="metric"><div class="value">{status_counts.get("selected_strong", 0):,}</div><div class="label">strong candidates</div></div>
      <div class="metric"><div class="value">{status_counts.get("selected_rescue", 0):,}</div><div class="label">rescue candidates</div></div>
      <div class="metric"><div class="value">{internal_stop_candidates:,}</div><div class="label">内部 stop 候选</div></div>
      <div class="metric"><div class="value">{terminal_stop_only_candidates:,}</div><div class="label">仅末端 stop 候选</div></div>
      <div class="metric"><div class="value">{sum(1 for row in manifest_rows if row.get("typing_status") == "present"):,}</div><div class="label">typing 结果存在物种</div></div>
      <div class="metric"><div class="value">{sum(1 for row in manifest_rows if row.get("nhmmer_status") == "present"):,}</div><div class="label">NHMMER 结果存在物种</div></div>
    </div>
    <p class="muted">筛选模式：{esc(", ".join(f"{k}={v}" for k, v in mode_counts.items()))}。当前报告只读取已有结果，不重新运行 miniprot、HMMER 或 pyGenomeTracks。</p>
  </section>

  <section>
    <h2>方法概览</h2>
    {workflow_svg()}
  </section>

  <section>
    <h2>流程 QC</h2>
    <div class="table-wrap">
    {table_html(step_summaries, [
        ("step", "step"),
        ("rows", "rows"),
        ("ok", "ok"),
        ("failed", "failed"),
        ("skipped", "skipped"),
        ("other", "other"),
    ])}
    </div>
  </section>

  <section>
    <h2>总体结果图表</h2>
    <div class="charts">
      {selection_chart}
      {species_chart}
      {type_chart}
      {svg_histogram(positive_values, "Positive 分布", "Positive", color="#2f6f8f")}
      {svg_histogram(coverage_values, "Query coverage 分布", "query coverage", color="#668f3f")}
      {svg_histogram(aa_values, "蛋白长度分布", "amino acid length", color="#c7633d")}
      {svg_histogram(exon_values, "外显子数量分布", "exon count", color="#9a6fb0")}
    </div>
    <p class="muted">
      Positive: min={pos_stats["min"]:.3f}, median={pos_stats["median"]:.3f}, max={pos_stats["max"]:.3f};
      coverage: min={cov_stats["min"]:.3f}, median={cov_stats["median"]:.3f}, max={cov_stats["max"]:.3f};
      aa length: min={aa_stats["min"]:.0f}, median={aa_stats["median"]:.0f}, max={aa_stats["max"]:.0f};
      exon count: min={exon_stats["min"]:.0f}, median={exon_stats["median"]:.0f}, max={exon_stats["max"]:.0f}.
    </p>
  </section>

  <section>
    <h2>代表性 selected MPID</h2>
    <div class="table-wrap">
    {table_html(representative_rows, [
        ("species", "species"),
        ("window_id", "window"),
        ("mpid", "MPID"),
        ("selection_status", "status"),
        ("typing_spidroin_type", "type"),
        ("positive", "positive"),
        ("query_coverage", "coverage"),
        ("protein_sequence_length", "aa"),
        ("terminal_stop_only", "terminal stop only"),
        ("internal_stop_count", "internal stops"),
        ("exon_count", "exons"),
        ("target_protein", "target protein"),
    ])}
    </div>
  </section>

  <section>
    <h2>物种汇总 Top {top_n_species}</h2>
    <div class="table-wrap">
    {table_html(manifest_display, [
        ("species", "species"),
        ("status", "status"),
        ("nhmmer_status", "NHMMER"),
        ("typing_status", "typing"),
        ("typing_locus_count", "typing loci"),
        ("typing_full_length_count", "full length loci"),
        ("selection_mode", "selection mode"),
        ("n_selected", "selected"),
        ("n_selected_proteins", "proteins"),
        ("n_plots", "plots in manifest"),
    ])}
    </div>
  </section>

  <section>
    <h2>Selected MPID 示例表</h2>
    <div class="table-wrap">
    {table_html(selected_display, [
        ("species", "species"),
        ("window_id", "window"),
        ("mpid", "MPID"),
        ("selection_status", "status"),
        ("typing_spidroin_id", "typing ID"),
        ("typing_spidroin_type", "type"),
        ("positive", "positive"),
        ("identity", "identity"),
        ("query_coverage", "coverage"),
        ("protein_sequence_length", "aa"),
        ("terminal_stop_only", "terminal stop only"),
        ("internal_stop_count", "internal stops"),
        ("exon_count", "exons"),
        ("intron_count", "introns"),
    ])}
    </div>
  </section>

  <section>
    <h2>001 物种 selected gene model 摘要图</h2>
    <p class="muted">从 selected GFF3 抽取 Positive 最高的 mRNA/CDS，绘制简化结构示意；完整坐标与属性请查看 GFF3。</p>
    {locus_svg_from_gff(locus_gff, "001.Allagelena_difficilis selected models")}
  </section>

  <section>
    <h2>已生成基因组轨道图</h2>
    <p class="muted">{esc(track_message)}</p>
    <div class="track-grid">
      {"".join(track_svgs) if track_svgs else '<p class="muted">暂无可嵌入轨道图。</p>'}
    </div>
  </section>

  <section>
    <h2>数据产物</h2>
    <div class="table-wrap">
    {table_html([
        {"name": "All selected MPID TSV", "path": selected_path},
        {"name": "All selected proteins FASTA", "path": proteins_path},
        {"name": "Species run manifest", "path": manifest_path},
        {"name": "Interim root", "path": interim_root},
    ], [("name", "name"), ("path", "path")])}
    </div>
  </section>

  <footer>
    Generated by <code>spider_silkome_module.miniprot_window_screen.build_html_report</code>.
  </footer>
</main>
</body>
</html>
"""
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(html_text)
    logger.success(
        f"Report written: {output} selected={len(selected_rows)} species={len(manifest_rows)} "
        f"embedded_tracks={len(embedded_plots)}/{len(plot_paths)}"
    )


@app.command()
def main(
    processed_root: Path = typer.Option(
        PROCESSED_DATA_DIR / "miniprot_window_spidroin_screen_20260526_115407",
        "--processed-root",
    ),
    interim_root: Path = typer.Option(
        INTERIM_DATA_DIR / "miniprot_window_spidroin_screen_20260526_115407",
        "--interim-root",
    ),
    output: Path = typer.Option(
        REPORTS_DIR / "miniprot_window_spidroin_screen_20260526_115407_report.html",
        "--output",
    ),
    max_track_plots: int = typer.Option(24, "--max-track-plots"),
    top_n_species: int = typer.Option(25, "--top-n-species"),
    top_n_examples: int = typer.Option(20, "--top-n-examples"),
) -> None:
    build_report(
        processed_root=processed_root,
        interim_root=interim_root,
        output=output,
        max_track_plots=max_track_plots,
        top_n_species=top_n_species,
        top_n_examples=top_n_examples,
    )


if __name__ == "__main__":
    app()
