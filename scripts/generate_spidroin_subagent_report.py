#!/usr/bin/env python3
"""Render the per-sequence subagent spidroin review report.

The subagent review JSONL is intentionally separate from deterministic feature
extraction. Each JSONL line must contain one reviewed FASTA record keyed by
``record_id``. By default the renderer fails unless every FASTA record has a
review, so missing subagent results are visible immediately.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime
import html
import json
import math
from pathlib import Path
import re
from statistics import mean, median
from typing import Any, Iterable

DEFAULT_FASTA = Path(
    "data/processed/miniprot_window_spidroin_screen_20260526_115407/"
    "all_species_selected_spidroin_proteins.faa"
)
DEFAULT_REVIEW_JSONL = Path(
    "reports/miniprot_window_spidroin_screen_20260526_115407_"
    "subagent_spidroin_review.jsonl"
)
DEFAULT_HTML = Path(
    "reports/miniprot_window_spidroin_screen_20260526_115407_"
    "subagent_spidroin_review.html"
)

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY*XUOBZJ")
MOTIF_PATTERNS: dict[str, str] = {
    "(GA)n": r"(?:GA){2,}",
    "(GS)n": r"(?:GS){2,}",
    "polyA": r"A{4,}",
    "polyS": r"S{4,}",
    "polyQ": r"Q{3,}",
    "GGX": r"GG[A-Z]",
    "GPGXX": r"GPG[A-Z]{2}",
    "QQ": r"Q{2,}",
}
MOTIF_COLORS: dict[str, str] = {
    "(GA)n": "#2f9c95",
    "(GS)n": "#7a77c8",
    "polyA": "#d96c4f",
    "polyS": "#58a05c",
    "polyQ": "#c59a2e",
    "GGX": "#4d84c4",
    "GPGXX": "#cf6fa8",
    "QQ": "#8f7a3d",
}
TYPE_COLORS = [
    "#2f9c95",
    "#d96c4f",
    "#4d84c4",
    "#58a05c",
    "#c59a2e",
    "#7a77c8",
    "#cf6fa8",
    "#8b6f47",
    "#5d8a72",
    "#b85d75",
    "#697a21",
    "#8f6bb2",
]


@dataclass(frozen=True)
class ProteinRecord:
    record_id: str
    header: str
    sequence: str
    species: str
    mpid: str
    declared_type: str
    aa_len_header: int | None
    positive: float | None
    query_coverage: float | None
    terminal_stop_only: bool | None
    internal_stop_count: int | None
    exon_count: int | None
    selection_mode: str
    typing_spidroin_id: str


def parse_bool(value: str | None) -> bool | None:
    if value is None:
        return None
    if value.lower() == "true":
        return True
    if value.lower() == "false":
        return False
    return None


def parse_int(value: str | None) -> int | None:
    if value in (None, ""):
        return None
    try:
        return int(float(value))
    except ValueError:
        return None


def parse_float(value: str | None) -> float | None:
    if value in (None, ""):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def parse_header(header: str, sequence: str) -> ProteinRecord:
    fields = dict(re.findall(r"([A-Za-z_]+)=([^\s]+)", header))
    first = header.split()[0]
    first_parts = first.split("|")
    species = first_parts[0] if first_parts else "unknown"
    mpid = first_parts[-1] if first_parts else first
    return ProteinRecord(
        record_id=first,
        header=header,
        sequence=sequence,
        species=species,
        mpid=mpid,
        declared_type=fields.get("typing_spidroin_type", "unknown"),
        aa_len_header=parse_int(fields.get("aa_len")),
        positive=parse_float(fields.get("positive")),
        query_coverage=parse_float(fields.get("query_coverage")),
        terminal_stop_only=parse_bool(fields.get("terminal_stop_only")),
        internal_stop_count=parse_int(fields.get("internal_stop_count")),
        exon_count=parse_int(fields.get("exon_count")),
        selection_mode=fields.get("selection_mode", "unknown"),
        typing_spidroin_id=fields.get("typing_spidroin_id", "unknown"),
    )


def read_fasta(path: Path) -> list[ProteinRecord]:
    records: list[ProteinRecord] = []
    header: str | None = None
    seq_parts: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    sequence = "".join(seq_parts).upper()
                    records.append(parse_header(header, sequence))
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        sequence = "".join(seq_parts).upper()
        records.append(parse_header(header, sequence))
    return records


def union_coverage(spans: Iterable[tuple[int, int]], length: int) -> float:
    if length <= 0:
        return 0.0
    merged: list[list[int]] = []
    for start, end in sorted(spans):
        if start >= end:
            continue
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    covered = sum(end - start for start, end in merged)
    return covered / length


def shannon_entropy(text: str) -> float:
    if not text:
        return 0.0
    counts = Counter(text)
    total = len(text)
    return -sum((count / total) * math.log2(count / total) for count in counts.values())


def low_complexity_fraction(seq: str, window: int = 50, step: int = 25) -> float:
    if not seq:
        return 0.0
    if len(seq) <= window:
        return 1.0 if shannon_entropy(seq) <= 2.25 else 0.0
    covered = [False] * len(seq)
    for start in range(0, len(seq) - window + 1, step):
        end = start + window
        frag = seq[start:end]
        most_common = Counter(frag).most_common(1)[0][1] / window
        if shannon_entropy(frag) <= 2.25 or most_common >= 0.36:
            for index in range(start, end):
                covered[index] = True
    return sum(covered) / len(seq)


def longest_run(seq: str, residue: str) -> int:
    best = current = 0
    for char in seq:
        if char == residue:
            current += 1
            best = max(best, current)
        else:
            current = 0
    return best


def calculate_features(record: ProteinRecord) -> dict[str, Any]:
    seq = record.sequence
    length = len(seq)
    aa_counts = Counter(seq)
    composition = {
        aa: aa_counts.get(aa, 0) / length if length else 0.0 for aa in ["A", "G", "S", "Q", "P"]
    }
    motif_counts: dict[str, int] = {}
    motif_spans: list[tuple[int, int]] = []
    longest_motif = 0
    for name, pattern in MOTIF_PATTERNS.items():
        matches = list(re.finditer(pattern, seq))
        motif_counts[name] = len(matches)
        for match in matches:
            motif_spans.append(match.span())
            longest_motif = max(longest_motif, match.end() - match.start())
    invalid_residues = sorted(set(seq) - VALID_AA)
    stop_count = aa_counts.get("*", 0)
    internal_stop_count = record.internal_stop_count if record.internal_stop_count is not None else 0
    motif_coverage = union_coverage(motif_spans, length)
    low_complexity = low_complexity_fraction(seq)
    ags_fraction = composition["A"] + composition["G"] + composition["S"]
    score = deterministic_score(
        length=length,
        ags_fraction=ags_fraction,
        motif_coverage=motif_coverage,
        low_complexity=low_complexity,
        query_coverage=record.query_coverage,
        positive=record.positive,
        internal_stop_count=internal_stop_count,
        invalid_residue_count=len(invalid_residues),
        motif_counts=motif_counts,
    )
    return {
        "length": length,
        "aa_len_header": record.aa_len_header,
        "composition": composition,
        "ags_fraction": ags_fraction,
        "motif_counts": motif_counts,
        "motif_coverage": motif_coverage,
        "low_complexity_fraction": low_complexity,
        "longest_motif_span": longest_motif,
        "longest_poly_a": longest_run(seq, "A"),
        "longest_poly_s": longest_run(seq, "S"),
        "longest_poly_q": longest_run(seq, "Q"),
        "stop_count": stop_count,
        "internal_stop_count": internal_stop_count,
        "invalid_residues": invalid_residues,
        "deterministic_score": score,
        "deterministic_verdict": score_to_verdict(score, internal_stop_count),
        "feature_type_hint": infer_type_hint(record.declared_type, motif_counts, composition),
    }


def deterministic_score(
    *,
    length: int,
    ags_fraction: float,
    motif_coverage: float,
    low_complexity: float,
    query_coverage: float | None,
    positive: float | None,
    internal_stop_count: int,
    invalid_residue_count: int,
    motif_counts: dict[str, int],
) -> float:
    score = 0.0
    if length >= 2500:
        score += 22
    elif length >= 1500:
        score += 17
    elif length >= 1000:
        score += 11
    else:
        score += 4

    if ags_fraction >= 0.48:
        score += 19
    elif ags_fraction >= 0.40:
        score += 14
    elif ags_fraction >= 0.32:
        score += 8

    if motif_coverage >= 0.30:
        score += 18
    elif motif_coverage >= 0.16:
        score += 13
    elif motif_coverage >= 0.06:
        score += 7

    if low_complexity >= 0.35:
        score += 15
    elif low_complexity >= 0.20:
        score += 10
    elif low_complexity >= 0.08:
        score += 5

    motif_total = sum(motif_counts.values())
    if motif_total >= 80:
        score += 8
    elif motif_total >= 30:
        score += 5
    elif motif_total >= 8:
        score += 3

    if query_coverage is not None:
        if query_coverage >= 0.95:
            score += 7
        elif query_coverage >= 0.80:
            score += 4
    if positive is not None:
        if positive >= 0.68:
            score += 7
        elif positive >= 0.58:
            score += 4
        elif positive >= 0.50:
            score += 2

    if internal_stop_count == 0:
        score += 4
    else:
        score -= min(25, 10 + internal_stop_count * 5)
    if invalid_residue_count:
        score -= 10
    if length < 900:
        score -= 10
    if motif_coverage < 0.03 and low_complexity < 0.08:
        score -= 8
    return round(max(0.0, min(100.0, score)), 1)


def score_to_verdict(score: float, internal_stop_count: int) -> str:
    if internal_stop_count > 0:
        return "needs_review"
    if score >= 70:
        return "strong_spidroin_like"
    if score >= 55:
        return "moderate_spidroin_like"
    if score >= 40:
        return "weak_or_partial"
    return "unlikely_or_fragmentary"


def infer_type_hint(
    declared_type: str,
    motif_counts: dict[str, int],
    composition: dict[str, float],
) -> str:
    declared = declared_type.lower()
    gly_ala = composition.get("G", 0.0) + composition.get("A", 0.0)
    ser = composition.get("S", 0.0)
    if "/" in declared_type:
        return "mixed_type_header"
    if "flag" in declared:
        return "flag_like" if motif_counts.get("GPGXX", 0) >= 5 else "flag_motif_sparse"
    if "masp" in declared:
        return "major_ampullate_like" if gly_ala >= 0.42 else "major_ampullate_low_AG"
    if "misp" in declared:
        return "minor_ampullate_like" if gly_ala + ser >= 0.48 else "minor_ampullate_low_AGS"
    if "acsp" in declared:
        return "aciniform_like" if ser + composition.get("A", 0.0) >= 0.28 else "aciniform_signal_sparse"
    if "cysp" in declared or "pysp" in declared:
        return "tubuliform_pyriform_like" if gly_ala + ser >= 0.38 else "type_signal_sparse"
    return "untyped"


def deterministic_agent_seed(record: ProteinRecord, features: dict[str, Any]) -> dict[str, Any]:
    matched: list[str] = []
    concerns: list[str] = []
    if features["length"] >= 1500:
        matched.append("长蛋白，符合典型 spidroin 大尺寸特征")
    if features["ags_fraction"] >= 0.40:
        matched.append("A/G/S 富集")
    if features["motif_coverage"] >= 0.10:
        matched.append("检测到丝蛋白常见重复 motif")
    if features["low_complexity_fraction"] >= 0.20:
        matched.append("低复杂度重复区明显")
    if record.query_coverage is not None and record.query_coverage >= 0.90:
        matched.append("miniprot query coverage 高")
    if record.positive is not None and record.positive >= 0.58:
        matched.append("同源比对 positive 支持")
    if features["internal_stop_count"]:
        concerns.append(f"存在 {features['internal_stop_count']} 个 internal stop")
    if features["motif_coverage"] < 0.05:
        concerns.append("典型重复 motif 覆盖较低")
    if features["low_complexity_fraction"] < 0.08:
        concerns.append("低复杂度特征较弱")
    if features["length"] < 1200:
        concerns.append("长度偏短，可能是片段")
    score = features["deterministic_score"]
    confidence = max(0.35, min(0.98, score / 100))
    return {
        "record_id": record.record_id,
        "species": record.species,
        "declared_type": record.declared_type,
        "is_spidroin_like": score >= 55 and not features["internal_stop_count"],
        "confidence": round(confidence, 2),
        "matched_features": matched,
        "concerns": concerns,
        "type_consistency": features["feature_type_hint"],
        "short_comment": (
            "程序特征摘要："
            + ("; ".join(matched[:3]) if matched else "蛛丝蛋白特征较少")
        ),
    }


def read_review_jsonl(path: Path) -> dict[str, dict[str, Any]]:
    reviews: dict[str, dict[str, Any]] = {}
    if not path.exists():
        return reviews
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            try:
                payload = json.loads(line)
            except json.JSONDecodeError as exc:
                raise SystemExit(f"{path}:{line_number}: invalid JSON: {exc}") from exc
            record_id = payload.get("record_id")
            if not isinstance(record_id, str) or not record_id:
                raise SystemExit(f"{path}:{line_number}: missing string record_id")
            if record_id in reviews:
                raise SystemExit(f"{path}:{line_number}: duplicate record_id {record_id}")
            reviews[record_id] = normalize_review(payload)
    return reviews


def normalize_review(payload: dict[str, Any]) -> dict[str, Any]:
    normalized = dict(payload)
    normalized["matched_features"] = normalize_string_list(payload.get("matched_features"))
    normalized["concerns"] = normalize_string_list(payload.get("concerns"))
    normalized["confidence"] = normalize_confidence(payload.get("confidence"))
    normalized["is_spidroin_like"] = bool(payload.get("is_spidroin_like"))
    normalized.setdefault("type_consistency", "not_reported")
    normalized.setdefault("short_comment", "")
    normalized.setdefault("declared_type", "unknown")
    normalized.setdefault("species", "unknown")
    return normalized


def normalize_string_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return [str(item) for item in value if str(item).strip()]
    if isinstance(value, str):
        return [item.strip() for item in re.split(r"[;；|]", value) if item.strip()]
    return [str(value)]


def normalize_confidence(value: Any) -> float:
    if isinstance(value, (int, float)):
        numeric = float(value)
        return max(0.0, min(1.0, numeric / 100 if numeric > 1 else numeric))
    if isinstance(value, str):
        lower = value.strip().lower()
        if lower in {"high", "高", "strong"}:
            return 0.9
        if lower in {"medium", "中", "moderate"}:
            return 0.65
        if lower in {"low", "低", "weak"}:
            return 0.4
        try:
            numeric = float(lower.rstrip("%"))
        except ValueError:
            return 0.0
        return max(0.0, min(1.0, numeric / 100 if numeric > 1 else numeric))
    return 0.0


def validate_inputs(
    records: list[ProteinRecord],
    reviews: dict[str, dict[str, Any]],
    expected_count: int | None,
    require_reviews: bool,
) -> None:
    if expected_count is not None and len(records) != expected_count:
        raise SystemExit(f"Expected {expected_count} FASTA records, observed {len(records)}")
    record_ids = [record.record_id for record in records]
    duplicates = [rid for rid, count in Counter(record_ids).items() if count > 1]
    if duplicates:
        raise SystemExit(f"Duplicate FASTA record_id values: {duplicates[:5]}")
    known = set(record_ids)
    unknown_reviews = sorted(set(reviews) - known)
    if unknown_reviews:
        raise SystemExit(f"Review JSONL contains record_id not in FASTA: {unknown_reviews[:5]}")
    if require_reviews:
        missing = sorted(known - set(reviews))
        if missing:
            raise SystemExit(
                f"Missing {len(missing)} subagent reviews. First missing record_id: {missing[0]}"
            )


def fmt_int(value: int | None) -> str:
    return "NA" if value is None else f"{value:,}"


def fmt_pct(value: float | None, digits: int = 1) -> str:
    if value is None:
        return "NA"
    return f"{value * 100:.{digits}f}%"


def fmt_float(value: float | None, digits: int = 2) -> str:
    if value is None:
        return "NA"
    return f"{value:.{digits}f}"


def esc(value: Any) -> str:
    return html.escape(str(value), quote=True)


def average(values: Iterable[float]) -> float:
    items = list(values)
    return mean(items) if items else 0.0


def build_rows(records: list[ProteinRecord], reviews: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for record in records:
        features = calculate_features(record)
        review = reviews.get(record.record_id) or deterministic_agent_seed(record, features)
        rows.append(
            {
                "record": record,
                "features": features,
                "review": review,
                "needs_review": needs_review(record, features, review),
            }
        )
    return rows


def needs_review(
    record: ProteinRecord,
    features: dict[str, Any],
    review: dict[str, Any],
) -> bool:
    consistency = str(review.get("type_consistency", "")).lower()
    if features["internal_stop_count"] > 0:
        return True
    if not review.get("is_spidroin_like", False):
        return True
    if review.get("confidence", 0.0) < 0.62:
        return True
    if "inconsistent" in consistency or "conflict" in consistency or "low_" in consistency:
        return True
    if record.query_coverage is not None and record.query_coverage < 0.85:
        return True
    return False


def summarize_species(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[row["record"].species].append(row)
    summary = []
    for species, items in grouped.items():
        type_counts = Counter(row["record"].declared_type for row in items)
        confidence_values = [row["review"]["confidence"] for row in items]
        scores = [row["features"]["deterministic_score"] for row in items]
        summary.append(
            {
                "species": species,
                "count": len(items),
                "like": sum(1 for row in items if row["review"].get("is_spidroin_like")),
                "needs_review": sum(1 for row in items if row["needs_review"]),
                "avg_confidence": average(confidence_values),
                "avg_score": average(scores),
                "main_types": ", ".join(
                    f"{name}({count})" for name, count in type_counts.most_common(3)
                ),
            }
        )
    return sorted(summary, key=lambda item: (-item["count"], item["species"]))


def summarize_types(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[row["record"].declared_type].append(row)
    summary = []
    for declared_type, items in grouped.items():
        lengths = [row["features"]["length"] for row in items]
        motif_coverages = [row["features"]["motif_coverage"] for row in items]
        low_complexities = [row["features"]["low_complexity_fraction"] for row in items]
        motif_totals = Counter()
        for row in items:
            motif_totals.update(row["features"]["motif_counts"])
        summary.append(
            {
                "type": declared_type,
                "count": len(items),
                "like": sum(1 for row in items if row["review"].get("is_spidroin_like")),
                "needs_review": sum(1 for row in items if row["needs_review"]),
                "median_length": median(lengths),
                "avg_length": average(lengths),
                "avg_motif_coverage": average(motif_coverages),
                "avg_low_complexity": average(low_complexities),
                "top_motifs": ", ".join(
                    f"{name}({count})" for name, count in motif_totals.most_common(4) if count
                )
                or "NA",
            }
        )
    return sorted(summary, key=lambda item: (-item["count"], item["type"]))


def svg_arc_path(cx: float, cy: float, radius: float, start: float, end: float) -> str:
    start_rad = math.radians(start - 90)
    end_rad = math.radians(end - 90)
    x1 = cx + radius * math.cos(start_rad)
    y1 = cy + radius * math.sin(start_rad)
    x2 = cx + radius * math.cos(end_rad)
    y2 = cy + radius * math.sin(end_rad)
    large_arc = 1 if end - start > 180 else 0
    return f"M {x1:.2f} {y1:.2f} A {radius} {radius} 0 {large_arc} 1 {x2:.2f} {y2:.2f}"


def render_type_donut(type_summary: list[dict[str, Any]], total: int) -> str:
    cx, cy, radius = 150, 148, 78
    angle = 0.0
    parts = [
        '<svg class="chart" viewBox="0 0 620 320" role="img" '
        'aria-label="蛛丝蛋白类型构成 donut 图">',
        '<text x="24" y="30" class="svg-title">蛛丝蛋白类型构成</text>',
    ]
    for index, item in enumerate(type_summary):
        value = item["count"]
        delta = 360 * value / total
        color = TYPE_COLORS[index % len(TYPE_COLORS)]
        path = svg_arc_path(cx, cy, radius, angle, angle + delta)
        parts.append(
            f'<path d="{path}" fill="none" stroke="{color}" stroke-width="34" '
            'stroke-linecap="butt"/>'
        )
        angle += delta
    parts.append(f'<text x="{cx}" y="{cy - 5}" class="donut-total" text-anchor="middle">{total}</text>')
    parts.append(f'<text x="{cx}" y="{cy + 20}" class="svg-note" text-anchor="middle">sequences</text>')
    y = 68
    for index, item in enumerate(type_summary[:12]):
        color = TYPE_COLORS[index % len(TYPE_COLORS)]
        pct = item["count"] / total
        parts.append(f'<rect x="310" y="{y - 10}" width="14" height="14" rx="3" fill="{color}"/>')
        parts.append(
            f'<text x="332" y="{y + 1}" class="svg-label">{esc(item["type"])}</text>'
            f'<text x="500" y="{y + 1}" class="svg-value">{item["count"]} · {pct:.1%}</text>'
        )
        y += 22
    parts.append("</svg>")
    return "\n".join(parts)


def render_species_bars(species_summary: list[dict[str, Any]]) -> str:
    top = species_summary[:20]
    width = 760
    row_h = 24
    height = 74 + row_h * len(top)
    max_count = max((item["count"] for item in top), default=1)
    parts = [
        f'<svg class="chart" viewBox="0 0 {width} {height}" role="img" '
        'aria-label="物种序列数 top 20 条形图">',
        '<text x="24" y="30" class="svg-title">物种序列数 Top 20</text>',
        '<text x="24" y="52" class="svg-note">颜色深段为需要人工复核的条目</text>',
    ]
    for i, item in enumerate(top):
        y = 74 + i * row_h
        label = item["species"].split(".", 1)[-1]
        bar_w = 430 * item["count"] / max_count
        risk_w = bar_w * (item["needs_review"] / item["count"] if item["count"] else 0)
        parts.append(f'<text x="24" y="{y + 13}" class="svg-label">{esc(label[:32])}</text>')
        parts.append(f'<rect x="270" y="{y}" width="{bar_w:.1f}" height="15" rx="4" fill="#91c7bb"/>')
        if risk_w:
            parts.append(f'<rect x="270" y="{y}" width="{risk_w:.1f}" height="15" rx="4" fill="#d96c4f"/>')
        parts.append(f'<text x="{280 + bar_w:.1f}" y="{y + 12}" class="svg-value">{item["count"]}</text>')
    parts.append("</svg>")
    return "\n".join(parts)


def render_score_histogram(rows: list[dict[str, Any]]) -> str:
    scores = [row["features"]["deterministic_score"] for row in rows]
    bins = [0] * 10
    for score in scores:
        index = min(9, int(score // 10))
        bins[index] += 1
    width = 620
    height = 300
    chart_x, chart_y = 56, 58
    chart_w, chart_h = 520, 190
    max_count = max(bins) or 1
    bar_w = chart_w / len(bins) - 8
    parts = [
        f'<svg class="chart" viewBox="0 0 {width} {height}" role="img" '
        'aria-label="综合特征评分分布图">',
        '<text x="24" y="30" class="svg-title">综合特征评分分布</text>',
        f'<line x1="{chart_x}" y1="{chart_y + chart_h}" '
        f'x2="{chart_x + chart_w}" y2="{chart_y + chart_h}" class="axis"/>',
    ]
    for i, count in enumerate(bins):
        bar_h = chart_h * count / max_count
        x = chart_x + i * (chart_w / len(bins)) + 4
        y = chart_y + chart_h - bar_h
        fill = "#d96c4f" if i < 5 else "#2f9c95"
        parts.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{bar_h:.1f}" rx="4" fill="{fill}"/>')
        parts.append(f'<text x="{x + bar_w / 2:.1f}" y="{chart_y + chart_h + 18}" class="svg-note" text-anchor="middle">{i * 10}</text>')
        if count:
            parts.append(f'<text x="{x + bar_w / 2:.1f}" y="{y - 5:.1f}" class="svg-value" text-anchor="middle">{count}</text>')
    parts.append('<text x="560" y="276" class="svg-note" text-anchor="end">score</text>')
    parts.append("</svg>")
    return "\n".join(parts)


def render_motif_legend() -> str:
    parts = [
        '<svg class="chart" viewBox="0 0 620 250" role="img" aria-label="典型 motif 图例">',
        '<text x="24" y="30" class="svg-title">典型重复 motif 图示</text>',
        '<text x="24" y="52" class="svg-note">每个彩色块表示常见丝蛋白重复信号；灰色为 spacer 或未匹配区域</text>',
        '<rect x="36" y="92" width="548" height="28" rx="14" fill="#e8edf0"/>',
    ]
    x = 52
    segments = [
        ("(GA)n", 70),
        ("polyA", 86),
        ("GGX", 52),
        ("(GS)n", 66),
        ("GPGXX", 70),
        ("polyQ", 58),
    ]
    for name, size in segments:
        parts.append(
            f'<rect x="{x}" y="96" width="{size}" height="20" rx="10" '
            f'fill="{MOTIF_COLORS[name]}"/>'
        )
        x += size + 18
    y = 160
    x = 36
    for index, (name, color) in enumerate(MOTIF_COLORS.items()):
        if index and index % 4 == 0:
            y += 26
            x = 36
        parts.append(f'<rect x="{x}" y="{y - 12}" width="16" height="16" rx="4" fill="{color}"/>')
        parts.append(f'<text x="{x + 24}" y="{y + 1}" class="svg-label">{esc(name)}</text>')
        x += 140
    parts.append("</svg>")
    return "\n".join(parts)


def render_metric_cards(rows: list[dict[str, Any]], species_count: int, type_count: int) -> str:
    total = len(rows)
    like = sum(1 for row in rows if row["review"].get("is_spidroin_like"))
    needs = sum(1 for row in rows if row["needs_review"])
    avg_conf = average(row["review"]["confidence"] for row in rows)
    avg_score = average(row["features"]["deterministic_score"] for row in rows)
    metrics = [
        ("序列总数", f"{total:,}", "subagent 逐条审阅"),
        ("物种数", f"{species_count:,}", "含选中序列的物种"),
        ("类型数", f"{type_count:,}", "来自 typing_spidroin_type"),
        ("Spidroin-like", f"{like / total:.1%}", f"{like:,} / {total:,}"),
        ("需复核", f"{needs:,}", f"{needs / total:.1%} of sequences"),
        ("平均置信度", f"{avg_conf:.2f}", f"特征评分均值 {avg_score:.1f}"),
    ]
    return "\n".join(
        '<div class="metric"><div class="value">{}</div><div class="label">{}</div>'
        '<div class="hint">{}</div></div>'.format(esc(value), esc(label), esc(hint))
        for label, value, hint in metrics
    )


def render_species_table(species_summary: list[dict[str, Any]]) -> str:
    rows = []
    for item in species_summary:
        rows.append(
            "<tr>"
            f"<td>{esc(item['species'])}</td>"
            f"<td>{item['count']}</td>"
            f"<td>{item['like']}</td>"
            f"<td>{item['needs_review']}</td>"
            f"<td>{item['avg_confidence']:.2f}</td>"
            f"<td>{item['avg_score']:.1f}</td>"
            f"<td>{esc(item['main_types'])}</td>"
            "</tr>"
        )
    return "\n".join(rows)


def render_type_table(type_summary: list[dict[str, Any]]) -> str:
    rows = []
    for item in type_summary:
        rows.append(
            "<tr>"
            f"<td>{esc(item['type'])}</td>"
            f"<td>{item['count']}</td>"
            f"<td>{item['like']}</td>"
            f"<td>{item['needs_review']}</td>"
            f"<td>{item['median_length']:.0f}</td>"
            f"<td>{item['avg_length']:.0f}</td>"
            f"<td>{fmt_pct(item['avg_motif_coverage'])}</td>"
            f"<td>{fmt_pct(item['avg_low_complexity'])}</td>"
            f"<td>{esc(item['top_motifs'])}</td>"
            "</tr>"
        )
    return "\n".join(rows)


def render_sequence_details(rows: list[dict[str, Any]]) -> str:
    rendered = []
    for row in rows:
        record: ProteinRecord = row["record"]
        features = row["features"]
        review = row["review"]
        status = "ok" if review.get("is_spidroin_like") and not row["needs_review"] else "review"
        matched = review.get("matched_features", [])
        concerns = review.get("concerns", [])
        matched_html = "".join(f"<li>{esc(item)}</li>" for item in matched) or "<li>未报告</li>"
        concerns_html = "".join(f"<li>{esc(item)}</li>" for item in concerns) or "<li>无显著疑点</li>"
        motif_text = ", ".join(
            f"{name}:{count}" for name, count in features["motif_counts"].items() if count
        ) or "none"
        composition = features["composition"]
        search_blob = " ".join(
            [
                record.record_id,
                record.species,
                record.declared_type,
                record.mpid,
                str(review.get("short_comment", "")),
                " ".join(matched),
                " ".join(concerns),
            ]
        )
        rendered.append(
            f'<details class="seq-detail {status}" data-record-id="{esc(record.record_id)}" '
            f'data-search="{esc(search_blob.lower())}">'
            "<summary>"
            f'<span class="badge {status}">{ "通过" if status == "ok" else "复核" }</span>'
            f'<span class="seq-title">{esc(record.species)} · {esc(record.mpid)} · {esc(record.declared_type)}</span>'
            f'<span class="seq-meta">score {features["deterministic_score"]:.1f} · conf {review["confidence"]:.2f}</span>'
            "</summary>"
            '<div class="seq-grid">'
            '<div>'
            f'<p class="comment">{esc(review.get("short_comment", ""))}</p>'
            '<h4>Subagent matched features</h4>'
            f"<ul>{matched_html}</ul>"
            '<h4>Concerns</h4>'
            f"<ul>{concerns_html}</ul>"
            "</div>"
            '<div class="facts">'
            f"<p><b>record_id</b><br>{esc(record.record_id)}</p>"
            f"<p><b>length</b> {features['length']:,} aa; header aa_len {fmt_int(record.aa_len_header)}</p>"
            f"<p><b>positive</b> {fmt_float(record.positive, 4)}; <b>coverage</b> {fmt_float(record.query_coverage, 4)}</p>"
            f"<p><b>A/G/S</b> {fmt_pct(features['ags_fraction'])}; "
            f"A {fmt_pct(composition['A'])}, G {fmt_pct(composition['G'])}, "
            f"S {fmt_pct(composition['S'])}</p>"
            f"<p><b>motif coverage</b> {fmt_pct(features['motif_coverage'])}; "
            f"<b>low complexity</b> {fmt_pct(features['low_complexity_fraction'])}</p>"
            f"<p><b>motifs</b><br>{esc(motif_text)}</p>"
            f"<p><b>type consistency</b><br>{esc(review.get('type_consistency', 'not_reported'))}</p>"
            f"<p><b>stops</b> terminal_only={esc(record.terminal_stop_only)}; "
            f"internal={features['internal_stop_count']}</p>"
            "</div>"
            "</div>"
            "</details>"
        )
    return "\n".join(rendered)


def render_html(
    *,
    rows: list[dict[str, Any]],
    species_summary: list[dict[str, Any]],
    type_summary: list[dict[str, Any]],
    review_jsonl: Path,
    fasta_path: Path,
) -> str:
    generated = datetime.now().strftime("%Y-%m-%d %H:%M")
    total = len(rows)
    species_total = len(species_summary)
    type_total = len(type_summary)
    type_donut = render_type_donut(type_summary, total)
    species_bars = render_species_bars(species_summary)
    score_hist = render_score_histogram(rows)
    motif_legend = render_motif_legend()
    metric_cards = render_metric_cards(rows, species_total, type_total)
    species_rows = render_species_table(species_summary)
    type_rows = render_type_table(type_summary)
    detail_rows = render_sequence_details(rows)
    return f"""<!doctype html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Subagent Spidroin Review - miniprot_window_spidroin_screen_20260526_115407</title>
<style>
:root {{
  --ink: #20262d;
  --muted: #66717d;
  --line: #dce4e8;
  --bg: #f7f9f8;
  --panel: #ffffff;
  --teal: #2f9c95;
  --orange: #d96c4f;
  --green: #58a05c;
  --blue: #4d84c4;
  --gold: #c59a2e;
}}
* {{ box-sizing: border-box; }}
body {{
  margin: 0;
  color: var(--ink);
  background: var(--bg);
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", "Noto Sans CJK SC",
    "Microsoft YaHei", Arial, sans-serif;
  line-height: 1.55;
}}
header {{
  background: #24313a;
  color: white;
  padding: 34px 42px 30px;
}}
header h1 {{ margin: 0 0 8px; font-size: 30px; font-weight: 760; }}
header p {{ margin: 0; color: #d9e5e8; max-width: 980px; }}
main {{ max-width: 1280px; margin: 0 auto; padding: 28px 28px 60px; }}
.band {{
  background: var(--panel);
  border: 1px solid var(--line);
  border-radius: 8px;
  padding: 22px;
  margin-bottom: 20px;
}}
h2 {{ margin: 0 0 14px; font-size: 22px; }}
h3 {{ margin: 20px 0 10px; font-size: 16px; }}
h4 {{ margin: 14px 0 6px; font-size: 13px; }}
.muted {{ color: var(--muted); }}
.metrics {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(168px, 1fr)); gap: 12px; }}
.metric {{ border-left: 4px solid var(--teal); background: #fbfcfc; padding: 13px 14px; border-radius: 8px; }}
.metric .value {{ font-size: 28px; font-weight: 800; color: var(--ink); }}
.metric .label {{ color: var(--muted); font-size: 13px; }}
.metric .hint {{ color: var(--muted); font-size: 12px; margin-top: 3px; }}
.charts {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(460px, 1fr)); gap: 16px; align-items: start; }}
.chart {{ width: 100%; height: auto; display: block; background: #fbfcfc; border: 1px solid var(--line); border-radius: 8px; }}
.svg-title {{ font-size: 18px; font-weight: 760; fill: var(--ink); }}
.svg-label {{ font-size: 12px; fill: var(--ink); }}
.svg-value, .svg-note {{ font-size: 12px; fill: var(--muted); }}
.donut-total {{ font-size: 28px; font-weight: 800; fill: var(--ink); }}
.axis {{ stroke: #aab7c1; stroke-width: 1; }}
table {{ width: 100%; border-collapse: collapse; font-size: 13px; }}
th, td {{ border-bottom: 1px solid var(--line); padding: 8px 9px; text-align: left; vertical-align: top; }}
th {{ background: #edf3f3; font-weight: 730; position: sticky; top: 0; z-index: 1; }}
.table-wrap {{ max-height: 520px; overflow: auto; border: 1px solid var(--line); border-radius: 8px; }}
.controls {{ display: flex; gap: 12px; align-items: center; margin-bottom: 12px; flex-wrap: wrap; }}
input[type="search"] {{
  width: min(560px, 100%);
  padding: 10px 12px;
  border: 1px solid var(--line);
  border-radius: 8px;
  font-size: 14px;
}}
.seq-detail {{
  border: 1px solid var(--line);
  border-radius: 8px;
  background: #fff;
  margin: 9px 0;
  overflow: hidden;
}}
.seq-detail.review {{ border-left: 4px solid var(--orange); }}
.seq-detail.ok {{ border-left: 4px solid var(--teal); }}
.seq-detail[hidden] {{ display: none; }}
summary {{
  cursor: pointer;
  list-style: none;
  display: grid;
  grid-template-columns: auto minmax(0, 1fr) auto;
  gap: 10px;
  align-items: center;
  padding: 11px 13px;
}}
summary::-webkit-details-marker {{ display: none; }}
.badge {{
  display: inline-flex;
  align-items: center;
  justify-content: center;
  min-width: 42px;
  height: 24px;
  border-radius: 999px;
  color: white;
  font-size: 12px;
  font-weight: 730;
}}
.badge.ok {{ background: var(--teal); }}
.badge.review {{ background: var(--orange); }}
.seq-title {{ overflow: hidden; text-overflow: ellipsis; white-space: nowrap; font-weight: 680; }}
.seq-meta {{ color: var(--muted); font-size: 12px; white-space: nowrap; }}
.seq-grid {{ display: grid; grid-template-columns: minmax(0, 1.3fr) minmax(280px, 0.7fr); gap: 18px; padding: 0 16px 16px; }}
.comment {{ margin: 0 0 8px; color: var(--ink); }}
.facts {{ background: #f7faf9; border: 1px solid var(--line); border-radius: 8px; padding: 12px; font-size: 12px; }}
.facts p {{ margin: 0 0 8px; overflow-wrap: anywhere; }}
ul {{ margin-top: 4px; padding-left: 20px; }}
footer {{ color: var(--muted); text-align: center; padding: 20px 0 36px; font-size: 12px; }}
@media (max-width: 720px) {{
  header {{ padding: 26px 22px; }}
  main {{ padding: 20px 14px 44px; }}
  .charts {{ grid-template-columns: 1fr; }}
  .seq-grid {{ grid-template-columns: 1fr; }}
  summary {{ grid-template-columns: auto minmax(0, 1fr); }}
  .seq-meta {{ grid-column: 2; }}
}}
</style>
</head>
<body>
<header>
  <h1>Subagent Spidroin Review</h1>
  <p>1054 条候选蛛丝蛋白逐条由 subagent 审阅；本页合并 deterministic 特征、agent 判断、物种聚合和类型聚合。生成时间：{esc(generated)}</p>
</header>
<main>
  <section class="band">
    <h2>总览</h2>
    <div class="metrics">{metric_cards}</div>
    <p class="muted">FASTA: {esc(fasta_path)} · Subagent JSONL: {esc(review_jsonl)}</p>
  </section>

  <section class="band">
    <h2>图示摘要</h2>
    <div class="charts">
      {type_donut}
      {species_bars}
      {score_hist}
      {motif_legend}
    </div>
  </section>

  <section class="band">
    <h2>按物种总结</h2>
    <div class="table-wrap">
      <table>
        <thead><tr><th>物种</th><th>序列</th><th>spidroin-like</th><th>需复核</th><th>平均置信度</th><th>平均评分</th><th>主要类型</th></tr></thead>
        <tbody>{species_rows}</tbody>
      </table>
    </div>
  </section>

  <section class="band">
    <h2>按蛛丝蛋白种类总结</h2>
    <div class="table-wrap">
      <table>
        <thead><tr><th>类型</th><th>序列</th><th>spidroin-like</th><th>需复核</th><th>长度中位数</th><th>平均长度</th><th>motif 覆盖</th><th>低复杂度</th><th>主要 motif</th></tr></thead>
        <tbody>{type_rows}</tbody>
      </table>
    </div>
  </section>

  <section class="band">
    <h2>逐序列 subagent 明细</h2>
    <div class="controls">
      <input id="sequenceSearch" type="search" placeholder="搜索 species / MPID / 类型 / comment">
      <span class="muted"><span id="visibleCount">{total}</span> / {total} 条</span>
    </div>
    {detail_rows}
  </section>
</main>
<footer>spider_silkome · standalone HTML with embedded SVG</footer>
<script>
const search = document.getElementById('sequenceSearch');
const visibleCount = document.getElementById('visibleCount');
const details = Array.from(document.querySelectorAll('.seq-detail'));
search.addEventListener('input', () => {{
  const q = search.value.trim().toLowerCase();
  let shown = 0;
  for (const item of details) {{
    const ok = !q || item.dataset.search.includes(q);
    item.hidden = !ok;
    if (ok) shown += 1;
  }}
  visibleCount.textContent = shown;
}});
</script>
</body>
</html>
"""


def write_agent_input_jsonl(path: Path, records: list[ProteinRecord]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for record in records:
            features = calculate_features(record)
            seed = deterministic_agent_seed(record, features)
            payload = {
                "record_id": record.record_id,
                "header": record.header,
                "sequence": record.sequence,
                "features": {
                    "length": features["length"],
                    "ags_fraction": features["ags_fraction"],
                    "motif_counts": features["motif_counts"],
                    "motif_coverage": features["motif_coverage"],
                    "low_complexity_fraction": features["low_complexity_fraction"],
                    "internal_stop_count": features["internal_stop_count"],
                    "deterministic_score": features["deterministic_score"],
                    "feature_type_hint": features["feature_type_hint"],
                },
                "suggested_review_shape": seed,
            }
            handle.write(json.dumps(payload, ensure_ascii=False, sort_keys=True) + "\n")


def write_seed_review_jsonl(path: Path, records: list[ProteinRecord]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for record in records:
            features = calculate_features(record)
            payload = deterministic_agent_seed(record, features)
            payload["short_comment"] = "DETERMINISTIC_SEED_ONLY: " + payload["short_comment"]
            handle.write(json.dumps(payload, ensure_ascii=False, sort_keys=True) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, default=DEFAULT_FASTA)
    parser.add_argument("--review-jsonl", type=Path, default=DEFAULT_REVIEW_JSONL)
    parser.add_argument("--html", type=Path, default=DEFAULT_HTML)
    parser.add_argument("--expected-count", type=int, default=1054)
    parser.add_argument(
        "--allow-partial-reviews",
        action="store_true",
        help="Render with deterministic fallback for records missing from review JSONL.",
    )
    parser.add_argument(
        "--agent-input-jsonl",
        type=Path,
        help="Optional path to write per-record subagent inputs, including full sequences.",
    )
    parser.add_argument(
        "--seed-review-jsonl",
        type=Path,
        help="Optional path to write deterministic seed reviews. Use only for debugging.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    records = read_fasta(args.fasta)
    if args.agent_input_jsonl:
        write_agent_input_jsonl(args.agent_input_jsonl, records)
    if args.seed_review_jsonl:
        write_seed_review_jsonl(args.seed_review_jsonl, records)
    reviews = read_review_jsonl(args.review_jsonl)
    validate_inputs(
        records=records,
        reviews=reviews,
        expected_count=args.expected_count,
        require_reviews=not args.allow_partial_reviews,
    )
    rows = build_rows(records, reviews)
    species_summary = summarize_species(rows)
    type_summary = summarize_types(rows)
    if sum(item["count"] for item in species_summary) != len(records):
        raise SystemExit("Species aggregation does not sum to FASTA record count")
    if sum(item["count"] for item in type_summary) != len(records):
        raise SystemExit("Type aggregation does not sum to FASTA record count")
    args.html.parent.mkdir(parents=True, exist_ok=True)
    args.html.write_text(
        render_html(
            rows=rows,
            species_summary=species_summary,
            type_summary=type_summary,
            review_jsonl=args.review_jsonl,
            fasta_path=args.fasta,
        ),
        encoding="utf-8",
    )
    print(f"Wrote {args.html} with {len(records)} sequence details")


if __name__ == "__main__":
    main()
