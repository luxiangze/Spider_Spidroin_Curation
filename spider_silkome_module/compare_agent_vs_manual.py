"""Compare automated spidroin typing agent results with manual curation table.

Pipeline:
  1. Load manual CSV and agent per-species TSVs; normalize species / chr / type fields
  2. For each (species, chr, strand) group, greedy-match loci by reciprocal interval
     overlap (>= min_overlap on both sides)
  3. Emit a long-form merged table (matched / manual_only / agent_only) plus several
     summary tables (species coverage, type confusion, boundary precision, Scoring-
     stratified recall, agent-confidence x Scoring crosstab) and a markdown report
"""

from __future__ import annotations

from pathlib import Path

from loguru import logger
import polars as pl
import typer

from spider_silkome_module.config import EXTERNAL_DATA_DIR, PROCESSED_DATA_DIR, PROJ_ROOT

app = typer.Typer()


# ---------------------------------------------------------------------------
# Normalization helpers
# ---------------------------------------------------------------------------

# Manual table uses both 'Unkown' (typo) and 'hypo'; agent uses 'unknown'.
_TYPE_ALIASES = {
    "unkown": "unknown",
    "hypo": "unknown",
}


def _norm_species(s: str | None) -> str | None:
    """Manual uses '_' separator, agent uses ' '; unify to space."""
    return None if s is None else s.replace("_", " ").strip()


def _norm_chr(s: str | None) -> str | None:
    """Strip whitespace + lowercase so 'Chr03', '\\nchr03' all collide."""
    return None if s is None else s.strip().lower()


def _norm_type(s: str | None) -> str | None:
    """Lowercase + dealias; keep compound types ('CySp/MiSp') as-is."""
    if s is None:
        return None
    key = s.strip().lower()
    return _TYPE_ALIASES.get(key, key)


def _norm_hint(s: str | None) -> str | None:
    return None if s is None else s.strip()


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------


def load_manual(csv_path: Path) -> pl.DataFrame:
    """Load + normalize the manual curation CSV."""
    df = pl.read_csv(csv_path)
    df = df.with_columns(
        pl.col("Species").map_elements(_norm_species, return_dtype=pl.String).alias("species"),
        pl.col("Chr").map_elements(_norm_chr, return_dtype=pl.String).alias("chr"),
        pl.col("Stran").alias("strand"),
        pl.col("Start").cast(pl.Int64, strict=False).alias("start"),
        pl.col("End").cast(pl.Int64, strict=False).alias("end"),
        pl.col("Spidroin_type").map_elements(_norm_type, return_dtype=pl.String).alias("type"),
        pl.col("Full_length").alias("full_length"),
        pl.col("Hint_type").map_elements(_norm_hint, return_dtype=pl.String).alias("hint_type"),
        pl.col("Scoring").alias("scoring"),
        pl.col("Spidroin_ID").alias("locus_id"),
    )
    return df.filter(
        pl.col("species").is_not_null()
        & pl.col("chr").is_not_null()
        & pl.col("start").is_not_null()
        & pl.col("end").is_not_null()
    )


def load_agent(typing_dir: Path) -> pl.DataFrame:
    """Concat all <species>/<species>.tsv files under typing_dir."""
    frames: list[pl.DataFrame] = []
    for tsv in sorted(typing_dir.glob("*/[0-9]*.tsv")):
        # Skip aux files: only main per-species TSV named like <species>.tsv inside a
        # folder of the same name (parent.name matches stem)
        if tsv.stem != tsv.parent.name:
            continue
        frames.append(pl.read_csv(tsv, separator="\t", infer_schema_length=10000))
    if not frames:
        raise FileNotFoundError(f"No agent TSVs found under {typing_dir}")
    df = pl.concat(frames, how="diagonal_relaxed")
    df = df.with_columns(
        pl.col("Species").map_elements(_norm_species, return_dtype=pl.String).alias("species"),
        pl.col("Chr").map_elements(_norm_chr, return_dtype=pl.String).alias("chr"),
        pl.col("Stran").alias("strand"),
        pl.col("Start").cast(pl.Int64).alias("start"),
        pl.col("End").cast(pl.Int64).alias("end"),
        pl.col("Spidroin_type").map_elements(_norm_type, return_dtype=pl.String).alias("type"),
        pl.col("Full_length").alias("full_length"),
        pl.col("Hint_type").map_elements(_norm_hint, return_dtype=pl.String).alias("hint_type"),
        pl.col("Spidroin_ID").alias("locus_id"),
    )
    return df


# ---------------------------------------------------------------------------
# Interval matching
# ---------------------------------------------------------------------------


def _reciprocal_overlap(
    a_start: int, a_end: int, b_start: int, b_end: int
) -> tuple[float, float]:
    """Return (a_in_b, b_in_a) = (overlap / len(a), overlap / len(b))."""
    ov = max(0, min(a_end, b_end) - max(a_start, b_start))
    la = max(1, a_end - a_start)
    lb = max(1, b_end - b_start)
    return ov / la, ov / lb


def _greedy_match_group(
    manual_rows: list[tuple], agent_rows: list[tuple], min_overlap: float
) -> tuple[list[tuple], set, set]:
    """Greedy max-weight 1-to-1 matching within one (species, chr, strand) group.

    Each row is (idx, start, end). Returns (matched_pairs, unmatched_m, unmatched_a)
    where matched_pairs = [(m_idx, a_idx, m_in_a, a_in_m), ...].
    """
    candidates: list[tuple[float, int, int, float, float]] = []
    for mi, ms, me in manual_rows:
        for ai, asx, ae in agent_rows:
            m_in_a, a_in_m = _reciprocal_overlap(ms, me, asx, ae)
            if min(m_in_a, a_in_m) >= min_overlap:
                candidates.append((min(m_in_a, a_in_m), mi, ai, m_in_a, a_in_m))
    candidates.sort(reverse=True)

    used_m: set[int] = set()
    used_a: set[int] = set()
    matched: list[tuple] = []
    for _score, mi, ai, m_in_a, a_in_m in candidates:
        if mi in used_m or ai in used_a:
            continue
        used_m.add(mi)
        used_a.add(ai)
        matched.append((mi, ai, m_in_a, a_in_m))
    unmatched_m = {r[0] for r in manual_rows} - used_m
    unmatched_a = {r[0] for r in agent_rows} - used_a
    return matched, unmatched_m, unmatched_a


def match_loci(
    manual: pl.DataFrame, agent: pl.DataFrame, min_overlap: float = 0.5
) -> pl.DataFrame:
    """Long-form merged table: one row per matched pair or single-side locus."""
    manual = manual.with_row_index("m_idx")
    agent = agent.with_row_index("a_idx")

    species_set = set(manual["species"].unique()) | set(agent["species"].unique())
    rows: list[dict] = []

    for sp in sorted(s for s in species_set if s is not None):
        m_sp = manual.filter(pl.col("species") == sp)
        a_sp = agent.filter(pl.col("species") == sp)

        keys_m = set(zip(m_sp["chr"].to_list(), m_sp["strand"].to_list()))
        keys_a = set(zip(a_sp["chr"].to_list(), a_sp["strand"].to_list()))
        for chr_, strand in sorted(keys_m | keys_a, key=lambda x: (str(x[0]), str(x[1]))):
            m_grp = m_sp.filter((pl.col("chr") == chr_) & (pl.col("strand") == strand))
            a_grp = a_sp.filter((pl.col("chr") == chr_) & (pl.col("strand") == strand))

            m_rows = list(zip(m_grp["m_idx"].to_list(), m_grp["start"].to_list(), m_grp["end"].to_list()))
            a_rows = list(zip(a_grp["a_idx"].to_list(), a_grp["start"].to_list(), a_grp["end"].to_list()))

            matched, un_m, un_a = _greedy_match_group(m_rows, a_rows, min_overlap)
            for mi, ai, m_in_a, a_in_m in matched:
                rows.append(_merged_row(manual, agent, mi, ai, m_in_a, a_in_m))
            for mi in un_m:
                rows.append(_merged_row(manual, agent, mi, None, None, None))
            for ai in un_a:
                rows.append(_merged_row(manual, agent, None, ai, None, None))

    if not rows:
        return pl.DataFrame()
    return pl.DataFrame(rows)


def _merged_row(
    manual: pl.DataFrame,
    agent: pl.DataFrame,
    m_idx: int | None,
    a_idx: int | None,
    m_in_a: float | None,
    a_in_m: float | None,
) -> dict:
    """Build one merged-row dict given optional manual + agent row indices."""
    m = manual.row(by_predicate=pl.col("m_idx") == m_idx, named=True) if m_idx is not None else None
    a = agent.row(by_predicate=pl.col("a_idx") == a_idx, named=True) if a_idx is not None else None

    status = (
        "matched" if m and a else ("manual_only" if m else "agent_only")
    )
    species = (m or a)["species"]
    chr_ = (m or a)["chr"]
    strand = (m or a)["strand"]

    def _get(d, k):
        return d[k] if d is not None else None

    start_diff = None
    end_diff = None
    if m and a:
        start_diff = a["start"] - m["start"]
        end_diff = a["end"] - m["end"]

    return {
        "species": species,
        "chr": chr_,
        "strand": strand,
        "match_status": status,
        "manual_id": _get(m, "locus_id"),
        "manual_start": _get(m, "start"),
        "manual_end": _get(m, "end"),
        "manual_type": _get(m, "type"),
        "manual_full_length": _get(m, "full_length"),
        "manual_hint_type": _get(m, "hint_type"),
        "manual_scoring": _get(m, "scoring"),
        "agent_id": _get(a, "locus_id"),
        "agent_start": _get(a, "start"),
        "agent_end": _get(a, "end"),
        "agent_type": _get(a, "type"),
        "agent_full_length": _get(a, "full_length"),
        "agent_hint_type": _get(a, "hint_type"),
        "agent_confidence": _get(a, "confidence"),
        "agent_needs_review": _get(a, "needs_review"),
        "agent_miniprot_identity": _get(a, "miniprot_identity"),
        "manual_in_agent_frac": m_in_a,
        "agent_in_manual_frac": a_in_m,
        "start_diff": start_diff,
        "end_diff": end_diff,
        "type_agree": (m and a and m["type"] == a["type"]) if m and a else None,
        "full_length_agree": (m and a and m["full_length"] == a["full_length"]) if m and a else None,
        "hint_type_agree": (m and a and m["hint_type"] == a["hint_type"]) if m and a else None,
    }


# ---------------------------------------------------------------------------
# Summary tables
# ---------------------------------------------------------------------------


def species_coverage(manual: pl.DataFrame, agent: pl.DataFrame) -> pl.DataFrame:
    """Per-species: manual_count, agent_count, side (both / manual_only / agent_only)."""
    m = manual.group_by("species").len().rename({"len": "manual_count"})
    a = agent.group_by("species").len().rename({"len": "agent_count"})
    joined = m.join(a, on="species", how="full", coalesce=True)
    joined = joined.with_columns(
        pl.col("manual_count").fill_null(0),
        pl.col("agent_count").fill_null(0),
    ).with_columns(
        pl.when((pl.col("manual_count") > 0) & (pl.col("agent_count") > 0))
        .then(pl.lit("both"))
        .when(pl.col("manual_count") > 0)
        .then(pl.lit("manual_only"))
        .otherwise(pl.lit("agent_only"))
        .alias("side")
    )
    return joined.sort("species")


def species_summary(merged: pl.DataFrame) -> pl.DataFrame:
    """Per-species TP / FP / FN / precision / recall / F1."""
    grp = merged.group_by("species").agg(
        tp=(pl.col("match_status") == "matched").sum(),
        fn=(pl.col("match_status") == "manual_only").sum(),
        fp=(pl.col("match_status") == "agent_only").sum(),
    )
    return grp.with_columns(
        precision=pl.when(pl.col("tp") + pl.col("fp") > 0)
        .then(pl.col("tp") / (pl.col("tp") + pl.col("fp")))
        .otherwise(None),
        recall=pl.when(pl.col("tp") + pl.col("fn") > 0)
        .then(pl.col("tp") / (pl.col("tp") + pl.col("fn")))
        .otherwise(None),
    ).with_columns(
        f1=pl.when((pl.col("precision").is_not_null()) & (pl.col("recall").is_not_null()) & (pl.col("precision") + pl.col("recall") > 0))
        .then(2 * pl.col("precision") * pl.col("recall") / (pl.col("precision") + pl.col("recall")))
        .otherwise(None)
    ).sort("recall")


def _confusion(merged: pl.DataFrame, manual_col: str, agent_col: str) -> pl.DataFrame:
    """Build a long-form confusion table (manual_value, agent_value, count) on matched rows."""
    matched = merged.filter(pl.col("match_status") == "matched")
    if matched.is_empty():
        return pl.DataFrame({"manual": [], "agent": [], "count": []})
    return (
        matched.group_by([manual_col, agent_col])
        .len()
        .rename({manual_col: "manual", agent_col: "agent", "len": "count"})
        .sort(["manual", "count"], descending=[False, True])
    )


def boundary_stats(merged: pl.DataFrame) -> pl.DataFrame:
    """Distribution of start_diff / end_diff (agent - manual) on matched pairs."""
    matched = merged.filter(pl.col("match_status") == "matched")
    if matched.is_empty():
        return pl.DataFrame()
    desc = matched.select(["start_diff", "end_diff"]).describe()
    return desc


def scoring_recall(merged: pl.DataFrame) -> pl.DataFrame:
    """For each manual Scoring level, n_total / n_matched / recall."""
    side = merged.filter(pl.col("match_status").is_in(["matched", "manual_only"]))
    return (
        side.group_by("manual_scoring")
        .agg(
            total=pl.len(),
            matched=(pl.col("match_status") == "matched").sum(),
        )
        .with_columns(recall=pl.col("matched") / pl.col("total"))
        .sort("manual_scoring")
    )


def confidence_correlation(merged: pl.DataFrame) -> pl.DataFrame:
    """Agent confidence x manual Scoring crosstab on matched pairs."""
    matched = merged.filter(pl.col("match_status") == "matched")
    if matched.is_empty():
        return pl.DataFrame()
    return (
        matched.group_by(["agent_confidence", "manual_scoring"])
        .len()
        .rename({"len": "count"})
        .sort(["agent_confidence", "manual_scoring"])
    )


# ---------------------------------------------------------------------------
# Report rendering
# ---------------------------------------------------------------------------


def render_report(
    merged: pl.DataFrame,
    coverage: pl.DataFrame,
    sp_summary: pl.DataFrame,
    type_conf: pl.DataFrame,
    fl_conf: pl.DataFrame,
    hint_conf: pl.DataFrame,
    boundary: pl.DataFrame,
    score_recall: pl.DataFrame,
    conf_corr: pl.DataFrame,
    min_overlap: float,
) -> str:
    """Produce a human-readable Markdown summary."""
    n_matched = (merged["match_status"] == "matched").sum()
    n_manual_only = (merged["match_status"] == "manual_only").sum()
    n_agent_only = (merged["match_status"] == "agent_only").sum()
    precision = n_matched / max(1, n_matched + n_agent_only)
    recall = n_matched / max(1, n_matched + n_manual_only)
    f1 = 2 * precision * recall / max(1e-9, precision + recall)

    type_agree = merged.filter(pl.col("type_agree").is_not_null())["type_agree"].sum()
    fl_agree = merged.filter(pl.col("full_length_agree").is_not_null())["full_length_agree"].sum()
    hint_agree = merged.filter(pl.col("hint_type_agree").is_not_null())["hint_type_agree"].sum()
    n_pair = max(1, n_matched)

    coverage_count = coverage.group_by("side").len()

    lines: list[str] = [
        "# Agent vs Manual Spidroin Typing Comparison",
        "",
        f"Reciprocal overlap threshold: **{min_overlap:.0%}**",
        "",
        "## Locus-level metrics (overall)",
        "",
        f"- True positives (matched): **{n_matched}**",
        f"- False negatives (manual_only): **{n_manual_only}**",
        f"- False positives (agent_only): **{n_agent_only}**",
        f"- Precision: **{precision:.3f}**",
        f"- Recall: **{recall:.3f}**",
        f"- F1: **{f1:.3f}**",
        "",
        "## Classification agreement on matched pairs",
        "",
        f"- Spidroin_type identical: {type_agree} / {n_pair} ({type_agree / n_pair:.1%})",
        f"- Full_length identical: {fl_agree} / {n_pair} ({fl_agree / n_pair:.1%})",
        f"- Hint_type identical: {hint_agree} / {n_pair} ({hint_agree / n_pair:.1%})",
        "",
        "## Species coverage",
        "",
        coverage_count.to_pandas().to_markdown(index=False),
        "",
        "## Recall stratified by manual Scoring",
        "",
        score_recall.to_pandas().to_markdown(index=False, floatfmt=".3f"),
        "",
        "## Boundary offsets (matched pairs, agent - manual)",
        "",
        boundary.to_pandas().to_markdown(index=False),
        "",
        "## Agent confidence vs manual Scoring crosstab",
        "",
        conf_corr.to_pandas().to_markdown(index=False) if not conf_corr.is_empty() else "_no matched pairs_",
        "",
        "## Bottom 10 species by Recall",
        "",
        sp_summary.head(10).to_pandas().to_markdown(index=False, floatfmt=".3f"),
        "",
        "## Top Spidroin_type confusions",
        "",
        type_conf.head(20).to_pandas().to_markdown(index=False),
        "",
        "## Full_length confusion matrix",
        "",
        fl_conf.to_pandas().to_markdown(index=False),
        "",
        "## Hint_type confusion matrix",
        "",
        hint_conf.to_pandas().to_markdown(index=False),
        "",
    ]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# CLI entry
# ---------------------------------------------------------------------------


@app.command()
def main(
    manual_csv: Path = EXTERNAL_DATA_DIR / "蛛丝蛋白鉴定_数据表_总表.csv",
    agent_dir: Path = PROCESSED_DATA_DIR / "typing_results",
    output_dir: Path = PROJ_ROOT / "results" / "comparison",
    min_overlap: float = 0.5,
):
    """Compare agent typing TSVs with the manual curation CSV."""
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Loading manual table: {manual_csv}")
    manual = load_manual(manual_csv)
    logger.info(f"  manual: {len(manual)} loci across {manual['species'].n_unique()} species")

    logger.info(f"Loading agent typing results: {agent_dir}")
    agent = load_agent(agent_dir)
    logger.info(f"  agent : {len(agent)} loci across {agent['species'].n_unique()} species")

    logger.info("Matching loci (reciprocal overlap >= %.2f)" % min_overlap)
    merged = match_loci(manual, agent, min_overlap=min_overlap)
    logger.info(
        "  matched={} | manual_only={} | agent_only={}".format(
            (merged["match_status"] == "matched").sum(),
            (merged["match_status"] == "manual_only").sum(),
            (merged["match_status"] == "agent_only").sum(),
        )
    )

    coverage = species_coverage(manual, agent)
    sp_summary = species_summary(merged)
    type_conf = _confusion(merged, "manual_type", "agent_type")
    fl_conf = _confusion(merged, "manual_full_length", "agent_full_length")
    hint_conf = _confusion(merged, "manual_hint_type", "agent_hint_type")
    boundary = boundary_stats(merged)
    score_recall = scoring_recall(merged)
    conf_corr = confidence_correlation(merged)

    discrepancies = merged.filter(
        (pl.col("match_status") != "matched")
        | (pl.col("type_agree") == False)  # noqa: E712
        | (pl.col("full_length_agree") == False)  # noqa: E712
        | (pl.col("hint_type_agree") == False)  # noqa: E712
    )

    for name, df in {
        "merged_loci.tsv": merged,
        "species_coverage.tsv": coverage,
        "species_summary.tsv": sp_summary,
        "type_confusion.tsv": type_conf,
        "full_length_confusion.tsv": fl_conf,
        "hint_type_confusion.tsv": hint_conf,
        "boundary_stats.tsv": boundary,
        "scoring_recall.tsv": score_recall,
        "confidence_correlation.tsv": conf_corr,
        "discrepancies.tsv": discrepancies,
    }.items():
        out = output_dir / name
        if df.is_empty():
            out.write_text("")
        else:
            df.write_csv(out, separator="\t")
        logger.info(f"  wrote {out} ({len(df)} rows)")

    report = render_report(
        merged, coverage, sp_summary, type_conf, fl_conf, hint_conf,
        boundary, score_recall, conf_corr, min_overlap,
    )
    (output_dir / "report.md").write_text(report)
    logger.success(f"Report written to {output_dir / 'report.md'}")


if __name__ == "__main__":
    app()
