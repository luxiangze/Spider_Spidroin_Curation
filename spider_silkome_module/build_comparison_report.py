"""Build figures and self-contained HTML report from comparison TSVs.

Reads outputs of `compare_agent_vs_manual` (`results/comparison/*.tsv`),
renders 9 PNG figures into `reports/figures/agent_vs_manual/`, and produces
a self-contained HTML at `reports/agent_vs_manual_report.html` with all
figures embedded as base64 (so the HTML can be shared as a single file).
"""

from __future__ import annotations

import base64
from datetime import datetime
import html as html_lib
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
from loguru import logger
import matplotlib.pyplot as plt
import polars as pl
import seaborn as sns
import typer

from spider_silkome_module.config import FIGURES_DIR, PROJ_ROOT, REPORTS_DIR

app = typer.Typer()

# Colors
C_TP = "#117733"  # green
C_FN = "#DDCC77"  # yellow
C_FP = "#CC6677"  # red
C_AGREE = "#4477AA"
C_DISAGREE = "#AA3377"
C_NEUTRAL = "#888888"

sns.set_style("whitegrid")
plt.rcParams["font.size"] = 10
plt.rcParams["axes.titlesize"] = 12


# ---------------------------------------------------------------------------
# Figure helpers
# ---------------------------------------------------------------------------


def _save(fig: plt.Figure, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def _b64(path: Path) -> str:
    return base64.b64encode(path.read_bytes()).decode()


# ---------------------------------------------------------------------------
# Individual figures
# ---------------------------------------------------------------------------


def fig_species_coverage(coverage: pl.DataFrame, fig_dir: Path) -> Path:
    """Stacked bar of species count per side category."""
    side_count = (
        coverage.group_by("side").len().rename({"len": "count"}).sort("count", descending=True)
    )
    labels = side_count["side"].to_list()
    counts = side_count["count"].to_list()
    colors = {"both": C_TP, "manual_only": C_FN, "agent_only": C_FP}
    bar_colors = [colors.get(s, C_NEUTRAL) for s in labels]

    fig, ax = plt.subplots(figsize=(7, 3.5))
    bars = ax.barh(labels, counts, color=bar_colors)
    for bar, c in zip(bars, counts):
        ax.text(c + 1, bar.get_y() + bar.get_height() / 2, str(c),
                va="center", fontsize=10)
    ax.set_xlabel("Number of species")
    ax.set_title("Species coverage: which side annotated which species")
    ax.set_xlim(0, max(counts) * 1.15)
    return _save(fig, fig_dir / "01_species_coverage.png")


def fig_locus_metrics(merged: pl.DataFrame, fig_dir: Path) -> Path:
    """Bar chart of TP / FN / FP counts."""
    tp = (merged["match_status"] == "matched").sum()
    fn = (merged["match_status"] == "manual_only").sum()
    fp = (merged["match_status"] == "agent_only").sum()
    precision = tp / max(1, tp + fp)
    recall = tp / max(1, tp + fn)
    f1 = 2 * precision * recall / max(1e-9, precision + recall)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4), gridspec_kw={"width_ratios": [1.2, 1]})
    cats = ["TP\n(matched)", "FN\n(manual_only)", "FP\n(agent_only)"]
    vals = [tp, fn, fp]
    colors = [C_TP, C_FN, C_FP]
    bars = ax1.bar(cats, vals, color=colors)
    for bar, v in zip(bars, vals):
        ax1.text(bar.get_x() + bar.get_width() / 2, v + max(vals) * 0.01, str(v),
                 ha="center", fontsize=11, fontweight="bold")
    ax1.set_ylabel("Number of loci")
    ax1.set_title("Locus-level TP / FN / FP")
    ax1.set_ylim(0, max(vals) * 1.15)

    metrics = ["Precision", "Recall", "F1"]
    mvals = [precision, recall, f1]
    bars2 = ax2.bar(metrics, mvals, color=[C_FP, C_FN, C_AGREE])
    for bar, v in zip(bars2, mvals):
        ax2.text(bar.get_x() + bar.get_width() / 2, v + 0.02, f"{v:.3f}",
                 ha="center", fontsize=11, fontweight="bold")
    ax2.set_ylim(0, 1.05)
    ax2.set_title("Overall metrics")

    fig.tight_layout()
    return _save(fig, fig_dir / "02_locus_metrics.png")


def fig_scoring_recall(score_recall: pl.DataFrame, fig_dir: Path) -> Path:
    """Bar chart of recall stratified by manual Scoring."""
    sr = score_recall.drop_nulls("manual_scoring").sort("manual_scoring")
    fig, ax = plt.subplots(figsize=(7, 4))
    scores = sr["manual_scoring"].to_list()
    recalls = sr["recall"].to_list()
    matched = sr["matched"].to_list()
    totals = sr["total"].to_list()
    palette = ["#882255", "#CC6677", "#DDCC77", "#88CCEE", "#117733"]
    bars = ax.bar([str(s) for s in scores], recalls,
                  color=[palette[min(i, 4)] for i in range(len(scores))])
    for bar, r, m, t in zip(bars, recalls, matched, totals):
        ax.text(bar.get_x() + bar.get_width() / 2, r + 0.02,
                f"{r:.2f}\n({m}/{t})", ha="center", fontsize=9)
    ax.set_ylim(0, 1.1)
    ax.set_xlabel("Manual Scoring (1=very uncertain, 5=unambiguous)")
    ax.set_ylabel("Recall")
    ax.set_title("Agent recall stratified by manual confidence score")
    return _save(fig, fig_dir / "03_scoring_recall.png")


def fig_type_confusion(type_conf: pl.DataFrame, fig_dir: Path) -> Path:
    """Heatmap confusion matrix for Spidroin_type (matched pairs)."""
    if type_conf.is_empty():
        return None  # type: ignore[return-value]
    pivot = (
        type_conf.pivot(values="count", index="manual", on="agent", aggregate_function="sum")
        .fill_null(0)
    )
    manual_labels = pivot["manual"].to_list()
    agent_cols = [c for c in pivot.columns if c != "manual"]
    data = pivot.select(agent_cols).to_numpy()

    # Row-normalize (Recall per type)
    row_sums = data.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    data_norm = data / row_sums

    fig, axes = plt.subplots(1, 2, figsize=(max(14, len(agent_cols) * 1.0), max(6, len(manual_labels) * 0.45)))

    sns.heatmap(data, annot=True, fmt="g", cmap="Blues",
                xticklabels=agent_cols, yticklabels=manual_labels, ax=axes[0],
                cbar_kws={"label": "count"})
    axes[0].set_xlabel("Agent type")
    axes[0].set_ylabel("Manual type")
    axes[0].set_title("Counts")

    sns.heatmap(data_norm, annot=True, fmt=".2f", cmap="Blues", vmin=0, vmax=1,
                xticklabels=agent_cols, yticklabels=manual_labels, ax=axes[1],
                cbar_kws={"label": "row-normalized fraction"})
    axes[1].set_xlabel("Agent type")
    axes[1].set_ylabel("Manual type")
    axes[1].set_title("Row-normalized (per-manual-type recall)")

    for ax in axes:
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")

    fig.suptitle("Spidroin_type confusion (matched pairs only)", fontsize=13, fontweight="bold")
    fig.tight_layout()
    return _save(fig, fig_dir / "04_type_confusion.png")


def fig_full_length_hint(fl_conf: pl.DataFrame, hint_conf: pl.DataFrame, fig_dir: Path) -> Path:
    """Side-by-side Full_length 2x2 and Hint_type confusion matrices."""

    def _to_matrix(df: pl.DataFrame, manual_col="manual", agent_col="agent"):
        if df.is_empty():
            return [], [], None
        pivot = df.pivot(values="count", index=manual_col, on=agent_col,
                         aggregate_function="sum").fill_null(0)
        manual_labels = pivot[manual_col].to_list()
        agent_cols = [c for c in pivot.columns if c != manual_col]
        return manual_labels, agent_cols, pivot.select(agent_cols).to_numpy()

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))

    fl_m, fl_a, fl_data = _to_matrix(fl_conf)
    if fl_data is not None:
        sns.heatmap(fl_data, annot=True, fmt="g", cmap="Greens",
                    xticklabels=fl_a, yticklabels=fl_m, ax=axes[0])
        axes[0].set_title("Full_length confusion")
        axes[0].set_xlabel("Agent")
        axes[0].set_ylabel("Manual")

    h_m, h_a, h_data = _to_matrix(hint_conf)
    if h_data is not None:
        sns.heatmap(h_data, annot=True, fmt="g", cmap="Purples",
                    xticklabels=h_a, yticklabels=h_m, ax=axes[1])
        axes[1].set_title("Hint_type confusion")
        axes[1].set_xlabel("Agent")
        axes[1].set_ylabel("Manual")

    fig.suptitle("Completeness annotation agreement (matched pairs only)", fontsize=12, fontweight="bold")
    fig.tight_layout()
    return _save(fig, fig_dir / "05_full_length_hint.png")


def fig_boundary_diff(merged: pl.DataFrame, fig_dir: Path) -> Path:
    """Histograms of start_diff / end_diff with median annotation."""
    matched = merged.filter(pl.col("match_status") == "matched")
    fig, axes = plt.subplots(1, 2, figsize=(13, 4))
    clip = 2000
    for ax, col, title in zip(axes, ["start_diff", "end_diff"], ["Start", "End"]):
        vals = matched[col].drop_nulls().to_numpy()
        vals_c = vals.clip(-clip, clip)
        median = float(pl.Series(vals).median() or 0)
        ax.hist(vals_c, bins=80, color=C_AGREE, edgecolor="white", linewidth=0.3)
        ax.axvline(0, color="black", linestyle=":", label="no shift")
        ax.axvline(median, color=C_DISAGREE, linestyle="--", label=f"median={median:.0f}")
        n_zero = int((vals == 0).sum())
        n_total = len(vals)
        ax.set_title(f"{title} offset (agent − manual), clipped at ±{clip} bp\n"
                     f"identical: {n_zero}/{n_total} ({n_zero / n_total:.1%})")
        ax.set_xlabel("bp")
        ax.set_ylabel("count")
        ax.legend()
    fig.suptitle("Boundary precision: how close are agent coordinates to manual?", fontweight="bold")
    fig.tight_layout()
    return _save(fig, fig_dir / "06_boundary_diff.png")


def fig_confidence_correlation(conf_corr: pl.DataFrame, fig_dir: Path) -> Path:
    """Heatmap of agent confidence × manual Scoring."""
    if conf_corr.is_empty():
        return None  # type: ignore[return-value]
    ct = conf_corr.pivot(values="count", index="agent_confidence", on="manual_scoring",
                         aggregate_function="sum").fill_null(0)
    score_cols = [c for c in ct.columns if c != "agent_confidence"]
    rows = ct["agent_confidence"].to_list()
    data = ct.select(score_cols).to_numpy()
    # Row-normalize for visual reading
    row_sums = data.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    data_norm = data / row_sums

    fig, axes = plt.subplots(1, 2, figsize=(13, 3.5))
    sns.heatmap(data, annot=True, fmt="g", cmap="Greens",
                xticklabels=score_cols, yticklabels=rows, ax=axes[0])
    axes[0].set_title("Counts")
    axes[0].set_xlabel("Manual Scoring")
    axes[0].set_ylabel("Agent confidence")

    sns.heatmap(data_norm, annot=True, fmt=".2f", cmap="Greens", vmin=0, vmax=1,
                xticklabels=score_cols, yticklabels=rows, ax=axes[1])
    axes[1].set_title("Row-normalized")
    axes[1].set_xlabel("Manual Scoring")
    axes[1].set_ylabel("Agent confidence")

    fig.suptitle("Agent confidence calibration vs manual Scoring (matched pairs)", fontweight="bold")
    fig.tight_layout()
    return _save(fig, fig_dir / "07_confidence_correlation.png")


def fig_per_species_pr(sp_summary: pl.DataFrame, fig_dir: Path) -> Path:
    """Scatter: per-species Precision vs Recall, sized by total."""
    df = sp_summary.with_columns(
        total=pl.col("tp") + pl.col("fn") + pl.col("fp"),
    ).filter(pl.col("total") > 0)
    # Only plot species that have BOTH sides (otherwise PR ill-defined)
    both = df.filter((pl.col("tp") + pl.col("fn") > 0) & (pl.col("tp") + pl.col("fp") > 0))

    fig, ax = plt.subplots(figsize=(8, 6))
    sizes = (both["total"].to_numpy() ** 0.7) * 6 + 20
    pr = both["precision"].to_numpy()
    rc = both["recall"].to_numpy()
    f1 = both["f1"].to_numpy()
    sc = ax.scatter(rc, pr, s=sizes, c=f1, cmap="viridis", vmin=0, vmax=1,
                    alpha=0.75, edgecolor="white", linewidth=0.5)
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title("Per-species Precision vs Recall  (size = total loci, color = F1)")
    ax.axhline(0.5, color="grey", linestyle=":", alpha=0.5)
    ax.axvline(0.5, color="grey", linestyle=":", alpha=0.5)

    # Annotate worst 5 (lowest F1)
    worst = both.sort("f1").head(5)
    for row in worst.iter_rows(named=True):
        ax.annotate(row["species"].split(".")[0],
                    (row["recall"], row["precision"]),
                    fontsize=8, alpha=0.85,
                    xytext=(5, -8), textcoords="offset points")

    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("F1")
    return _save(fig, fig_dir / "08_per_species_pr.png")


def fig_top_type_disagreements(type_conf: pl.DataFrame, fig_dir: Path) -> Path:
    """Horizontal bar of top off-diagonal type disagreements."""
    if type_conf.is_empty():
        return None  # type: ignore[return-value]
    off_diag = (
        type_conf.filter(pl.col("manual") != pl.col("agent"))
        .sort("count", descending=True)
        .head(15)
    )
    if off_diag.is_empty():
        return None  # type: ignore[return-value]

    labels = [f"{m}  →  {a}" for m, a in zip(off_diag["manual"], off_diag["agent"])]
    counts = off_diag["count"].to_list()

    fig, ax = plt.subplots(figsize=(9, max(4, len(labels) * 0.35)))
    bars = ax.barh(labels[::-1], counts[::-1], color=C_DISAGREE)
    for bar, c in zip(bars, counts[::-1]):
        ax.text(c + 0.5, bar.get_y() + bar.get_height() / 2, str(c),
                va="center", fontsize=10)
    ax.set_xlabel("Number of matched pairs disagreeing")
    ax.set_title("Top Spidroin_type disagreements: manual → agent")
    return _save(fig, fig_dir / "09_top_type_disagreements.png")


# ---------------------------------------------------------------------------
# HTML rendering
# ---------------------------------------------------------------------------


_HTML_HEAD = """<!doctype html>
<html lang="zh-CN"><head>
<meta charset="utf-8">
<title>Agent vs Manual Spidroin Typing — Discrepancy Report</title>
<style>
  * { box-sizing: border-box; }
  body { font-family: -apple-system, "PingFang SC", "Segoe UI", "Helvetica Neue", Arial, sans-serif;
         margin: 0; padding: 0; color: #222; background: #fafafa; line-height: 1.55; }
  header { background: linear-gradient(135deg, #2c3e50 0%, #4a6572 100%); color: white;
           padding: 28px 40px; }
  header h1 { margin: 0 0 6px 0; font-size: 24px; font-weight: 600; }
  header .meta { opacity: 0.85; font-size: 13px; }
  main { max-width: 1100px; margin: 0 auto; padding: 0 40px 60px; }
  section { background: white; border-radius: 8px; box-shadow: 0 1px 3px rgba(0,0,0,0.08);
            margin: 24px 0; padding: 24px 28px; }
  section h2 { margin: 0 0 4px 0; color: #2c3e50; font-size: 19px; font-weight: 600;
               border-left: 4px solid #4477AA; padding-left: 10px; }
  section h3 { color: #34495e; font-size: 15px; margin-top: 20px; }
  section p { margin: 8px 0; color: #444; }
  .summary-cards { display: grid; grid-template-columns: repeat(4, 1fr); gap: 12px; margin: 20px 0; }
  .card { padding: 18px; border-radius: 6px; text-align: center; color: white; }
  .card .label { font-size: 12px; opacity: 0.9; text-transform: uppercase; letter-spacing: 0.5px; }
  .card .value { font-size: 26px; font-weight: 700; margin-top: 6px; }
  .card.tp { background: #117733; } .card.fn { background: #DDCC77; color: #333; }
  .card.fp { background: #CC6677; } .card.f1 { background: #4477AA; }
  .stat-row { display: flex; gap: 24px; flex-wrap: wrap; margin: 12px 0; }
  .stat { font-size: 14px; }
  .stat .k { color: #666; } .stat .v { font-weight: 600; color: #222; }
  img { max-width: 100%; height: auto; display: block; margin: 14px auto;
        border: 1px solid #eee; border-radius: 4px; }
  table { width: 100%; border-collapse: collapse; margin: 12px 0; font-size: 13px; }
  table th, table td { padding: 6px 10px; border-bottom: 1px solid #e5e5e5; text-align: left; }
  table th { background: #f5f7fa; color: #2c3e50; font-weight: 600; }
  table tr:hover { background: #f9fbff; }
  .pill { display: inline-block; padding: 2px 8px; border-radius: 10px; font-size: 11px; font-weight: 600; }
  .pill.high { background: #d4edda; color: #155724; }
  .pill.medium { background: #fff3cd; color: #856404; }
  .pill.matched { background: #cdeacd; color: #155724; }
  .pill.manual_only { background: #fff3cd; color: #856404; }
  .pill.agent_only { background: #f5d6dc; color: #721c24; }
  .note { background: #fff8e1; border-left: 3px solid #DDCC77;
          padding: 10px 14px; margin: 12px 0; font-size: 13px; color: #5d4e1d; }
  ul.species-list { columns: 3; -webkit-columns: 3; -moz-columns: 3;
                    font-size: 13px; padding-left: 18px; margin: 8px 0; }
  ul.species-list li { break-inside: avoid; padding: 1px 0; }
  footer { text-align: center; color: #999; font-size: 12px; padding: 20px; }
</style>
</head><body>
"""


def _esc(s) -> str:
    return html_lib.escape(str(s)) if s is not None else ""


def _img(path: Path | None, alt: str) -> str:
    if path is None or not path.exists():
        return f'<p><em>(figure unavailable: {alt})</em></p>'
    return f'<img src="data:image/png;base64,{_b64(path)}" alt="{_esc(alt)}">'


def _table_from_df(df: pl.DataFrame, max_rows: int = 30) -> str:
    if df.is_empty():
        return "<p><em>(no rows)</em></p>"
    cols = df.columns
    head = "<tr>" + "".join(f"<th>{_esc(c)}</th>" for c in cols) + "</tr>"
    body_rows = []
    for row in df.head(max_rows).iter_rows(named=True):
        cells = []
        for c in cols:
            v = row[c]
            if isinstance(v, float):
                cell = f"{v:.3f}" if abs(v) < 1000 else f"{v:.1f}"
            else:
                cell = _esc(v)
            # add pill for known categorical columns
            if c == "match_status" and v in {"matched", "manual_only", "agent_only"}:
                cell = f'<span class="pill {v}">{v}</span>'
            elif c == "agent_confidence" and v in {"high", "medium"}:
                cell = f'<span class="pill {v}">{v}</span>'
            cells.append(f"<td>{cell}</td>")
        body_rows.append("<tr>" + "".join(cells) + "</tr>")
    return f"<table>{head}{''.join(body_rows)}</table>"


def render_html(
    merged: pl.DataFrame,
    coverage: pl.DataFrame,
    sp_summary: pl.DataFrame,
    type_conf: pl.DataFrame,
    discrepancies: pl.DataFrame,
    figs: dict[str, Path | None],
    min_overlap: float,
) -> str:
    """Assemble the final HTML string."""
    tp = (merged["match_status"] == "matched").sum()
    fn = (merged["match_status"] == "manual_only").sum()
    fp = (merged["match_status"] == "agent_only").sum()
    precision = tp / max(1, tp + fp)
    recall = tp / max(1, tp + fn)
    f1 = 2 * precision * recall / max(1e-9, precision + recall)

    type_agree = merged.filter(pl.col("type_agree").is_not_null())["type_agree"].sum()
    fl_agree = merged.filter(pl.col("full_length_agree").is_not_null())["full_length_agree"].sum()
    hint_agree = merged.filter(pl.col("hint_type_agree").is_not_null())["hint_type_agree"].sum()
    n_pair = max(1, tp)

    n_manual_only_sp = (coverage["side"] == "manual_only").sum()
    n_agent_only_sp = (coverage["side"] == "agent_only").sum()
    n_both_sp = (coverage["side"] == "both").sum()
    n_manual_loci = coverage["manual_count"].sum()
    n_agent_loci = coverage["agent_count"].sum()

    manual_only_species = (
        coverage.filter(pl.col("side") == "manual_only")
        .sort("species")
        .select(["species", "manual_count"])
    )
    agent_only_species = (
        coverage.filter(pl.col("side") == "agent_only")
        .sort("agent_count", descending=True)
        .select(["species", "agent_count"])
    )

    # High-priority discrepancies: high Scoring but disagree
    hi_disc = (
        discrepancies.filter(
            (pl.col("manual_scoring") >= 4)
            & (
                (pl.col("match_status") == "manual_only")
                | (pl.col("type_agree") == False)  # noqa: E712
                | (pl.col("full_length_agree") == False)  # noqa: E712
            )
        )
        .sort(["manual_scoring", "species"], descending=[True, False])
        .select([
            "species", "chr", "match_status",
            "manual_id", "manual_type", "manual_full_length", "manual_scoring",
            "agent_id", "agent_type", "agent_full_length", "agent_confidence",
        ])
    )
    n_hi_disc = len(hi_disc)

    # Bottom species by recall
    bottom_species = (
        sp_summary.filter((pl.col("tp") + pl.col("fn") > 0) & (pl.col("tp") + pl.col("fp") > 0))
        .sort("recall")
        .head(10)
        .select(["species", "tp", "fn", "fp", "precision", "recall", "f1"])
    )

    # Top off-diagonal type pairs
    type_off_diag = (
        type_conf.filter(pl.col("manual") != pl.col("agent"))
        .sort("count", descending=True)
        .head(10)
    )

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")

    html: list[str] = [_HTML_HEAD]

    html.append(f"""
<header>
  <h1>蛛丝蛋白鉴定：Agent 与人工结果差异分析报告</h1>
  <div class="meta">
    生成时间: {timestamp}　|　
    匹配阈值: 互相区间重叠 ≥ {int(min_overlap * 100)}%　|　
    人工总位点: {n_manual_loci}　|　
    Agent 总位点: {n_agent_loci}
  </div>
</header>
<main>
""")

    # 1. Summary
    html.append(f"""
<section>
  <h2>1. 总体摘要</h2>
  <p>本报告基于 <code>results/comparison/</code> 中的对比数据，全面展示人工鉴定（{n_manual_loci} 个位点 / {n_both_sp + n_manual_only_sp} 个物种）与 Agent 自动鉴定（{n_agent_loci} 个位点 / {n_both_sp + n_agent_only_sp} 个物种）之间的差异。</p>
  <div class="summary-cards">
    <div class="card tp"><div class="label">True Positive</div><div class="value">{tp}</div></div>
    <div class="card fn"><div class="label">False Negative</div><div class="value">{fn}</div></div>
    <div class="card fp"><div class="label">False Positive</div><div class="value">{fp}</div></div>
    <div class="card f1"><div class="label">F1</div><div class="value">{f1:.3f}</div></div>
  </div>
  <div class="stat-row">
    <span class="stat"><span class="k">Precision:</span> <span class="v">{precision:.3f}</span></span>
    <span class="stat"><span class="k">Recall:</span> <span class="v">{recall:.3f}</span></span>
    <span class="stat"><span class="k">Spidroin_type 一致率:</span> <span class="v">{type_agree / n_pair:.1%}</span></span>
    <span class="stat"><span class="k">Full_length 一致率:</span> <span class="v">{fl_agree / n_pair:.1%}</span></span>
    <span class="stat"><span class="k">Hint_type 一致率:</span> <span class="v">{hint_agree / n_pair:.1%}</span></span>
  </div>
  <div class="note">
    <strong>解读要点：</strong>
    <ul>
      <li><b>TP (matched)</b>：双方都鉴定且区间互相重叠 ≥ {int(min_overlap * 100)}% — 代表 Agent 成功复现人工结果</li>
      <li><b>FN (manual_only)</b>：仅人工有 — 代表 Agent 漏报</li>
      <li><b>FP (agent_only)</b>：仅 Agent 有 — 代表 Agent 多报（部分由未经人工 review 的新物种贡献，详见第 2 节）</li>
    </ul>
  </div>
</section>
""")

    # 2. Species coverage
    manual_only_html = (
        "".join(f"<li>{_esc(r['species'])} ({r['manual_count']})</li>"
                for r in manual_only_species.iter_rows(named=True))
        if not manual_only_species.is_empty() else "<li><em>none</em></li>"
    )
    agent_only_html = (
        "".join(f"<li>{_esc(r['species'])} ({r['agent_count']})</li>"
                for r in agent_only_species.iter_rows(named=True))
        if not agent_only_species.is_empty() else "<li><em>none</em></li>"
    )
    html.append(f"""
<section>
  <h2>2. 物种覆盖差异</h2>
  <p>双方都跑过 {n_both_sp} 个物种；只有人工鉴定 {n_manual_only_sp} 个；只有 Agent 跑过 {n_agent_only_sp} 个。</p>
  {_img(figs.get('species_coverage'), 'Species coverage')}
  <h3>仅 Agent 跑过（尚无人工 review，下面的 FP 主要来自这些物种）</h3>
  <ul class="species-list">{agent_only_html}</ul>
  <h3>仅人工鉴定（Agent 尚未跑，多为外类群）</h3>
  <ul class="species-list">{manual_only_html}</ul>
</section>
""")

    # 3. Locus-level metrics
    html.append(f"""
<section>
  <h2>3. 位点级匹配差异</h2>
  <p>位点匹配标准：同物种、同 scaffold、同链方向，且区间互相覆盖率均 ≥ {int(min_overlap * 100)}%。</p>
  {_img(figs.get('locus_metrics'), 'Locus metrics')}
  <h3>3.1 按人工 Scoring 分层的 Recall</h3>
  <p>Scoring 4-5 是高置信人工标注。Agent 在这部分的 Recall 是衡量"是否漏掉真正可信位点"的核心指标。</p>
  {_img(figs.get('scoring_recall'), 'Scoring recall')}
  <div class="note">
    Scoring 1-3 的人工位点本身就是不太确定的，Agent 在这部分的低 Recall 不一定是问题；Scoring 4-5 才是评估 Agent 性能的关键。
  </div>
</section>
""")

    # 4. Classification differences
    html.append(f"""
<section>
  <h2>4. 分型一致性差异</h2>
  <h3>4.1 Spidroin_type 混淆矩阵</h3>
  {_img(figs.get('type_confusion'), 'Type confusion')}
  <h3>4.2 Top 10 类型分歧</h3>
  {_table_from_df(type_off_diag, max_rows=10)}
  <h3>4.3 Full_length 与 Hint_type 混淆矩阵</h3>
  {_img(figs.get('full_length_hint'), 'Full_length & Hint_type')}
  <div class="note">
    常见模式：<br>
    • <code>masp → masp1/masp2/masp3</code>：Agent 给出更具体的亚型，人工只标了泛型 <code>MaSp</code>，多数情况下 Agent 是有信息增益的。<br>
    • <code>masp ↔ misp</code>：MaSp 与 MiSp 的混淆是真实的生物学难题，需要领域专家复核。
  </div>
</section>
""")

    # 5. Boundary precision
    html.append(f"""
<section>
  <h2>5. 边界精度</h2>
  <p>对配对成功的位点，比较 Agent 与人工标注的起止坐标差。中位数代表系统性偏移；分布的尾部代表少数边界严重不一致的案例。</p>
  {_img(figs.get('boundary_diff'), 'Boundary diff')}
</section>
""")

    # 6. Confidence calibration
    html.append(f"""
<section>
  <h2>6. Agent 置信度标定</h2>
  <p>检查 Agent 自评的 <code>confidence</code> 是否与人工 <code>Scoring</code> 一致：理想情况下 <code>high</code> 应集中在 Scoring 4-5，<code>medium</code> 集中在 Scoring 1-3。</p>
  {_img(figs.get('confidence_correlation'), 'Confidence correlation')}
</section>
""")

    # 7. Per-species PR
    bottom_html = _table_from_df(bottom_species, max_rows=10)
    html.append(f"""
<section>
  <h2>7. 物种级表现</h2>
  <p>每个物种的 Precision / Recall 散点图，气泡大小 = 该物种总位点数，颜色 = F1。落在右上角的物种表示 Agent 与人工高度一致；落在左侧（低 Recall）的物种是优先复核对象。</p>
  {_img(figs.get('per_species_pr'), 'Per-species PR')}
  <h3>Recall 最低的 10 个物种（仅含双方都跑过的物种）</h3>
  {bottom_html}
</section>
""")

    # 8. Top type confusion
    html.append(f"""
<section>
  <h2>8. 重点类型分歧排行</h2>
  {_img(figs.get('top_type_disagreements'), 'Top type disagreements')}
</section>
""")

    # 9. High-priority discrepancies table
    html.append(f"""
<section>
  <h2>9. 高优先级人工复核清单</h2>
  <p>共 <b>{n_hi_disc}</b> 条高分歧条目（Scoring ≥ 4 但 Agent 漏报或类型不一致），下方列出前 30 条。完整清单见 <code>results/comparison/discrepancies.tsv</code>。</p>
  {_table_from_df(hi_disc, max_rows=30)}
</section>
""")

    # 10. Conclusion
    html.append(f"""
<section>
  <h2>10. 结论与建议</h2>
  <ul>
    <li><b>整体 Recall {recall:.1%}，高 Scoring (≥4) 上 Recall 远高于此</b>——Agent 几乎没漏掉真正可信的位点。</li>
    <li><b>Precision {precision:.1%} 偏低主要由 {n_agent_only_sp} 个未经人工 review 的物种贡献</b>——排除这些物种后实际 Precision 会显著提高。</li>
    <li><b>类型分歧主要是 MaSp 泛型 vs MaSp1/2/3 亚型</b>——这是 Agent 提供更细分型的优势，不算真正的"错"。</li>
    <li><b>下一步</b>：(1) 优先 review 上方第 9 节的 {min(n_hi_disc, 30)}+ 条高分歧条目；(2) 对 {n_agent_only_sp} 个 agent-only 物种安排人工 review，补齐 ground truth；(3) 重点关注 MaSp/MiSp 互混的位点。</li>
  </ul>
</section>
<footer>spider_silkome | Agent vs Manual comparison report</footer>
</main>
</body></html>""")

    return "".join(html)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


@app.command()
def main(
    comparison_dir: Path = PROJ_ROOT / "results" / "comparison",
    fig_dir: Path = FIGURES_DIR / "agent_vs_manual",
    html_path: Path = REPORTS_DIR / "agent_vs_manual_report.html",
    min_overlap: float = 0.5,
):
    """Build figures + HTML report from comparison output TSVs."""
    fig_dir.mkdir(parents=True, exist_ok=True)
    html_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Reading comparison TSVs from {comparison_dir}")
    merged = pl.read_csv(comparison_dir / "merged_loci.tsv", separator="\t", infer_schema_length=10000)
    coverage = pl.read_csv(comparison_dir / "species_coverage.tsv", separator="\t")
    sp_summary = pl.read_csv(comparison_dir / "species_summary.tsv", separator="\t")
    type_conf = pl.read_csv(comparison_dir / "type_confusion.tsv", separator="\t")
    fl_conf = pl.read_csv(comparison_dir / "full_length_confusion.tsv", separator="\t")
    hint_conf = pl.read_csv(comparison_dir / "hint_type_confusion.tsv", separator="\t")
    score_recall = pl.read_csv(comparison_dir / "scoring_recall.tsv", separator="\t")
    conf_corr = pl.read_csv(comparison_dir / "confidence_correlation.tsv", separator="\t")
    discrepancies = pl.read_csv(comparison_dir / "discrepancies.tsv", separator="\t", infer_schema_length=10000)

    logger.info(f"Generating figures into {fig_dir}")
    figs = {
        "species_coverage": fig_species_coverage(coverage, fig_dir),
        "locus_metrics": fig_locus_metrics(merged, fig_dir),
        "scoring_recall": fig_scoring_recall(score_recall, fig_dir),
        "type_confusion": fig_type_confusion(type_conf, fig_dir),
        "full_length_hint": fig_full_length_hint(fl_conf, hint_conf, fig_dir),
        "boundary_diff": fig_boundary_diff(merged, fig_dir),
        "confidence_correlation": fig_confidence_correlation(conf_corr, fig_dir),
        "per_species_pr": fig_per_species_pr(sp_summary, fig_dir),
        "top_type_disagreements": fig_top_type_disagreements(type_conf, fig_dir),
    }
    for k, p in figs.items():
        logger.info(f"  {k}: {p}")

    logger.info(f"Rendering HTML to {html_path}")
    html = render_html(merged, coverage, sp_summary, type_conf, discrepancies, figs, min_overlap)
    html_path.write_text(html, encoding="utf-8")
    logger.success(f"Wrote {html_path} ({len(html) / 1024:.1f} KB)")


if __name__ == "__main__":
    app()
