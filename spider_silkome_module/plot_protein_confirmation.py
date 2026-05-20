"""
plot_protein_confirmation.py — five publication-ready figures summarising
the protein confirmation pipeline.  Each subcommand writes both PNG and SVG
with a consistent CCDS-style colour palette.

Subcommands
-----------
  funnel                  Stage-by-stage retention of spidroin candidates
  score-by-type           Validation score distribution per type (boxplot)
  identity-coverage       Pairwise identity vs coverage scatter
  species-type-heatmap    Validated count per species × type
  gene-structure          Representative gene tracks per type
  all                     Run every subcommand
"""

from pathlib import Path

from loguru import logger
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import polars as pl
import typer

from spider_silkome_module.config import FIGURES_DIR, PROCESSED_DATA_DIR

app = typer.Typer()

DEFAULT_TASK_NAME = "protein_confirmation"

# Consistent palette across figures
_STATUS_COLORS = {
    "validated": "#4CAF50",
    "partial":   "#FFC107",
    "failed":    "#9E9E9E",
}
plt.rcParams.update({
    "font.size": 11,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


# ── data loading ─────────────────────────────────────────────────────────────


def load_summary(summary_tsv: Path) -> pl.DataFrame:
    """Cross-species summary written by the notebook's Step 8."""
    if not summary_tsv.exists():
        logger.error(f"Summary file not found: {summary_tsv}")
        raise typer.Exit(1)
    return pl.read_csv(summary_tsv, separator="\t", infer_schema_length=2000)


def _save(fig, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path.with_suffix(".png"), dpi=200, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


# ── 1) funnel ────────────────────────────────────────────────────────────────


@app.command()
def funnel(
    summary_tsv: Path = PROCESSED_DATA_DIR / DEFAULT_TASK_NAME / "summary_all_species.tsv",
    output_dir: Path = FIGURES_DIR / DEFAULT_TASK_NAME,
):
    """Horizontal bar funnel: how many candidates survive each stage."""
    df = load_summary(summary_tsv)

    n_total = df.height
    n_translated = df.filter(pl.col("predicted_protein_length").fill_null(0) > 0).height
    n_pairwise = df.filter(pl.col("pairwise_identity").fill_null(0) >= 0.5).height
    n_diamond = df.filter(pl.col("diamond_evalue").fill_null(1.0) <= 1e-10).height
    n_hmm = df.filter(
        (pl.col("hmm_ntd_evalue").fill_null(1.0) <= 1e-10)
        | (pl.col("hmm_ctd_evalue").fill_null(1.0) <= 1e-10)
    ).height
    n_validated = df.filter(pl.col("validation_status") == "validated").height
    n_partial = df.filter(pl.col("validation_status") == "partial").height

    stages = [
        ("Loci (total)",          n_total),
        ("Translated",            n_translated),
        ("Pairwise ≥50%",         n_pairwise),
        ("DIAMOND e≤1e-10",       n_diamond),
        ("HMM hit (NTD or CTD)",  n_hmm),
        ("Partial+Validated",     n_partial + n_validated),
        ("Validated (≥80)",       n_validated),
    ]
    labels = [s[0] for s in stages][::-1]
    counts = [s[1] for s in stages][::-1]
    pct = [c / max(1, n_total) * 100 for c in counts]

    fig, ax = plt.subplots(figsize=(8, 4.5))
    bars = ax.barh(labels, counts, color="#1f77b4", alpha=0.85)
    for bar, c, p in zip(bars, counts, pct):
        ax.text(bar.get_width() + max(counts) * 0.01,
                bar.get_y() + bar.get_height() / 2,
                f"{c}  ({p:.1f}%)", va="center", fontsize=10)
    ax.set_xlabel("Number of spidroin loci")
    ax.set_title("Pipeline funnel: spidroin loci surviving each filter")
    ax.set_xlim(0, max(counts) * 1.18)

    _save(fig, output_dir / "pipeline_funnel")
    logger.success(f"Saved {output_dir / 'pipeline_funnel.png'}")


# ── 2) score by type ─────────────────────────────────────────────────────────


@app.command()
def score_by_type(
    summary_tsv: Path = PROCESSED_DATA_DIR / DEFAULT_TASK_NAME / "summary_all_species.tsv",
    output_dir: Path = FIGURES_DIR / DEFAULT_TASK_NAME,
):
    """Boxplot of validation_score per spidroin type, ordered by median."""
    df = load_summary(summary_tsv)
    df = df.filter(pl.col("validation_score").is_not_null())
    if df.is_empty():
        logger.warning("No rows with validation_score; skip")
        return

    types = (df.group_by("spidroin_type")
               .agg(pl.col("validation_score").median().alias("med"))
               .sort("med", descending=True))
    type_order = types.get_column("spidroin_type").to_list()

    data = [df.filter(pl.col("spidroin_type") == t)
              .get_column("validation_score").to_numpy()
            for t in type_order]
    counts = [len(d) for d in data]

    fig, ax = plt.subplots(figsize=(max(7, len(type_order) * 0.6), 5))
    bp = ax.boxplot(data, tick_labels=type_order, patch_artist=True, widths=0.6,
                    medianprops=dict(color="black", linewidth=1.5))
    cmap = plt.get_cmap("tab20")
    for patch, color in zip(bp["boxes"], cmap.colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    for i, n in enumerate(counts, 1):
        ax.text(i, 102, f"n={n}", ha="center", fontsize=9, color="dimgrey")

    ax.axhline(80, color=_STATUS_COLORS["validated"], linestyle="--",
               linewidth=1, label="validated ≥80")
    ax.axhline(50, color=_STATUS_COLORS["partial"], linestyle="--",
               linewidth=1, label="partial ≥50")
    ax.set_ylabel("Validation score (0-100)")
    ax.set_xlabel("Spidroin type")
    ax.set_title("Validation score distribution per spidroin type")
    ax.set_ylim(-5, 110)
    ax.legend(loc="lower left", frameon=False)
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right")

    _save(fig, output_dir / "validation_score_by_type")
    logger.success(f"Saved {output_dir / 'validation_score_by_type.png'}")


# ── 3) identity vs coverage ──────────────────────────────────────────────────


@app.command()
def identity_coverage(
    summary_tsv: Path = PROCESSED_DATA_DIR / DEFAULT_TASK_NAME / "summary_all_species.tsv",
    output_dir: Path = FIGURES_DIR / DEFAULT_TASK_NAME,
):
    """Scatter plot of pairwise identity vs coverage, coloured by status."""
    df = load_summary(summary_tsv)
    df = df.filter(
        pl.col("pairwise_identity").is_not_null()
        & pl.col("pairwise_coverage_query").is_not_null()
    )
    if df.is_empty():
        logger.warning("No pairwise data; skip")
        return

    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    for status, color in _STATUS_COLORS.items():
        sub = df.filter(pl.col("validation_status") == status)
        if sub.is_empty():
            continue
        ax.scatter(
            sub.get_column("pairwise_identity").to_numpy(),
            sub.get_column("pairwise_coverage_query").to_numpy(),
            s=24, alpha=0.65, color=color, label=f"{status} (n={sub.height})",
        )
    ax.axvline(0.8, color="grey", linestyle=":", linewidth=1)
    ax.axhline(0.8, color="grey", linestyle=":", linewidth=1)
    ax.set_xlim(0, 1.05)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("Pairwise identity (predicted vs NCBI)")
    ax.set_ylabel("Coverage of NCBI reference")
    ax.set_title("Pairwise alignment: identity vs coverage")
    ax.legend(frameon=False, loc="lower right")

    _save(fig, output_dir / "identity_vs_coverage")
    logger.success(f"Saved {output_dir / 'identity_vs_coverage.png'}")


# ── 4) species × type heatmap ────────────────────────────────────────────────


@app.command()
def species_type_heatmap(
    summary_tsv: Path = PROCESSED_DATA_DIR / DEFAULT_TASK_NAME / "summary_all_species.tsv",
    output_dir: Path = FIGURES_DIR / DEFAULT_TASK_NAME,
    top_species: int = 30,
):
    """Heatmap of validated-spidroin counts per species × type."""
    df = load_summary(summary_tsv).filter(
        pl.col("validation_status") == "validated"
    )
    if df.is_empty():
        logger.warning("No validated spidroins for heatmap; skip")
        return

    pivot = (df.group_by(["species", "spidroin_type"])
               .agg(pl.len().alias("n"))
               .pivot(values="n", index="species", on="spidroin_type", aggregate_function="sum")
               .fill_null(0))
    species_order = (pivot.with_columns(
        pl.sum_horizontal(pl.exclude("species")).alias("_total"))
                          .sort("_total", descending=True)
                          .head(top_species)
                          .get_column("species").to_list())
    pivot = pivot.filter(pl.col("species").is_in(species_order))
    type_cols = [c for c in pivot.columns if c != "species"]
    matrix = pivot.select(type_cols).to_numpy()

    fig, ax = plt.subplots(figsize=(max(7, len(type_cols) * 0.6),
                                    max(4, len(species_order) * 0.28)))
    im = ax.imshow(matrix, aspect="auto", cmap="Blues")
    ax.set_xticks(range(len(type_cols)))
    ax.set_xticklabels(type_cols, rotation=30, ha="right")
    ax.set_yticks(range(len(species_order)))
    ax.set_yticklabels(species_order, fontsize=8)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            v = int(matrix[i, j])
            if v > 0:
                color = "white" if v > matrix.max() * 0.55 else "black"
                ax.text(j, i, str(v), ha="center", va="center",
                        color=color, fontsize=8)
    fig.colorbar(im, ax=ax, label="# validated spidroins", fraction=0.04)
    ax.set_title(f"Validated spidroins per species × type (top {len(species_order)})")
    ax.set_xlabel("Spidroin type")
    ax.set_ylabel("Species")

    _save(fig, output_dir / "species_type_heatmap")
    logger.success(f"Saved {output_dir / 'species_type_heatmap.png'}")


# ── 5) representative gene structures ────────────────────────────────────────


def _parse_fullprot_gff(gff_path: Path) -> dict[str, list[tuple[str, int, int]]]:
    """Return {spidroin_id: [(feature, start, end), ...]} for exon/intron rows."""
    out: dict[str, list[tuple[str, int, int]]] = {}
    if not gff_path.exists():
        return out
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 9 or parts[2] not in ("exon", "intron"):
                continue
            sid = parts[0]
            out.setdefault(sid, []).append((parts[2], int(parts[3]), int(parts[4])))
    return out


@app.command()
def gene_structure(
    summary_tsv: Path = PROCESSED_DATA_DIR / DEFAULT_TASK_NAME / "summary_all_species.tsv",
    confirmation_dir: Path = PROCESSED_DATA_DIR / DEFAULT_TASK_NAME / "confirmation",
    output_dir: Path = FIGURES_DIR / DEFAULT_TASK_NAME,
):
    """Pick one validated representative per type and draw its exon/intron tracks."""
    df = load_summary(summary_tsv).filter(
        (pl.col("validation_status") == "validated")
        & pl.col("validation_score").is_not_null()
    )
    if df.is_empty():
        logger.warning("No validated spidroins to draw; skip")
        return

    # Best-scoring representative per type.
    df = df.sort("validation_score", descending=True)
    reps = df.group_by("spidroin_type").first()
    reps = reps.sort("spidroin_type")

    rows = reps.to_dicts()
    fig, ax = plt.subplots(figsize=(10, max(3, len(rows) * 0.55)))

    max_len = 0
    for i, row in enumerate(rows):
        sid = row["spidroin_id"]
        sp = row["species"]
        sp_dir = confirmation_dir / sp
        feats = _parse_fullprot_gff(sp_dir / "fullprot_miniprot.gff").get(sid, [])
        if not feats:
            continue
        coords = [c for _, s, e in feats for c in (s, e)]
        gmin, gmax = min(coords), max(coords)
        max_len = max(max_len, gmax - gmin + 1)
        # Backbone line.
        ax.hlines(i, 0, gmax - gmin, color="#888", linewidth=1)
        for feat, s, e in feats:
            x0 = s - gmin
            w = max(1, e - s + 1)
            if feat == "exon":
                ax.add_patch(mpatches.Rectangle(
                    (x0, i - 0.18), w, 0.36, color="#1f77b4"
                ))
            elif feat == "intron":
                ax.add_patch(mpatches.Rectangle(
                    (x0, i - 0.04), w, 0.08, color="#d62728", alpha=0.6
                ))
        ax.text(-0.02 * max_len, i,
                f"{row['spidroin_type']}  ({sp[:20]})",
                ha="right", va="center", fontsize=9)

    ax.set_yticks([])
    ax.set_xlim(-0.05 * max_len, max_len * 1.02)
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.set_xlabel("Position on spidroin locus (bp, mRNA orientation)")
    ax.set_title("Representative gene structures per spidroin type")
    legend = [
        mpatches.Patch(color="#1f77b4", label="Exon"),
        mpatches.Patch(color="#d62728", alpha=0.6, label="Intron"),
    ]
    ax.legend(handles=legend, loc="upper right", frameon=False)

    _save(fig, output_dir / "gene_structure_examples")
    logger.success(f"Saved {output_dir / 'gene_structure_examples.png'}")


# ── meta ─────────────────────────────────────────────────────────────────────


@app.command(name="all")
def all_(
    summary_tsv: Path = PROCESSED_DATA_DIR / DEFAULT_TASK_NAME / "summary_all_species.tsv",
    confirmation_dir: Path = PROCESSED_DATA_DIR / DEFAULT_TASK_NAME / "confirmation",
    output_dir: Path = FIGURES_DIR / DEFAULT_TASK_NAME,
    top_species: int = 30,
):
    """Generate all five figures."""
    funnel(summary_tsv, output_dir)
    score_by_type(summary_tsv, output_dir)
    identity_coverage(summary_tsv, output_dir)
    species_type_heatmap(summary_tsv, output_dir, top_species)
    gene_structure(summary_tsv, confirmation_dir, output_dir)


if __name__ == "__main__":
    app()
