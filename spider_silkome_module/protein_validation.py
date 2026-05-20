"""
protein_validation.py — multi-evidence validation of predicted spidroin proteins.

Per species, this module ingests predicted_proteins.fa (from
protein_confirmation.py) and the typing_results TSV, and outputs:

  diamond_vs_silkome.tsv     # B. DIAMOND blastp vs silkome database
  hmmscan_domains.tsv        # C. hmmscan vs combined NTD/CTD HMM library
  motif_check.tsv            # D. signal-peptide hint + CTD cysteine motif
  protein_confirmation.tsv   # combined table with validation_score + status

The final TSV merges:
  - typing_results metadata (species, spidroin_type)
  - pairwise_vs_ncbi.tsv (dimension A)
  - dimension B / C / D / E rows
  - per-row validation_score (0-100) and status (validated / partial / failed)
"""

from pathlib import Path
import re

from Bio import SeqIO
from loguru import logger
import polars as pl
from tqdm import tqdm
import typer

from spider_silkome_module.config import INTERIM_DATA_DIR, PROCESSED_DATA_DIR
from spider_silkome_module.utils.run_cmd import run_cmd

app = typer.Typer()

DEFAULT_TASK_NAME = "protein_confirmation"

# Spidroin family groups (used for cross-checking type calls).
FAMILY_GROUPS: dict[str, str] = {
    "MaSp": "MaSp_grp", "MaSp1": "MaSp_grp", "MaSp2": "MaSp_grp",
    "MaSp2B": "MaSp_grp", "MaSp3": "MaSp_grp", "MaSp3B": "MaSp_grp",
    "AgSp1": "AgSp_grp", "AgSp2": "AgSp_grp",
    "Flag": "Flag_grp", "Pflag": "Flag_grp",
    "MiSp": "MiSp_grp", "AcSp": "AcSp_grp", "CySp": "CySp_grp",
    "CrSp": "CrSp_grp", "PySp": "PySp_grp", "TuSp": "TuSp_grp",
    "Spidroin": "Generic_grp", "hypo": "Generic_grp",
}

# Conserved CTD cysteine motif: most spidroin CTDs end with two cysteines
# spaced 30-50 aa apart in the C-terminal ~80 aa.
_CTD_CYS_RE = re.compile(r"C.{20,80}C")
# A simple signal-peptide proxy: ≥4 hydrophobic residues in positions 5-25.
_HYDRO_AA = set("AILMFWVCY")


# ── DIAMOND ──────────────────────────────────────────────────────────────────


def run_diamond_blastp(
    query_fa: Path,
    diamond_db: Path,
    out_tsv: Path,
    threads: int = 8,
    evalue: float = 1e-5,
    force: bool = False,
) -> None:
    """blastp -> tabular: qid sid pident length evalue bitscore qlen slen stitle."""
    cmd = (
        f"diamond blastp --quiet -p {threads} -d {diamond_db} -q {query_fa} "
        f"-o {out_tsv} -e {evalue} -k 1 --outfmt 6 qseqid sseqid pident length "
        f"evalue bitscore qlen slen stitle"
    )
    run_cmd(cmd, [out_tsv], force=force)


def parse_diamond_tsv(path: Path) -> pl.DataFrame:
    """Parse diamond blastp tabular output (custom outfmt above)."""
    if not path.exists() or path.stat().st_size == 0:
        return pl.DataFrame(schema={
            "spidroin_id": pl.Utf8, "diamond_best_hit_silkome": pl.Utf8,
            "diamond_pident": pl.Float64, "diamond_length": pl.Int64,
            "diamond_evalue": pl.Float64, "diamond_bitscore": pl.Float64,
            "diamond_best_hit_type": pl.Utf8,
        })
    df = pl.read_csv(
        path, separator="\t", has_header=False,
        new_columns=["spidroin_id", "diamond_best_hit_silkome", "diamond_pident",
                     "diamond_length", "diamond_evalue", "diamond_bitscore",
                     "qlen", "slen", "stitle"],
        schema_overrides={"diamond_pident": pl.Float64,
                          "diamond_evalue": pl.Float64,
                          "diamond_bitscore": pl.Float64},
    )
    # silkome IDs look like '1242_AcSp_CTD_610aa'  or  '1_MiSp_CTD'.
    df = df.with_columns(
        pl.col("diamond_best_hit_silkome").str.split("_").list.get(1, null_on_oob=True)
        .alias("diamond_best_hit_type")
    ).drop(["qlen", "slen", "stitle"])
    return df


# ── hmmscan ──────────────────────────────────────────────────────────────────


def run_hmmscan(
    query_fa: Path,
    hmm_db: Path,
    out_tbl: Path,
    threads: int = 8,
    evalue: float = 1e-5,
    force: bool = False,
) -> None:
    cmd = (
        f"hmmscan --noali --cpu {threads} -E {evalue} "
        f"--domtblout {out_tbl} {hmm_db} {query_fa} > /dev/null"
    )
    run_cmd(cmd, [out_tbl], force=force)


def parse_hmmscan_domtbl(path: Path) -> pl.DataFrame:
    """
    Parse hmmscan --domtblout into a tidy DataFrame.
    Field order (HMMER 3.x domtblout):
      0 target_name  3 query_name  6 full_evalue  7 full_score
      11 i_evalue   12 dom_score  17 ali_from   18 ali_to
    """
    rows: list[dict] = []
    if not path.exists():
        return pl.DataFrame(rows)
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        if len(f) < 23:
            continue
        rows.append({
            "spidroin_id": f[3],
            "domain_name": f[0],
            "full_evalue": float(f[6]),
            "i_evalue": float(f[12]),
            "ali_from": int(f[17]),
            "ali_to": int(f[18]),
            "dom_score": float(f[13]),
        })
    if not rows:
        return pl.DataFrame(schema={
            "spidroin_id": pl.Utf8, "domain_name": pl.Utf8,
            "full_evalue": pl.Float64, "i_evalue": pl.Float64,
            "ali_from": pl.Int64, "ali_to": pl.Int64,
            "dom_score": pl.Float64,
        })
    df = pl.DataFrame(rows)
    df = df.with_columns([
        pl.col("domain_name").str.split("_").list.get(0).alias("hmm_type"),
        pl.col("domain_name").str.split("_").list.get(1).alias("hmm_kind"),
    ])
    return df


def best_hmm_per_spidroin(df: pl.DataFrame) -> pl.DataFrame:
    """Pick the lowest-evalue NTD and CTD hit per spidroin_id."""
    if df.is_empty():
        return pl.DataFrame(schema={
            "spidroin_id": pl.Utf8,
            "hmm_ntd_best": pl.Utf8, "hmm_ntd_evalue": pl.Float64,
            "hmm_ntd_from": pl.Int64, "hmm_ntd_to": pl.Int64, "hmm_ntd_type": pl.Utf8,
            "hmm_ctd_best": pl.Utf8, "hmm_ctd_evalue": pl.Float64,
            "hmm_ctd_from": pl.Int64, "hmm_ctd_to": pl.Int64, "hmm_ctd_type": pl.Utf8,
            "hmm_domain_order_correct": pl.Boolean,
        })
    ntd = (df.filter(pl.col("hmm_kind") == "NTD")
             .sort("i_evalue")
             .group_by("spidroin_id")
             .first()
             .rename({
                 "domain_name": "hmm_ntd_best",
                 "i_evalue": "hmm_ntd_evalue",
                 "ali_from": "hmm_ntd_from",
                 "ali_to": "hmm_ntd_to",
                 "hmm_type": "hmm_ntd_type",
             })
             .select([
                 "spidroin_id", "hmm_ntd_best", "hmm_ntd_evalue",
                 "hmm_ntd_from", "hmm_ntd_to", "hmm_ntd_type",
             ]))
    ctd = (df.filter(pl.col("hmm_kind") == "CTD")
             .sort("i_evalue")
             .group_by("spidroin_id")
             .first()
             .rename({
                 "domain_name": "hmm_ctd_best",
                 "i_evalue": "hmm_ctd_evalue",
                 "ali_from": "hmm_ctd_from",
                 "ali_to": "hmm_ctd_to",
                 "hmm_type": "hmm_ctd_type",
             })
             .select([
                 "spidroin_id", "hmm_ctd_best", "hmm_ctd_evalue",
                 "hmm_ctd_from", "hmm_ctd_to", "hmm_ctd_type",
             ]))
    merged = ntd.join(ctd, on="spidroin_id", how="full", coalesce=True)
    merged = merged.with_columns(
        (pl.col("hmm_ntd_to").fill_null(0) < pl.col("hmm_ctd_from").fill_null(10**9))
        .alias("hmm_domain_order_correct")
    )
    return merged


# ── motif & length checks ────────────────────────────────────────────────────


def has_signal_peptide_motif(seq: str) -> bool:
    """Heuristic: ≥4 hydrophobic residues among AA positions 5–25 (1-based)."""
    if len(seq) < 25:
        return False
    window = seq[4:25]
    return sum(1 for aa in window if aa in _HYDRO_AA) >= 4


def has_ctd_cysteines(seq: str) -> tuple[bool, str]:
    """
    Look for two cysteines in the C-terminal 120 aa, spaced 20-80 aa apart.
    Returns (found, "pos1,pos2" or "").
    """
    if len(seq) < 60:
        return False, ""
    tail = seq[-min(150, len(seq)):]
    m = _CTD_CYS_RE.search(tail)
    if not m:
        return False, ""
    offset = len(seq) - len(tail)
    p1 = offset + m.start() + 1
    p2 = offset + m.end()  # 1-based position of the second Cys
    return True, f"{p1},{p2}"


def motif_table(proteins: dict[str, str]) -> pl.DataFrame:
    rows = []
    for sid, seq in proteins.items():
        sp_has = has_signal_peptide_motif(seq)
        cys_ok, cys_pos = has_ctd_cysteines(seq)
        rows.append({
            "spidroin_id": sid,
            "has_signal_peptide_motif": sp_has,
            "has_ctd_cysteines": cys_ok,
            "ctd_cysteine_positions": cys_pos,
        })
    return pl.DataFrame(rows) if rows else pl.DataFrame(schema={
        "spidroin_id": pl.Utf8,
        "has_signal_peptide_motif": pl.Boolean,
        "has_ctd_cysteines": pl.Boolean,
        "ctd_cysteine_positions": pl.Utf8,
    })


def length_zscore(df: pl.DataFrame) -> pl.DataFrame:
    """
    Compute length z-score within each spidroin_type using non-failed
    (predicted_protein_length > 0) rows.  Types with < 5 samples → NA.
    """
    if df.is_empty() or "spidroin_type" not in df.columns:
        return df.with_columns([
            pl.lit(None, dtype=pl.Float64).alias("length_zscore_within_type"),
            pl.lit(None, dtype=pl.Boolean).alias("length_outlier"),
        ])
    stats = (df.filter(pl.col("predicted_protein_length") > 0)
               .group_by("spidroin_type")
               .agg([
                   pl.col("predicted_protein_length").mean().alias("_mean"),
                   pl.col("predicted_protein_length").std().alias("_std"),
                   pl.len().alias("_n"),
               ]))
    df = df.join(stats, on="spidroin_type", how="left")
    df = df.with_columns(
        pl.when((pl.col("_n") >= 5) & (pl.col("_std") > 0) & (pl.col("predicted_protein_length") > 0))
        .then((pl.col("predicted_protein_length") - pl.col("_mean")) / pl.col("_std"))
        .otherwise(None)
        .alias("length_zscore_within_type")
    )
    df = df.with_columns(
        pl.when(pl.col("length_zscore_within_type").is_not_null())
        .then(pl.col("length_zscore_within_type").abs() > 3)
        .otherwise(None)
        .alias("length_outlier")
    ).drop(["_mean", "_std", "_n"])
    return df


# ── scoring ──────────────────────────────────────────────────────────────────


def compute_validation_score(df: pl.DataFrame) -> pl.DataFrame:
    """
    Combine all evidence dimensions into a 0-100 validation_score.
        A pairwise           30 pts (identity ≥0.8 cov ≥0.8 -> full 30)
        B diamond            30 pts (e<1e-20 + type consistent -> 30)
        C hmm                25 pts (NTD+CTD match type, e<1e-10 each)
        D motif              5  pts
        E length             5  pts
        F bonus              5  pts (premature stop -10 penalty)
    """
    expected = {"pairwise_identity", "pairwise_coverage_query",
                "diamond_evalue", "diamond_best_hit_type",
                "hmm_ntd_evalue", "hmm_ctd_evalue",
                "hmm_ntd_type", "hmm_ctd_type", "hmm_domain_order_correct",
                "spidroin_type", "has_signal_peptide_motif",
                "has_ctd_cysteines", "length_outlier",
                "predicted_protein_length", "has_premature_stop"}
    missing = expected - set(df.columns)
    for col in missing:
        df = df.with_columns(pl.lit(None).alias(col))

    df = df.with_columns([
        # A: pairwise vs NCBI
        (pl.min_horizontal(
            pl.col("pairwise_identity").fill_null(0) / 0.8,
            pl.lit(1.0)
         ) * 15
         + pl.min_horizontal(
            pl.col("pairwise_coverage_query").fill_null(0) / 0.8,
            pl.lit(1.0)
         ) * 15
        ).alias("_score_A"),
        # B: DIAMOND vs silkome (e-value contribution + type match)
        (pl.when(pl.col("diamond_evalue").fill_null(1e3) <= 1e-20).then(15)
            .when(pl.col("diamond_evalue").fill_null(1e3) <= 1e-10).then(10)
            .when(pl.col("diamond_evalue").fill_null(1e3) <= 1e-5).then(5)
            .otherwise(0)
         + pl.when(
             (pl.col("diamond_best_hit_type").fill_null("") == pl.col("spidroin_type").fill_null(""))
             & (pl.col("diamond_best_hit_type").fill_null("") != "")
         ).then(15).otherwise(0)
        ).alias("_score_B"),
        # C: HMM evidence — NTD and CTD found and matching type
        (pl.when(pl.col("hmm_ntd_evalue").fill_null(1.0) <= 1e-10).then(8).otherwise(0)
         + pl.when(pl.col("hmm_ctd_evalue").fill_null(1.0) <= 1e-10).then(8).otherwise(0)
         + pl.when(pl.col("hmm_ntd_type").fill_null("") == pl.col("spidroin_type").fill_null("")).then(4).otherwise(0)
         + pl.when(pl.col("hmm_ctd_type").fill_null("") == pl.col("spidroin_type").fill_null("")).then(4).otherwise(0)
         + pl.when(pl.col("hmm_domain_order_correct").fill_null(False)).then(1).otherwise(0)
        ).alias("_score_C"),
        # D: motif checks
        (pl.when(pl.col("has_signal_peptide_motif").fill_null(False)).then(2).otherwise(0)
         + pl.when(pl.col("has_ctd_cysteines").fill_null(False)).then(3).otherwise(0)
        ).alias("_score_D"),
        # E: length sanity
        (pl.when(pl.col("predicted_protein_length").fill_null(0) >= 100).then(3).otherwise(0)
         + pl.when(pl.col("length_outlier").fill_null(False).not_()).then(2).otherwise(0)
        ).alias("_score_E"),
        # bonus / penalty
        (pl.when(pl.col("has_premature_stop").fill_null(False)).then(-10).otherwise(5))
            .alias("_score_F"),
    ])

    df = df.with_columns(
        (pl.col("_score_A") + pl.col("_score_B") + pl.col("_score_C")
         + pl.col("_score_D") + pl.col("_score_E") + pl.col("_score_F"))
        .clip(0, 100)
        .round(1)
        .alias("validation_score")
    ).drop(["_score_A", "_score_B", "_score_C", "_score_D", "_score_E", "_score_F"])

    df = df.with_columns(
        pl.when(pl.col("validation_score") >= 80).then(pl.lit("validated"))
          .when(pl.col("validation_score") >= 50).then(pl.lit("partial"))
          .otherwise(pl.lit("failed"))
          .alias("validation_status")
    )
    return df


# ── per-species orchestration ────────────────────────────────────────────────


def load_typing_meta(typing_dir: Path, species_name: str) -> pl.DataFrame:
    """Load species-specific spidroin metadata (id, type)."""
    tsv = typing_dir / species_name / f"{species_name}.tsv"
    if not tsv.exists():
        return pl.DataFrame(schema={
            "spidroin_id": pl.Utf8, "spidroin_type": pl.Utf8})
    df = pl.read_csv(tsv, separator="\t", infer_schema_length=1000)
    keep = [c for c in df.columns if c.lower() in ("spidroin_id",
                                                    "spidroin_type",
                                                    "hint_type",
                                                    "full_length",
                                                    "stran",
                                                    "strand")]
    df = df.select(keep)
    rename = {c: c.lower() for c in df.columns if c.lower() != c}
    if rename:
        df = df.rename(rename)
    if "stran" in df.columns:
        df = df.rename({"stran": "strand"})
    df = df.with_columns(pl.col("spidroin_id").cast(pl.Utf8))
    return df


def add_family_group(df: pl.DataFrame) -> pl.DataFrame:
    return df.with_columns(
        pl.col("spidroin_type").map_elements(
            lambda t: FAMILY_GROUPS.get(t, "Other_grp"),
            return_dtype=pl.Utf8,
        ).alias("family_group")
    )


def process_species(
    species_name: str,
    confirmation_species_dir: Path,
    typing_dir: Path,
    diamond_db: Path,
    hmm_db: Path,
    threads: int,
    force: bool,
) -> pl.DataFrame:
    """Run B/C/D/E for one species and write per-dimension TSVs.
    Returns the combined per-spidroin DataFrame for cross-species concat.
    """
    pred_fa = confirmation_species_dir / "predicted_proteins.fa"
    pairwise_tsv = confirmation_species_dir / "pairwise_vs_ncbi.tsv"
    if not pred_fa.exists() or pred_fa.stat().st_size == 0:
        logger.warning(f"[{species_name}] no predicted proteins, skipping")
        return pl.DataFrame()

    diamond_tsv = confirmation_species_dir / "diamond_vs_silkome.tsv"
    domtbl = confirmation_species_dir / "hmmscan_domains.tsv"

    # Run external tools.
    run_diamond_blastp(pred_fa, diamond_db, diamond_tsv, threads=threads, force=force)
    run_hmmscan(pred_fa, hmm_db, domtbl, threads=threads, force=force)

    # Parse external tool outputs.
    diamond_df = parse_diamond_tsv(diamond_tsv)
    hmm_df = parse_hmmscan_domtbl(domtbl)
    hmm_best = best_hmm_per_spidroin(hmm_df)

    # Sequence-derived motif table.
    proteins = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(pred_fa, "fasta")}
    motif_df = motif_table(proteins)
    motif_df.write_csv(confirmation_species_dir / "motif_check.tsv", separator="\t")

    # Pairwise table from upstream step (dimension A).
    if pairwise_tsv.exists():
        pairwise_df = pl.read_csv(pairwise_tsv, separator="\t", infer_schema_length=1000)
    else:
        pairwise_df = pl.DataFrame(schema={"spidroin_id": pl.Utf8})

    typing_df = load_typing_meta(typing_dir, species_name)

    # Merge everything.
    merged = (
        pairwise_df
        .join(typing_df, on="spidroin_id", how="left")
        .join(diamond_df, on="spidroin_id", how="left")
        .join(hmm_best, on="spidroin_id", how="left")
        .join(motif_df, on="spidroin_id", how="left")
    )
    merged = add_family_group(merged)
    merged = length_zscore(merged)
    merged = compute_validation_score(merged)

    merged.write_csv(confirmation_species_dir / "protein_confirmation.tsv", separator="\t")
    return merged


# ── CLI ──────────────────────────────────────────────────────────────────────


@app.command()
def main(
    confirmation_dir: Path = PROCESSED_DATA_DIR / DEFAULT_TASK_NAME / "confirmation",
    typing_dir: Path = PROCESSED_DATA_DIR / "typing_results",
    diamond_db: Path = INTERIM_DATA_DIR / DEFAULT_TASK_NAME / "diamond_db" / "silkome.dmnd",
    hmm_db: Path = INTERIM_DATA_DIR / DEFAULT_TASK_NAME / "hmm_library" / "all_domains.hmm",
    threads: int = 8,
    species: str | None = None,
    force: bool = False,
):
    """
    Validate predicted proteins for one or all species.
    """
    species_dirs = sorted(d for d in confirmation_dir.iterdir() if d.is_dir())
    if species is not None:
        species_dirs = [d for d in species_dirs if d.name == species]
    if not species_dirs:
        logger.error(f"No species dirs under {confirmation_dir} matching filter")
        raise typer.Exit(1)

    logger.info(f"Validating {len(species_dirs)} species")
    all_rows: list[pl.DataFrame] = []
    for sp_dir in tqdm(species_dirs, desc="Validating proteins"):
        df = process_species(
            sp_dir.name, sp_dir, typing_dir,
            diamond_db, hmm_db, threads, force,
        )
        if not df.is_empty():
            all_rows.append(df)

    logger.success(f"Validation complete for {len(species_dirs)} species")


if __name__ == "__main__":
    app()
