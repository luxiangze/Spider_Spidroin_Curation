from pathlib import Path

from loguru import logger
from tqdm import tqdm
import typer
import polars as pl

from spider_silkome_module.config import INTERIM_DATA_DIR

app = typer.Typer()

# Column names for nhmmer tblout format
TBL_COLUMNS = [
    "target_name",
    "query_name",
    "query_accession",
    "hmm_from",
    "hmm_to",
    "ali_from",
    "ali_to",
    "env_from",
    "env_to",
    "sq_len",
    "strand",
    "e_value",
    "score",
    "bias",
    "description",
]


def clean_query_name(query_name: str) -> str:
    """
    Remove _trimmed suffix from HMM model name.
    """
    if query_name.endswith("_trimmed"):
        return query_name[:-8]
    return query_name


def parse_tbl_file(tbl_path: Path, species_name: str = "") -> pl.DataFrame:
    """
    Parse a single nhmmer tblout file into a DataFrame.
    Removes _trimmed suffix from query_name.
    """
    rows = []
    with open(tbl_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 15:
                continue
            row = {
                "target_name": parts[0],
                "query_name": clean_query_name(parts[2]),
                "query_accession": parts[3],
                "hmm_from": int(parts[4]),
                "hmm_to": int(parts[5]),
                "ali_from": int(parts[6]),
                "ali_to": int(parts[7]),
                "env_from": int(parts[8]),
                "env_to": int(parts[9]),
                "sq_len": int(parts[10]),
                "strand": parts[11],
                "e_value": float(parts[12]),
                "score": float(parts[13]),
                "bias": float(parts[14]),
                "description": " ".join(parts[15:]) if len(parts) > 15 else "",
            }
            rows.append(row)

    if not rows:
        return pl.DataFrame(schema={col: pl.Utf8 for col in TBL_COLUMNS[:11]} |
                           {"e_value": pl.Float64, "score": pl.Float64, "bias": pl.Float64, "description": pl.Utf8})

    return pl.DataFrame(rows)


def get_hmm_length(df: pl.DataFrame) -> int:
    """
    Estimate HMM profile length from the max hmm_to value.
    """
    if df.is_empty():
        return 0
    return df["hmm_to"].max()


def filter_hits(
    df: pl.DataFrame,
    e_value_threshold: float = 1e-10,
    coverage_threshold: float = 0.9,
) -> pl.DataFrame:
    """
    Filter nhmmer hits based on:
    - E-value < threshold
    - HMM profile coverage >= threshold
    - score > bias
    """
    if df.is_empty():
        return df

    hmm_length = get_hmm_length(df)
    if hmm_length == 0:
        return df.head(0)

    filtered = df.filter(
        (pl.col("e_value") < e_value_threshold) &
        ((pl.col("hmm_to") - pl.col("hmm_from") + 1) / hmm_length >= coverage_threshold) &
        (pl.col("score") > pl.col("bias"))
    )

    return filtered


def to_gff(df: pl.DataFrame, species: str, source: str = "nhmmer") -> str:
    """
    Convert filtered hits to GFF3 format.
    Uses env_from and env_to for coordinates.
    """
    gff_lines = ["##gff-version 3"]

    for row in df.iter_rows(named=True):
        seqid = row["target_name"]
        start = min(row["env_from"], row["env_to"])
        end = max(row["env_from"], row["env_to"])
        score = row["score"]
        strand = row["strand"]
        query_name = row["query_name"]
        e_value = row["e_value"]

        # Extract feature type (CTD or NTD) from query_name
        if "_CTD" in query_name:
            feature_type = "CTD"
        elif "_NTD" in query_name:
            feature_type = "NTD"
        else:
            feature_type = "nhmmer_hit"

        attributes = f"ID={query_name}_{seqid}_{start}_{end};Name={query_name};E-value={e_value};Species={species}"

        gff_line = f"{seqid}\t{source}\t{feature_type}\t{start}\t{end}\t{score}\t{strand}\t.\t{attributes}"
        gff_lines.append(gff_line)

    return "\n".join(gff_lines)


def process_species(
    species_dir: Path,
    output_dir: Path,
    e_value_threshold: float,
    coverage_threshold: float,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """
    Process all tbl files for a species.
    Returns (all_hits_df, filtered_hits_df).
    """
    species_name = species_dir.name
    tbl_files = list(species_dir.glob("*.tbl"))

    all_hits = []
    filtered_hits = []

    for tbl_file in tbl_files:
        df = parse_tbl_file(tbl_file, species_name)
        if df.is_empty():
            continue

        df = df.with_columns(pl.lit(species_name).alias("species"))
        df = df.with_columns(pl.lit(tbl_file.stem).alias("hmm_profile"))
        all_hits.append(df)

        filtered_df = filter_hits(df, e_value_threshold, coverage_threshold)
        if not filtered_df.is_empty():
            filtered_hits.append(filtered_df)

    all_df = pl.concat(all_hits) if all_hits else pl.DataFrame()
    filtered_df = pl.concat(filtered_hits) if filtered_hits else pl.DataFrame()

    return all_df, filtered_df


@app.command()
def main(
    input_path: Path = INTERIM_DATA_DIR / "automated_spidroin_annotation" / "nhmmer_search",
    output_path: Path = INTERIM_DATA_DIR / "automated_spidroin_annotation" / "nhmmer_filtered",
    e_value_threshold: float = 1e-10,
    coverage_threshold: float = 0.9,
):
    """
    Parse and filter nhmmer search results.

    Filtering criteria:
    - E-value < threshold (default: 1e-10)
    - HMM profile coverage >= threshold (default: 90%)
    - score > bias

    Outputs:
    - all_hits.tsv: All hits summary
    - filtered_hits.tsv: Filtered hits
    - {species}/{species}.gff: GFF file for each species
    """
    logger.info(f"Input path: {input_path}")
    logger.info(f"Output path: {output_path}")
    logger.info(f"E-value threshold: {e_value_threshold}")
    logger.info(f"Coverage threshold: {coverage_threshold}")

    output_path.mkdir(parents=True, exist_ok=True)

    species_dirs = [d for d in input_path.iterdir() if d.is_dir()]
    logger.info(f"Found {len(species_dirs)} species directories")

    all_species_hits = []
    all_species_filtered = []

    for species_dir in tqdm(species_dirs, desc="Processing species"):
        species_name = species_dir.name
        logger.info(f"Processing {species_name}")

        all_df, filtered_df = process_species(
            species_dir, output_path, e_value_threshold, coverage_threshold
        )

        if not all_df.is_empty():
            all_species_hits.append(all_df)

        if not filtered_df.is_empty():
            all_species_filtered.append(filtered_df)

            species_output_dir = output_path / species_name
            species_output_dir.mkdir(parents=True, exist_ok=True)

            gff_content = to_gff(filtered_df, species_name)
            gff_file = species_output_dir / f"{species_name}.gff"
            gff_file.write_text(gff_content)
            logger.info(f"  Wrote {len(filtered_df)} filtered hits to {gff_file}")

    if all_species_hits:
        all_hits_df = pl.concat(all_species_hits)
        all_hits_file = output_path / "all_hits.tsv"
        all_hits_df.write_csv(all_hits_file, separator="\t")
        logger.info(f"Wrote {len(all_hits_df)} total hits to {all_hits_file}")

    if all_species_filtered:
        filtered_df = pl.concat(all_species_filtered)
        filtered_file = output_path / "filtered_hits.tsv"
        filtered_df.write_csv(filtered_file, separator="\t")
        logger.info(f"Wrote {len(filtered_df)} filtered hits to {filtered_file}")

    logger.success("Completed parsing nhmmer results")


if __name__ == "__main__":
    app()
