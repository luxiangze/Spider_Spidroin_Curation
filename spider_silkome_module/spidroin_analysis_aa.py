#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Spidroin Analysis Pipeline (Amino Acid Version)

Description:
    Analyze spider silk protein data from amino acid HMMER search results (hmmsearch).
    Parse domtblout files to identify full-length spidroins and partial matches.

Input:
    - HMMER search result directories (containing .domtbl files)
    - Species subdirectories with NTD/CTD domain search results

Output:
    - Per-target classification table (TSV)
    - Per-gene summary table (TSV)
    - Statistics summary

Author: Generated for spider_silkome project
Date: 2024
"""

import argparse
import glob
import os
import re
import subprocess
from dataclasses import dataclass, field
from enum import Enum
from typing import Optional

import matplotlib.pyplot as plt
import polars as pl
import seaborn as sns
import typer
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from loguru import logger
from tqdm import tqdm

app = typer.Typer()


class SpidClassification(str, Enum):
    """Classification of spidroin based on NTD/CTD domain hits."""
    FULL_LENGTH = "full_length"
    NTD_ONLY = "ntd_only"
    CTD_ONLY = "ctd_only"
    BOTH_HITS_NOT_FULL = "both_hits_not_full"
    NONE = "none"


@dataclass
class DomainHit:
    """Represents a single domain hit from domtblout file."""
    target_id: str
    tlen: int
    query_id: str
    qlen: int
    full_evalue: float
    full_score: float
    full_bias: float
    dom_index: int
    dom_count: int
    c_evalue: float
    i_evalue: float
    dom_score: float
    dom_bias: float
    hmm_from: int
    hmm_to: int
    ali_from: int
    ali_to: int
    env_from: int
    env_to: int
    acc: float
    species: str = ""
    source_file: str = ""
    domain_type: str = ""  # NTD or CTD
    spidroin_type: str = ""  # e.g., MaSp, AcSp, etc.

    @property
    def env_len(self) -> int:
        """Length of envelope coordinates."""
        return abs(self.env_to - self.env_from) + 1

    @property
    def env_coverage(self) -> float:
        """Envelope coverage as fraction of qlen."""
        return self.env_len / self.qlen if self.qlen > 0 else 0.0

    @property
    def is_first_dom(self) -> bool:
        """Check if this is the first domain."""
        return self.dom_index == 1

    @property
    def is_last_dom(self) -> bool:
        """Check if this is the last domain."""
        return self.dom_index == self.dom_count

    @property
    def ntd_end_ok(self) -> bool:
        """Check if NTD endpoint condition is satisfied (env_from == 1)."""
        return self.env_from == 1

    @property
    def ctd_end_ok(self) -> bool:
        """Check if CTD endpoint condition is satisfied (env_to == tlen)."""
        return self.env_to == self.tlen


@dataclass
class TargetSummary:
    """Summary of domain hits for a single target sequence."""
    target_id: str
    tlen: int
    species: str
    ntd_hits: list = field(default_factory=list)
    ctd_hits: list = field(default_factory=list)
    best_ntd: Optional[DomainHit] = None
    best_ctd: Optional[DomainHit] = None
    ntd_end_ok: bool = False
    ctd_end_ok: bool = False
    classification: SpidClassification = SpidClassification.NONE
    ntd_profile: str = ""
    ctd_profile: str = ""
    ntd_spidroin_type: str = ""
    ctd_spidroin_type: str = ""


def parse_domtblout(
    filepath: str,
    species: str,
    c_evalue_threshold: float = 1e-10,
    require_score_gt_bias: bool = True,
    min_env_coverage: Optional[float] = None,
) -> list[DomainHit]:
    """
    Parse a domtblout file and filter hits.

    Args:
        filepath: Path to domtblout file
        species: Species name
        c_evalue_threshold: Maximum c-Evalue threshold
        require_score_gt_bias: Require domain score > domain bias
        min_env_coverage: Minimum env coverage as fraction of qlen (optional)

    Returns:
        List of filtered DomainHit objects
    """
    hits = []
    filename = os.path.basename(filepath)

    # Determine domain type and spidroin type from filename
    # Pattern: {spidroin_type}_{domain_type}.domtbl
    match = re.match(r"(.+)_(NTD|CTD)\.domtbl", filename)
    if not match:
        logger.warning(f"Cannot parse domain type from filename: {filename}")
        return hits

    spidroin_type = match.group(1)
    domain_type = match.group(2)

    try:
        with open(filepath, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue

                parts = line.strip().split()
                if len(parts) < 22:
                    continue

                try:
                    hit = DomainHit(
                        target_id=parts[0],
                        tlen=int(parts[2]),
                        query_id=parts[3],
                        qlen=int(parts[5]),
                        full_evalue=float(parts[6]),
                        full_score=float(parts[7]),
                        full_bias=float(parts[8]),
                        dom_index=int(parts[9]),
                        dom_count=int(parts[10]),
                        c_evalue=float(parts[11]),
                        i_evalue=float(parts[12]),
                        dom_score=float(parts[13]),
                        dom_bias=float(parts[14]),
                        hmm_from=int(parts[15]),
                        hmm_to=int(parts[16]),
                        ali_from=int(parts[17]),
                        ali_to=int(parts[18]),
                        env_from=int(parts[19]),
                        env_to=int(parts[20]),
                        acc=float(parts[21]),
                        species=species,
                        source_file=filepath,
                        domain_type=domain_type,
                        spidroin_type=spidroin_type,
                    )

                    # Apply filters
                    if hit.c_evalue >= c_evalue_threshold:
                        continue
                    if require_score_gt_bias and hit.dom_score <= hit.dom_bias:
                        continue
                    if min_env_coverage is not None and hit.env_coverage < min_env_coverage:
                        continue

                    hits.append(hit)

                except (ValueError, IndexError) as e:
                    logger.debug(f"Failed to parse line in {filepath}: {e}")
                    continue

    except Exception as e:
        logger.error(f"Error reading {filepath}: {e}")

    return hits


def collect_all_hits(
    input_dir: str,
    c_evalue_threshold: float = 1e-10,
    require_score_gt_bias: bool = True,
    min_env_coverage: Optional[float] = None,
) -> list[DomainHit]:
    """
    Collect all domain hits from input directory.

    Args:
        input_dir: Directory containing species subdirectories with domtbl files
        c_evalue_threshold: Maximum c-Evalue threshold
        require_score_gt_bias: Require domain score > domain bias
        min_env_coverage: Minimum env coverage as fraction of qlen (optional)

    Returns:
        List of all filtered DomainHit objects
    """
    all_hits = []

    # Find all species directories
    species_dirs = [d for d in glob.glob(f"{input_dir}/*") if os.path.isdir(d)]

    if not species_dirs:
        # Maybe input_dir itself contains domtbl files
        domtbl_files = glob.glob(f"{input_dir}/*.domtbl")
        if domtbl_files:
            species = os.path.basename(input_dir)
            for filepath in tqdm(domtbl_files, desc=f"Parsing {species}"):
                hits = parse_domtblout(
                    filepath, species, c_evalue_threshold,
                    require_score_gt_bias, min_env_coverage
                )
                all_hits.extend(hits)
        return all_hits

    for species_dir in tqdm(species_dirs, desc="Processing species"):
        species = os.path.basename(species_dir)
        domtbl_files = glob.glob(f"{species_dir}/*.domtbl")

        for filepath in domtbl_files:
            hits = parse_domtblout(
                filepath, species, c_evalue_threshold,
                require_score_gt_bias, min_env_coverage
            )
            all_hits.extend(hits)

    return all_hits


def select_best_terminal_hit(
    hits: list[DomainHit],
    domain_type: str,
) -> Optional[DomainHit]:
    """
    Select the best terminal domain hit.

    For NTD: select from hits where dom_index == 1, prefer smaller c_evalue
    For CTD: select from hits where dom_index == dom_count, prefer smaller c_evalue

    Args:
        hits: List of hits for one domain type
        domain_type: 'NTD' or 'CTD'

    Returns:
        Best hit or None
    """
    if not hits:
        return None

    if domain_type == "NTD":
        # Filter to first domains only
        terminal_hits = [h for h in hits if h.is_first_dom]
    else:  # CTD
        # Filter to last domains only
        terminal_hits = [h for h in hits if h.is_last_dom]

    if not terminal_hits:
        return None

    # Sort by c_evalue (lower is better), then by dom_score (higher is better)
    terminal_hits.sort(key=lambda h: (h.c_evalue, -h.dom_score))
    return terminal_hits[0]


def aggregate_target_hits(
    all_hits: list[DomainHit],
) -> dict[tuple[str, str], TargetSummary]:
    """
    Aggregate hits by target and species, select best NTD/CTD.

    Args:
        all_hits: All filtered domain hits

    Returns:
        Dictionary mapping (species, target_id) to TargetSummary
    """
    # Group hits by (species, target_id)
    target_hits: dict[tuple[str, str], dict] = {}

    for hit in all_hits:
        key = (hit.species, hit.target_id)
        if key not in target_hits:
            target_hits[key] = {
                "target_id": hit.target_id,
                "tlen": hit.tlen,
                "species": hit.species,
                "ntd_hits": [],
                "ctd_hits": [],
            }

        if hit.domain_type == "NTD":
            target_hits[key]["ntd_hits"].append(hit)
        elif hit.domain_type == "CTD":
            target_hits[key]["ctd_hits"].append(hit)

    # Create summaries
    summaries = {}
    for key, data in target_hits.items():
        summary = TargetSummary(
            target_id=data["target_id"],
            tlen=data["tlen"],
            species=data["species"],
            ntd_hits=data["ntd_hits"],
            ctd_hits=data["ctd_hits"],
        )

        # Select best NTD (from first domains across all profiles)
        summary.best_ntd = select_best_terminal_hit(data["ntd_hits"], "NTD")
        # Select best CTD (from last domains across all profiles)
        summary.best_ctd = select_best_terminal_hit(data["ctd_hits"], "CTD")

        # Check endpoint conditions
        if summary.best_ntd:
            summary.ntd_end_ok = summary.best_ntd.ntd_end_ok
            summary.ntd_profile = summary.best_ntd.query_id
            summary.ntd_spidroin_type = summary.best_ntd.spidroin_type

        if summary.best_ctd:
            summary.ctd_end_ok = summary.best_ctd.ctd_end_ok
            summary.ctd_profile = summary.best_ctd.query_id
            summary.ctd_spidroin_type = summary.best_ctd.spidroin_type

        # Classify
        summary.classification = classify_target(summary)

        summaries[key] = summary

    return summaries


def classify_target(summary: TargetSummary) -> SpidClassification:
    """
    Classify a target based on NTD/CTD hits and endpoint conditions.

    Args:
        summary: TargetSummary object

    Returns:
        SpidClassification enum value
    """
    has_ntd = summary.best_ntd is not None
    has_ctd = summary.best_ctd is not None

    if has_ntd and has_ctd:
        if summary.ntd_end_ok and summary.ctd_end_ok:
            return SpidClassification.FULL_LENGTH
        else:
            return SpidClassification.BOTH_HITS_NOT_FULL
    elif has_ntd and summary.ntd_end_ok:
        return SpidClassification.NTD_ONLY
    elif has_ctd and summary.ctd_end_ok:
        return SpidClassification.CTD_ONLY
    elif has_ntd or has_ctd:
        # Has hits but endpoint conditions not met
        return SpidClassification.BOTH_HITS_NOT_FULL
    else:
        return SpidClassification.NONE


def assign_spidroin_type(summary: TargetSummary) -> str:
    """
    Assign final spidroin type based on NTD and CTD types.

    Args:
        summary: TargetSummary object

    Returns:
        Assigned spidroin type string
    """
    ntd_type = summary.ntd_spidroin_type
    ctd_type = summary.ctd_spidroin_type

    if not ntd_type and not ctd_type:
        return "Unknown"
    if not ntd_type:
        return ctd_type
    if not ctd_type:
        return ntd_type

    # Both present
    if ntd_type == ctd_type:
        return ntd_type

    # Handle special cases
    maj_ampullates = ["MaSp", "MaSp1", "MaSp2", "MaSp2B", "MaSp3", "MaSp3B"]

    if ntd_type == "Other" or ctd_type == "Other":
        return ctd_type if ntd_type == "Other" else ntd_type

    if ntd_type == "MaSp" and ctd_type in maj_ampullates:
        return ctd_type
    if ctd_type == "MaSp" and ntd_type in maj_ampullates:
        return ntd_type

    # Discordant
    return f"{ntd_type}/{ctd_type}"


def summaries_to_dataframe(
    summaries: dict[tuple[str, str], TargetSummary],
) -> pl.DataFrame:
    """
    Convert target summaries to a polars DataFrame.

    Args:
        summaries: Dictionary of TargetSummary objects

    Returns:
        DataFrame with one row per target
    """
    rows = []
    for (species, target_id), summary in summaries.items():
        spidroin_type = assign_spidroin_type(summary)

        row = {
            "species": species,
            "target_id": target_id,
            "tlen": summary.tlen,
            "classification": summary.classification.value,
            "spidroin_type": spidroin_type,
            "ntd_profile": summary.ntd_profile,
            "ctd_profile": summary.ctd_profile,
            "ntd_spidroin_type": summary.ntd_spidroin_type,
            "ctd_spidroin_type": summary.ctd_spidroin_type,
            "ntd_end_ok": summary.ntd_end_ok,
            "ctd_end_ok": summary.ctd_end_ok,
            "ntd_c_evalue": summary.best_ntd.c_evalue if summary.best_ntd else None,
            "ntd_dom_score": summary.best_ntd.dom_score if summary.best_ntd else None,
            "ntd_env_from": summary.best_ntd.env_from if summary.best_ntd else None,
            "ntd_env_to": summary.best_ntd.env_to if summary.best_ntd else None,
            "ntd_env_coverage": summary.best_ntd.env_coverage if summary.best_ntd else None,
            "ctd_c_evalue": summary.best_ctd.c_evalue if summary.best_ctd else None,
            "ctd_dom_score": summary.best_ctd.dom_score if summary.best_ctd else None,
            "ctd_env_from": summary.best_ctd.env_from if summary.best_ctd else None,
            "ctd_env_to": summary.best_ctd.env_to if summary.best_ctd else None,
            "ctd_env_coverage": summary.best_ctd.env_coverage if summary.best_ctd else None,
            "num_ntd_hits": len(summary.ntd_hits),
            "num_ctd_hits": len(summary.ctd_hits),
        }
        rows.append(row)

    if not rows:
        return pl.DataFrame()

    df = pl.DataFrame(rows)
    df = df.sort(["species", "target_id"])
    return df


def generate_statistics(df: pl.DataFrame) -> pl.DataFrame:
    """
    Generate classification statistics by species.

    Args:
        df: Target results DataFrame

    Returns:
        Statistics DataFrame
    """
    if df.height == 0:
        return pl.DataFrame()

    # Group by species and classification, count occurrences
    stats = (
        df.group_by(["species", "classification"])
        .len()
        .pivot(on="classification", index="species", values="len")
        .fill_null(0)
    )

    # Add total column
    classification_cols = [c for c in stats.columns if c != "species"]

    # Cast all numeric columns to Int64 for consistency
    stats = stats.with_columns([
        pl.col(c).cast(pl.Int64) for c in classification_cols
    ])

    stats = stats.with_columns(
        pl.sum_horizontal(classification_cols).alias("total")
    )

    # Add totals row
    totals = {"species": "TOTAL"}
    for col in classification_cols + ["total"]:
        totals[col] = int(stats[col].sum())
    totals_df = pl.DataFrame([totals])
    stats = pl.concat([stats, totals_df])

    return stats


def find_protein_fasta(genome_dir: str, species: str) -> Optional[str]:
    """
    Find protein FASTA file for a species.

    Args:
        genome_dir: Base directory containing species subdirectories
        species: Species name

    Returns:
        Path to protein FASTA file or None if not found
    """
    patterns = [
        f"{genome_dir}/{species}/*.protein.fasta",
        f"{genome_dir}/{species}/*.protein.fa",
        f"{genome_dir}/{species}/*_protein.fasta",
        f"{genome_dir}/{species}/*.pep.fa",
        f"{genome_dir}/{species}/*.pep.fasta",
    ]

    for pattern in patterns:
        files = glob.glob(pattern)
        if files:
            return files[0]

    return None


def find_gff_file(genome_dir: str, species: str) -> Optional[str]:
    """
    Find GFF/GTF annotation file for a species.

    Args:
        genome_dir: Base directory containing species subdirectories
        species: Species name

    Returns:
        Path to GFF file or None if not found
    """
    patterns = [
        f"{genome_dir}/{species}/*.gff",
        f"{genome_dir}/{species}/*.gff3",
        f"{genome_dir}/{species}/*.gtf",
    ]

    for pattern in patterns:
        files = glob.glob(pattern)
        if files:
            # Prefer non-fixed files first
            for f in files:
                if "fixed" not in f:
                    return f
            return files[0]

    return None


def parse_gff_coordinates(gff_path: str) -> dict[str, dict]:
    """
    Parse GFF file to get transcript/mRNA coordinates.

    Args:
        gff_path: Path to GFF file

    Returns:
        Dictionary mapping transcript ID to {scaffold, strand, start, end}
    """
    coords = {}

    with open(gff_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            feature_type = parts[2]
            if feature_type not in ("mRNA", "transcript"):
                continue

            scaffold = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attributes = parts[8]

            # Parse ID from attributes
            transcript_id = None
            for attr in attributes.split(";"):
                if attr.startswith("ID="):
                    transcript_id = attr[3:]
                    break

            if transcript_id:
                coords[transcript_id] = {
                    "scaffold": scaffold,
                    "strand": strand,
                    "start": start,
                    "end": end,
                }

    return coords


def load_protein_sequences(fasta_path: str) -> dict[str, SeqRecord]:
    """
    Load protein sequences from FASTA file into a dictionary.

    Args:
        fasta_path: Path to protein FASTA file

    Returns:
        Dictionary mapping sequence ID to SeqRecord
    """
    sequences = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        # Use first part of ID (before any space)
        seq_id = record.id.split()[0]
        sequences[seq_id] = record
    return sequences


def extract_full_length_sequences(
    df: pl.DataFrame,
    genome_dir: str,
    output_path: str,
) -> int:
    """
    Extract protein sequences for full-length spidroins and write to FASTA.

    Args:
        df: DataFrame with full-length spidroin classifications
        genome_dir: Directory containing species protein FASTA files
        output_path: Output FASTA file path

    Returns:
        Number of sequences extracted
    """
    if df.height == 0:
        logger.warning("No full-length spidroins to extract")
        return 0

    # Group by species
    species_list = df["species"].unique().to_list()

    # Track spidroin type counts per species for naming
    type_counters: dict[str, dict[str, int]] = {}

    extracted_records = []
    missing_sequences = []

    for species in tqdm(species_list, desc="Extracting sequences"):
        # Find protein FASTA for this species
        fasta_path = find_protein_fasta(genome_dir, species)
        if not fasta_path:
            logger.warning(f"No protein FASTA found for {species}")
            continue

        # Load sequences
        try:
            sequences = load_protein_sequences(fasta_path)
        except Exception as e:
            logger.error(f"Failed to load {fasta_path}: {e}")
            continue

        # Find and parse GFF for coordinates
        gff_path = find_gff_file(genome_dir, species)
        coords = {}
        if gff_path:
            try:
                coords = parse_gff_coordinates(gff_path)
            except Exception as e:
                logger.warning(f"Failed to parse GFF {gff_path}: {e}")

        # Get targets for this species
        species_df = df.filter(pl.col("species") == species)

        # Initialize counter for this species
        if species not in type_counters:
            type_counters[species] = {}

        for row in species_df.iter_rows(named=True):
            target_id = row["target_id"]
            spidroin_type = row["spidroin_type"]

            # Get sequence
            if target_id not in sequences:
                missing_sequences.append(f"{species}/{target_id}")
                continue

            seq_record = sequences[target_id]

            # Update counter
            if spidroin_type not in type_counters[species]:
                type_counters[species][spidroin_type] = 0
            type_counters[species][spidroin_type] += 1
            count = type_counters[species][spidroin_type]

            # Create new ID: {species}_{spidroin_type}_{count}
            new_id = f"{species}_{spidroin_type}_{count}"

            # Create description: {scaffold}_{strand}_{start}_{end}
            if target_id in coords:
                c = coords[target_id]
                description = f"{c['scaffold']}_{c['strand']}_{c['start']}_{c['end']}"
            else:
                description = target_id

            # Create new record
            new_record = SeqRecord(
                seq_record.seq,
                id=new_id,
                description=description,
            )
            extracted_records.append(new_record)

    if missing_sequences:
        logger.warning(f"Missing {len(missing_sequences)} sequences: {missing_sequences[:5]}...")

    # Write output
    if extracted_records:
        SeqIO.write(extracted_records, output_path, "fasta")
        logger.info(f"Extracted {len(extracted_records)} sequences to {output_path}")

    return len(extracted_records)


def find_cds_fasta(genome_dir: str, species: str) -> Optional[str]:
    """
    Find CDS FASTA file for a species.

    Args:
        genome_dir: Base directory containing species subdirectories
        species: Species name

    Returns:
        Path to CDS FASTA file or None if not found
    """
    patterns = [
        f"{genome_dir}/{species}/*.cds.fasta",
        f"{genome_dir}/{species}/*.cds.fa",
        f"{genome_dir}/{species}/*_cds.fasta",
        f"{genome_dir}/{species}/*.transcript.fasta",
        f"{genome_dir}/{species}/*.mrna.fasta",
    ]

    for pattern in patterns:
        files = glob.glob(pattern)
        if files:
            return files[0]

    return None


def find_genome_fasta(genome_dir: str, species: str) -> Optional[str]:
    """
    Find genome FASTA file for a species.

    Args:
        genome_dir: Base directory containing species subdirectories
        species: Species name

    Returns:
        Path to genome FASTA file or None if not found
    """
    patterns = [
        f"{genome_dir}/{species}/{species}.fa",
        f"{genome_dir}/{species}/{species}.fasta",
        f"{genome_dir}/{species}/*.genome.fa",
        f"{genome_dir}/{species}/*.fa",
    ]

    for pattern in patterns:
        files = glob.glob(pattern)
        if files:
            # Exclude protein and other non-genome files
            for f in files:
                if "protein" not in f and "pep" not in f and "cds" not in f:
                    return f

    return None


def extract_cds_with_gffread(
    genome_fasta: str,
    gff_path: str,
    output_cds: str,
    transcript_ids: list[str],
) -> bool:
    """
    Extract CDS sequences using gffread.

    Args:
        genome_fasta: Path to genome FASTA file
        gff_path: Path to GFF annotation file
        output_cds: Output CDS FASTA file path
        transcript_ids: List of transcript IDs to extract

    Returns:
        True if successful, False otherwise
    """
    try:
        # Create temp file with transcript IDs
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            for tid in transcript_ids:
                f.write(f"{tid}\n")
            id_file = f.name

        # Run gffread to extract CDS
        cmd = [
            "gffread", "-x", output_cds,
            "-g", genome_fasta,
            gff_path,
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)

        os.unlink(id_file)

        if result.returncode != 0:
            logger.warning(f"gffread failed: {result.stderr}")
            return False

        return True
    except Exception as e:
        logger.warning(f"Failed to run gffread: {e}")
        return False


def extract_sequences_to_fasta(
    df: pl.DataFrame,
    genome_dir: str,
    output_dir: str,
    classification_filter: str = "full_length",
) -> tuple[int, int]:
    """
    Extract protein and CDS sequences for spidroins and write to FASTA files.

    Args:
        df: DataFrame with spidroin classifications
        genome_dir: Directory containing species genome files
        output_dir: Output directory
        classification_filter: Classification to filter (default: full_length)

    Returns:
        Tuple of (protein_count, cds_count)
    """
    filtered_df = df.filter(pl.col("classification") == classification_filter)

    if filtered_df.height == 0:
        logger.warning(f"No {classification_filter} spidroins to extract")
        return 0, 0

    species_list = filtered_df["species"].unique().to_list()

    protein_records = []
    cds_records = []
    missing_proteins = []

    for species in tqdm(species_list, desc="Extracting sequences"):
        # Find protein FASTA
        fasta_path = find_protein_fasta(genome_dir, species)
        if not fasta_path:
            logger.warning(f"No protein FASTA found for {species}")
            continue

        # Load protein sequences
        try:
            proteins = load_protein_sequences(fasta_path)
        except Exception as e:
            logger.error(f"Failed to load {fasta_path}: {e}")
            continue

        # Try to find CDS FASTA
        cds_path = find_cds_fasta(genome_dir, species)
        cds_sequences = {}
        if cds_path:
            try:
                cds_sequences = load_protein_sequences(cds_path)
            except Exception as e:
                logger.warning(f"Failed to load CDS {cds_path}: {e}")

        # If no CDS file, try gffread extraction
        if not cds_sequences:
            gff_path = find_gff_file(genome_dir, species)
            genome_fasta = find_genome_fasta(genome_dir, species)
            if gff_path and genome_fasta:
                temp_cds = os.path.join(output_dir, f"{species}_temp_cds.fa")
                species_df = filtered_df.filter(pl.col("species") == species)
                transcript_ids = species_df["target_id"].to_list()
                if extract_cds_with_gffread(genome_fasta, gff_path, temp_cds, transcript_ids):
                    try:
                        cds_sequences = load_protein_sequences(temp_cds)
                    except Exception:
                        pass
                    finally:
                        if os.path.exists(temp_cds):
                            os.unlink(temp_cds)

        # Get targets for this species
        species_df = filtered_df.filter(pl.col("species") == species)

        for row in species_df.iter_rows(named=True):
            target_id = row["target_id"]
            spidroin_type = row["spidroin_type"]

            # Get protein sequence
            if target_id not in proteins:
                missing_proteins.append(f"{species}/{target_id}")
                continue

            seq_record = proteins[target_id]

            # Create new record with gene ID as ID, spidroin type as description
            # Format: >{gene_id} {spidroin_type}
            new_protein = SeqRecord(
                seq_record.seq,
                id=target_id,
                description=spidroin_type,
            )
            protein_records.append(new_protein)

            # Get CDS sequence if available
            if target_id in cds_sequences:
                cds_record = cds_sequences[target_id]
                new_cds = SeqRecord(
                    cds_record.seq,
                    id=target_id,
                    description=spidroin_type,
                )
                cds_records.append(new_cds)

    if missing_proteins:
        logger.warning(f"Missing {len(missing_proteins)} protein sequences")

    # Write output files
    protein_count = 0
    cds_count = 0

    if protein_records:
        protein_file = os.path.join(output_dir, f"{classification_filter}_spidroins.protein.fasta")
        SeqIO.write(protein_records, protein_file, "fasta")
        protein_count = len(protein_records)
        logger.info(f"Extracted {protein_count} protein sequences to {protein_file}")

    if cds_records:
        cds_file = os.path.join(output_dir, f"{classification_filter}_spidroins.cds.fasta")
        SeqIO.write(cds_records, cds_file, "fasta")
        cds_count = len(cds_records)
        logger.info(f"Extracted {cds_count} CDS sequences to {cds_file}")

    return protein_count, cds_count


def plot_spidroin_heatmap(
    df: pl.DataFrame,
    output_prefix: str,
    spidroin_order: list[str] = None,
) -> None:
    """
    Generate heatmap of spidroin counts by species and type.

    Args:
        df: DataFrame with spidroin classifications
        output_prefix: Output file prefix
        spidroin_order: Optional list to specify spidroin type order
    """
    if df.height == 0:
        logger.warning("No data for heatmap")
        return

    # Create pivot table
    pivot_data = (
        df.group_by(["species", "spidroin_type"])
        .len()
        .pivot(on="spidroin_type", index="species", values="len")
        .fill_null(0)
    )

    # Convert to pandas for seaborn
    pivot_pd = pivot_data.to_pandas()
    pivot_pd = pivot_pd.set_index("species")

    # Reorder columns if specified
    if spidroin_order:
        existing_cols = [c for c in spidroin_order if c in pivot_pd.columns]
        other_cols = [c for c in pivot_pd.columns if c not in spidroin_order]
        pivot_pd = pivot_pd[existing_cols + other_cols]

    # Save pivot table
    pivot_pd.to_csv(f"{output_prefix}_pivot_table.tsv", sep="\t")

    # Create heatmap
    sns.set_theme()
    plt.figure(figsize=(12, 8))
    sns.heatmap(
        pivot_pd, cmap='YlGnBu', linewidths=0.5,
        annot=True, fmt='d', cbar=True
    )
    plt.tight_layout()
    plt.xlabel('Spidroin Type')
    plt.ylabel('Species')
    plt.title('Spidroin Counts by Species and Type')

    # Save as PDF and PNG
    plt.savefig(f"{output_prefix}_heatmap.pdf", bbox_inches="tight")
    plt.savefig(f"{output_prefix}_heatmap.png", bbox_inches="tight", dpi=300)
    plt.close()

    logger.info(f"Saved heatmap to {output_prefix}_heatmap.pdf/png")


def plot_classification_heatmap(
    df: pl.DataFrame,
    output_prefix: str,
) -> None:
    """
    Generate heatmap of classification counts by species.

    Args:
        df: DataFrame with spidroin classifications
        output_prefix: Output file prefix
    """
    if df.height == 0:
        logger.warning("No data for classification heatmap")
        return

    # Create pivot table
    pivot_data = (
        df.group_by(["species", "classification"])
        .len()
        .pivot(on="classification", index="species", values="len")
        .fill_null(0)
    )

    # Convert to pandas for seaborn
    pivot_pd = pivot_data.to_pandas()
    pivot_pd = pivot_pd.set_index("species")

    # Order columns
    col_order = ["full_length", "ntd_only", "ctd_only", "both_hits_not_full", "none"]
    existing_cols = [c for c in col_order if c in pivot_pd.columns]
    pivot_pd = pivot_pd[existing_cols]

    # Create heatmap
    sns.set_theme()
    plt.figure(figsize=(10, 8))
    sns.heatmap(
        pivot_pd, cmap='RdYlGn_r', linewidths=0.5,
        annot=True, fmt='d', cbar=True
    )
    plt.tight_layout()
    plt.xlabel('Classification')
    plt.ylabel('Species')
    plt.title('Spidroin Classification by Species')

    # Save as PDF and PNG
    plt.savefig(f"{output_prefix}_classification_heatmap.pdf", bbox_inches="tight")
    plt.savefig(f"{output_prefix}_classification_heatmap.png", bbox_inches="tight", dpi=300)
    plt.close()

    logger.info(f"Saved classification heatmap to {output_prefix}_classification_heatmap.pdf/png")


def plot_spidroin_barplot(
    df: pl.DataFrame,
    output_prefix: str,
) -> None:
    """
    Generate stacked bar plot of spidroin types by species.

    Args:
        df: DataFrame with spidroin classifications
        output_prefix: Output file prefix
    """
    if df.height == 0:
        logger.warning("No data for barplot")
        return

    # Create pivot table
    pivot_data = (
        df.group_by(["species", "spidroin_type"])
        .len()
        .pivot(on="spidroin_type", index="species", values="len")
        .fill_null(0)
    )

    # Convert to pandas
    pivot_pd = pivot_data.to_pandas()
    pivot_pd = pivot_pd.set_index("species")

    # Create stacked bar plot
    sns.set_theme()
    pivot_pd.plot(kind="bar", stacked=True, figsize=(12, 6), colormap="tab20")
    plt.xlabel("Species")
    plt.ylabel("Count")
    plt.title("Spidroin Types by Species")
    plt.legend(title="Spidroin Type", bbox_to_anchor=(1.02, 1), loc='upper left')
    plt.tight_layout()

    # Save as PDF and PNG
    plt.savefig(f"{output_prefix}_barplot.pdf", bbox_inches="tight")
    plt.savefig(f"{output_prefix}_barplot.png", bbox_inches="tight", dpi=300)
    plt.close()

    logger.info(f"Saved barplot to {output_prefix}_barplot.pdf/png")


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Spidroin Analysis Pipeline - Analyze spider silk protein data from amino acid HMMER search results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic usage
    python spidroin_analysis_aa.py -i hmmer_search_output -o results

    # With env coverage filter
    python spidroin_analysis_aa.py -i hmmer_search_output -o results --min-env-coverage 0.9

    # Custom E-value threshold
    python spidroin_analysis_aa.py -i hmmer_search_output -o results -e 1e-5
        """,
    )

    # Required arguments
    parser.add_argument(
        "-i",
        "--input-dir",
        required=True,
        help="Directory containing HMMER search results (includes .domtbl files)",
    )

    # Optional arguments
    parser.add_argument(
        "-o",
        "--output-dir",
        default="spidroin_analysis_aa",
        help="Output directory (default: spidroin_analysis_aa)",
    )
    parser.add_argument(
        "-e",
        "--c-evalue",
        type=float,
        default=1e-10,
        help="c-Evalue threshold for filtering (default: 1e-10)",
    )
    parser.add_argument(
        "--min-env-coverage",
        type=float,
        default=None,
        help="Minimum env coverage as fraction of qlen (optional, e.g., 0.9)",
    )
    parser.add_argument(
        "--no-score-bias-filter",
        action="store_true",
        help="Disable the score > bias filter",
    )

    return parser.parse_args()


@app.command()
def main(
    input_dir: str = typer.Option(..., "-i", "--input-dir", help="Input directory with domtbl files"),
    output_dir: str = typer.Option("spidroin_analysis_aa", "-o", "--output-dir", help="Output directory"),
    genome_dir: Optional[str] = typer.Option(None, "-g", "--genome-dir", help="Genome directory with protein FASTA files"),
    c_evalue: float = typer.Option(1e-10, "-e", "--c-evalue", help="c-Evalue threshold"),
    min_env_coverage: Optional[float] = typer.Option(None, "--min-env-coverage", help="Min env coverage (0-1)"),
    no_score_bias_filter: bool = typer.Option(False, "--no-score-bias-filter", help="Disable score>bias filter"),
):
    """Analyze spider silk protein data from amino acid HMMER search results."""

    # Setup output directory
    os.makedirs(output_dir, exist_ok=True)

    logger.info(f"Input directory: {input_dir}")
    logger.info(f"Output directory: {output_dir}")
    if genome_dir:
        logger.info(f"Genome directory: {genome_dir}")
    logger.info(f"c-Evalue threshold: {c_evalue}")
    if min_env_coverage:
        logger.info(f"Min env coverage: {min_env_coverage}")

    # Collect and filter all hits
    logger.info("Parsing domtblout files...")
    all_hits = collect_all_hits(
        input_dir,
        c_evalue_threshold=c_evalue,
        require_score_gt_bias=not no_score_bias_filter,
        min_env_coverage=min_env_coverage,
    )
    logger.info(f"Found {len(all_hits)} filtered domain hits")

    if len(all_hits) == 0:
        logger.warning("No hits found after filtering. Check input files and thresholds.")
        return

    # Aggregate by target
    logger.info("Aggregating hits by target...")
    summaries = aggregate_target_hits(all_hits)
    logger.info(f"Found {len(summaries)} unique targets with hits")

    # Convert to DataFrame
    df = summaries_to_dataframe(summaries)

    # Save results
    output_file = os.path.join(output_dir, "target_classification.tsv")
    df.write_csv(output_file, separator="\t")
    logger.info(f"Saved target classification to {output_file}")

    # Generate and save statistics
    stats = generate_statistics(df)
    stats_file = os.path.join(output_dir, "classification_stats.tsv")
    stats.write_csv(stats_file, separator="\t")
    logger.info(f"Saved statistics to {stats_file}")

    # Print summary
    logger.info("\n=== Classification Summary ===")
    for cls in SpidClassification:
        count = df.filter(pl.col("classification") == cls.value).height
        logger.info(f"  {cls.value}: {count}")

    # Save full-length and partial lists
    full_length = df.filter(pl.col("classification") == SpidClassification.FULL_LENGTH.value)
    if full_length.height > 0:
        full_file = os.path.join(output_dir, "full_length_spidroins.tsv")
        full_length.write_csv(full_file, separator="\t")
        logger.info(f"Saved {full_length.height} full-length spidroins to {full_file}")

    partial = df.filter(
        pl.col("classification").is_in([
            SpidClassification.NTD_ONLY.value,
            SpidClassification.CTD_ONLY.value,
        ])
    )
    if partial.height > 0:
        partial_file = os.path.join(output_dir, "partial_spidroins.tsv")
        partial.write_csv(partial_file, separator="\t")
        logger.info(f"Saved {partial.height} partial spidroins to {partial_file}")

    # Extract protein and CDS sequences for full-length spidroins
    if genome_dir and full_length.height > 0:
        logger.info("Extracting sequences...")
        extract_sequences_to_fasta(df, genome_dir, output_dir, "full_length")

    # Generate plots
    logger.info("Generating plots...")
    plot_prefix = os.path.join(output_dir, "spidroin")

    # Spidroin type heatmap
    spidroin_order = [
        "MaSp1", "MaSp2", "MaSp3", "MaSp", "MiSp", "Flag", "AcSp",
        "TuSp", "CySp", "TuSp_CySp", "PySp", "AgSp", "CrSp", "Other", "Unknown"
    ]
    plot_spidroin_heatmap(df, plot_prefix, spidroin_order)

    # Classification heatmap
    plot_classification_heatmap(df, plot_prefix)

    # Spidroin barplot
    plot_spidroin_barplot(df, plot_prefix)

    logger.info("Analysis complete!")


def run_cli():
    """Entry point for command-line interface."""
    app()


if __name__ == "__main__":
    app()
