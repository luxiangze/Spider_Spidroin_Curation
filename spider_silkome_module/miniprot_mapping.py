"""
Spidroin Curation Pipeline

This module provides the complete pipeline for spidroin gene curation:
1. Cluster spidroin sequences using MMseqs2
2. Align sequences to spider genomes using miniprot
3. Extract and combine gene predictions
"""

import os
from pathlib import Path
import re
import subprocess
import time
from typing import List, Optional

from loguru import logger
from natsort import natsort_keygen, ns
import pandas as pd

from spider_silkome_module import (
    EXTERNAL_DATA_DIR,
    INTERIM_DATA_DIR,
    PROCESSED_DATA_DIR,
    RAW_DATA_DIR,
    Attributes,
    GFFData,
    extract_positions_from_gff,
    positions_export,
)


def cluster_spidroin_sequences(
    spidroin_fasta_file: str,
    output_dir: str,
    min_seq_id: float = 0.9,
    coverage: float = 0.8,
) -> str:
    """
    Redundant clustering of Spidroin protein sequences using MMseqs2

    Parameters
    ----------
    spidroin_fasta_file : str
        Input spidroin FASTA file path
    output_dir : str
        Output directory
    min_seq_id : float
        Minimum sequence similarity threshold (default: 0.9)
    coverage : float
        Coverage threshold (default: 0.8)

    Returns
    -------
    str
        Representative sequences FASTA file path
    """
    os.makedirs(output_dir, exist_ok=True)

    base_name = Path(spidroin_fasta_file).stem
    output_prefix = f"{output_dir}/{base_name}"
    rep_seq_file = f"{output_prefix}_rep_seq.fasta"

    if os.path.exists(rep_seq_file):
        logger.info(f"Representative sequences already exist: {rep_seq_file}")
        return rep_seq_file

    logger.info(f"Running MMseqs2 clustering on {spidroin_fasta_file}")
    cmd = (
        f"pixi run --environment mmseqs mmseqs easy-cluster "
        f"{spidroin_fasta_file} {output_prefix} {output_dir}/tmp "
        f"--min-seq-id {min_seq_id} -c {coverage} --cov-mode 1"
    )

    subprocess.run(cmd, shell=True, check=True)
    logger.info(f"Clustering completed: {rep_seq_file}")

    return rep_seq_file


def index_genome(genome_file: str, output_dir: str, threads: int = 70) -> str:
    """
    Indexing the genome using miniprot

    Parameters
    ----------
    genome_file : str
        Input genome FASTA file path
    output_dir : str
        Output directory
    threads : int
        Number of threads (default: 70)

    Returns
    -------
    str
        Index file path
    """
    os.makedirs(output_dir, exist_ok=True)

    # Please make sure the genome_file end with 'fa.gz'
    if not genome_file.endswith(".fa.gz"):
        raise ValueError("Genome file must end with 'fa.gz'")

    genome_name = Path(genome_file).stem.replace(".fa.gz", "")
    index_file = f"{output_dir}/{genome_name}.mpi"

    if os.path.exists(index_file):
        logger.info(f"Index already exists: {index_file}")
        return index_file

    logger.info(f"Indexing genome: {genome_name}")
    start_time = time.time()

    cmd = f"miniprot -t{threads} -d {index_file} {genome_file}"
    subprocess.run(cmd, shell=True, check=True)

    elapsed = time.time() - start_time
    logger.info(f"Indexing {genome_name} completed in {elapsed:.2f} seconds")

    return index_file


def align_to_genome(
    index_file: str,
    query_fasta: str,
    output_dir: str,
    threads: int = 70,
) -> str:
    """
    Aligning sequences to genome using miniprot

    Parameters
    ----------
    index_file : str
        Genome index file path
    query_fasta : str
        Query sequence FASTA file path
    output_dir : str
        Output directory
    threads : int
        Number of threads (default: 70)

    Returns
    -------
    str
        GFF format alignment result file path
    """
    os.makedirs(output_dir, exist_ok=True)

    genome_name = Path(index_file).stem.replace(".mpi", "")
    gff_file = f"{output_dir}/{genome_name}.gff"

    if os.path.exists(gff_file):
        logger.info(f"Alignment already exists: {gff_file}")
        return gff_file

    logger.info(f"Aligning sequences to {genome_name}")
    start_time = time.time()

    cmd = f"miniprot -t {threads} -I --gff {index_file} {query_fasta} > {gff_file}"
    subprocess.run(cmd, shell=True, check=True)

    elapsed = time.time() - start_time
    logger.info(f"Alignment to {genome_name} completed in {elapsed:.2f} seconds")

    return gff_file


def extract_mrna_gff(gff_file: str) -> str:
    """
    Extracting mRNA records from GFF file and sorting

    Parameters
    ----------
    gff_file : str
        Input GFF file path

    Returns
    -------
    str
        mRNA GFF file path
    """
    mrna_gff_file = gff_file.replace(".gff", ".mRNA.gff")

    if os.path.exists(mrna_gff_file):
        logger.info(f"mRNA GFF already exists: {mrna_gff_file}")
        return mrna_gff_file

    logger.info(f"Extracting mRNA records from {gff_file}")
    cmd = f"grep 'mRNA' {gff_file} | sort -k1,1V -k4,4n > {mrna_gff_file}"
    subprocess.run(cmd, shell=True, check=True)

    return mrna_gff_file


def split_by_spidroin_type(mrna_gff_file: str) -> List[str]:
    """
    Splitting GFF file by spidroin type

    Parameters
    ----------
    mrna_gff_file : str
        mRNA GFF file path

    Returns
    -------
    List[str]
        List of GFF files split by spidroin type
    """
    output_dir = Path(mrna_gff_file).parent
    genome_name = Path(mrna_gff_file).stem.replace(".mRNA", "")

    # Reading GFF file and extracting spidroin types
    mRNA_gff = pd.read_csv(mrna_gff_file, sep="\t", header=None)
    gff_header = [
        "seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attribute",
    ]
    mRNA_gff.columns = gff_header

    spidroins = list(
        set([row["attribute"].split(";")[-1].split("|")[-2] for _, row in mRNA_gff.iterrows()])
    )

    logger.info(f"Found {len(spidroins)} spidroin types in {genome_name}: {spidroins}")

    output_files = []
    for spidroin in spidroins:
        gff_spidroin_output = f"{output_dir}/{genome_name}.mRNA.{spidroin}.gff"

        if os.path.exists(gff_spidroin_output):
            logger.debug(f"{genome_name}.mRNA.{spidroin}.gff already exists")
        else:
            grep_cmd = f"grep '|{spidroin}|' {mrna_gff_file} > {gff_spidroin_output}"
            subprocess.run(grep_cmd, shell=True, check=True)

        output_files.append(gff_spidroin_output)

    return output_files


def parse_gff_file(gff_file: str) -> List[GFFData]:
    """
    Parsing GFF file into GFFData objects

    Parameters
    ----------
    gff_file : str
        GFF file path

    Returns
    -------
    List[GFFData]
        List of GFFData objects
    """
    gff_data = []

    with open(gff_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")

            # Parsing attributes field
            attr_dict = {}
            for attr in fields[8].split(";"):
                if "=" in attr:
                    key, value = attr.split("=", 1)
                    attr_dict[key] = value

            # Creating attributes object
            attr_obj = Attributes(
                ID=attr_dict["ID"],
                Rank=int(attr_dict["Rank"]),
                Identity=float(attr_dict["Identity"]),
                Positive=float(attr_dict["Positive"]),
                Target=attr_dict["Target"].split("|"),
            )

            # Creating GFFData object
            gff_data.append(
                GFFData(
                    seqid=fields[0],
                    source=fields[1],
                    type=fields[2],
                    start=int(fields[3]),
                    end=int(fields[4]),
                    score=float(fields[5]),
                    strand=fields[6],
                    frame=fields[7],
                    attributes=attr_obj,
                )
            )

    return gff_data


def process_spidroin_type(
    spidroin_gff_file: str,
    positive_threshold: float = 0.75,
    min_length: int = 1000,
    max_length: int = 100000,
    extension_length: int = 10000,
) -> tuple:
    """
    Processing single spidroin type GFF file

    Parameters
    ----------
    spidroin_gff_file : str
        spidroin type GFF file path
    positive_threshold : float
        Positive threshold (default: 0.75)
    min_length : int
        Minimum gene length threshold (default: 1000)
    max_length : int
        Maximum gene length threshold (default: 100000)
    extension_length : int
        Extension length when start or end is missing (default: 10000)

    Returns
    -------
    tuple
        (DataFrame, List[dict]) - CSV data and GFF records
    """
    spidroin = Path(spidroin_gff_file).stem.split(".")[-1]

    # Parsing GFF files
    spidroin_gff_data = parse_gff_file(spidroin_gff_file)

    # Sort by Positive in descending order
    spidroin_gff_data.sort(key=lambda x: x.attributes.Positive, reverse=True)

    # Extract positions
    positions = extract_positions_from_gff(
        spidroin_gff_data, positive_threshold=positive_threshold
    )

    # Export CSV data
    csv_output = spidroin_gff_file.replace(".gff", ".csv")
    df = positions_export(positions, csv_output, format="csv")

    # Export GFF data
    gff_output = spidroin_gff_file.replace(".gff", ".combined.gff")
    gff_records = positions_export(
        positions,
        gff_output,
        format="gff",
        spidroin=spidroin,
        min_length=min_length,
        max_length=max_length,
        extension_length=extension_length,
    )

    return df, gff_records


def combine_all_spidroins(
    spider_name: str,
    spidroin_gff_files: List[str],
    output_dir: str,
    positive_threshold: float = 0.75,
    min_length: int = 1000,
    max_length: int = 100000,
    extension_length: int = 10000,
) -> None:
    """
    Combine all spidroin type predictions

    Parameters
    ----------
    spider_name : str
        Spider species name
    spidroin_gff_files : List[str]
        List of spidroin type GFF file paths
    output_dir : str
        Output directory
    positive_threshold : float
        Positive threshold (default: 0.75)
    min_length : int
        Minimum gene length threshold (default: 1000)
    max_length : int
        Maximum gene length threshold (default: 100000)
    extension_length : int
        Extension length when start or end is missing (default: 10000)
    """
    os.makedirs(output_dir, exist_ok=True)

    all_gff_records = []
    all_csv = pd.DataFrame()

    logger.info(f"Processing {len(spidroin_gff_files)} spidroin types for {spider_name}")

    for spidroin_gff_file in spidroin_gff_files:
        df, gff_records = process_spidroin_type(
            spidroin_gff_file,
            positive_threshold=positive_threshold,
            min_length=min_length,
            max_length=max_length,
            extension_length=extension_length,
        )

        all_csv = pd.concat([all_csv, df], ignore_index=True)
        all_gff_records.extend(gff_records)

    # Save combined GFF file
    if all_gff_records:
        df_combined = pd.DataFrame(all_gff_records)
        # Sort by chromosome and start position with natural sorting
        nkey = natsort_keygen(alg=ns.IGNORECASE)

        def get_sort_key(seqid):
            """Generate sort key with priority: chr > ctg > others"""
            if re.match(r"^chr", seqid, flags=re.IGNORECASE):
                priority = 0
            elif re.match(r"^ctg", seqid, flags=re.IGNORECASE):
                priority = 1
            else:
                priority = 2
            return (priority, nkey(seqid))

        df_combined["sort_key"] = df_combined["seqid"].apply(lambda x: get_sort_key(x))
        df_combined = df_combined.sort_values(["sort_key", "start"]).drop("sort_key", axis=1)

        # Write combined GFF file
        gff_output = f"{output_dir}/{spider_name}.gff"
        with open(gff_output, "w") as f:
            f.write("##gff-version 3\n")
            df_combined.to_csv(f, sep="\t", header=False, index=False)

        logger.info(f"Combined GFF saved to {gff_output} ({len(all_gff_records)} genes)")
    else:
        logger.warning(f"No GFF records to combine for {spider_name}")

    # Save combined CSV file
    if not all_csv.empty:
        csv_output = f"{output_dir}/{spider_name}.csv"
        all_csv.to_csv(csv_output, index=False)
        logger.info(f"Combined CSV saved to {csv_output}")


def run_miniprot_mapping_pipeline(
    spider_genome_path: Optional[str] = None,
    spidroin_fasta_file: Optional[str] = None,
    positive_threshold: float = 0.75,
    min_length: int = 1000,
    max_length: int = 100000,
    extension_length: int = 1000,
) -> None:
    """
    Run the complete miniprot mapping pipeline

    Parameters
    ----------
    spider_genome_path : str, optional
        Spider genome file directory (default: data/raw/spider_genome)
    spidroin_fasta_file : str, optional
        Spidroin protein sequence FASTA file (default: data/external/spider-silkome-database.v1.prot.fixed.fasta)
    positive_threshold : float
        Positive threshold (default: 0.75)
    min_length : int
        Minimum gene length threshold (default: 1000)
    max_length : int
        Maximum gene length threshold (default: 100000)
    extension_length : int
        Extension length when start or end is missing (default: 10000)
    """
    # Set default paths
    if spider_genome_path is None:
        spider_genome_path = f"{RAW_DATA_DIR}/spider_genome"

    if spidroin_fasta_file is None:
        spidroin_fasta_file = f"{EXTERNAL_DATA_DIR}/spider-silkome-database.v1.prot.fixed.fasta"

    # Step 1: Cluster spidroin sequences with MMseqs2
    logger.info("=" * 80)
    logger.info("Step 1: Clustering spidroin sequences with MMseqs2")
    logger.info("=" * 80)

    mmseqs_output_dir = f"{INTERIM_DATA_DIR}/mmseqs"
    spidroin_fasta_file_rep = cluster_spidroin_sequences(
        spidroin_fasta_file,
        mmseqs_output_dir,
        min_seq_id=0.9,
        coverage=0.8,
    )

    # Manually curated representative sequence file
    spidroin_fasta_file_rep_manually = spidroin_fasta_file_rep.replace(".fasta", "_manually.fasta")

    if not os.path.exists(spidroin_fasta_file_rep_manually):
        logger.warning(
            f"Manually curated file not found: {spidroin_fasta_file_rep_manually}\n"
            f"If you want to obtain more accurate prediction results, you can manually edit {spidroin_fasta_file_rep}, then rename the file to {spidroin_fasta_file_rep_manually}, and rerun the script."
            f"Using automated clustering result: {spidroin_fasta_file_rep}"
        )
        spidroin_fasta_file_rep_manually = spidroin_fasta_file_rep

    # Step 2: Process each spider genome
    logger.info("=" * 80)
    logger.info("Step 2: Processing spider genomes")
    logger.info("=" * 80)

    spider_genomes = [f for f in os.listdir(spider_genome_path) if f.endswith(".fa.gz")]
    logger.info(f"Found {len(spider_genomes)} spider genomes")

    genome_mpi_path = f"{INTERIM_DATA_DIR}/genome_mpi"
    miniprot_output_path = f"{INTERIM_DATA_DIR}/miniprot"
    final_output_dir = f"{PROCESSED_DATA_DIR}/01.miniprot_mapping"

    for spider_genome in spider_genomes:
        spider_name = spider_genome.split(".")[0]
        logger.info(f"\n{'=' * 80}")
        logger.info(f"Processing: {spider_name}")
        logger.info(f"{'=' * 80}")

        # Step 2.1: Index genome
        genome_file = f"{spider_genome_path}/{spider_genome}"
        index_file = index_genome(genome_file, genome_mpi_path)

        # Step 2.2: Align to genome
        output_dir = f"{miniprot_output_path}/{spider_name}_all"
        gff_file = align_to_genome(index_file, spidroin_fasta_file_rep_manually, output_dir)

        # Step 2.3: Extract mRNA records
        mrna_gff_file = extract_mrna_gff(gff_file)

        # Step 2.4: Split by spidroin type
        spidroin_gff_files = split_by_spidroin_type(mrna_gff_file)

        # Step 2.5: Combine all spidroin type predictions
        combine_all_spidroins(
            spider_name,
            spidroin_gff_files,
            final_output_dir,
            positive_threshold=positive_threshold,
            min_length=min_length,
            max_length=max_length,
            extension_length=extension_length,
        )

    logger.info("=" * 80)
    logger.info("Spidroin curation pipeline completed!")
    logger.info(f"Final results saved to: {final_output_dir}")
    logger.info("=" * 80)


def main():
    """Main function entry point"""
    import argparse

    parser = argparse.ArgumentParser(
        description="Spidroin miniprot mapping pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--genome-path", type=str, default=None, help="Path to spider genome directory"
    )

    parser.add_argument(
        "--spidroin-fasta", type=str, default=None, help="Path to spidroin protein FASTA file"
    )

    parser.add_argument(
        "--positive-threshold",
        type=float,
        default=0.75,
        help="Positive score threshold for filtering alignments",
    )

    parser.add_argument(
        "--min-length", type=int, default=1000, help="Minimum gene length threshold"
    )

    parser.add_argument(
        "--max-length", type=int, default=100000, help="Maximum gene length threshold"
    )

    parser.add_argument(
        "--extension-length",
        type=int,
        default=1000,
        help="Extension length when start or end is missing",
    )

    args = parser.parse_args()

    run_miniprot_mapping_pipeline(
        spider_genome_path=args.genome_path,
        spidroin_fasta_file=args.spidroin_fasta,
        positive_threshold=args.positive_threshold,
        min_length=args.min_length,
        max_length=args.max_length,
        extension_length=args.extension_length,
    )


if __name__ == "__main__":
    main()
