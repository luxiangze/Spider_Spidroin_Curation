#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Spidroin Analysis Pipeline

Description:
    This script analyzes spidroin (spider silk protein) data from HMMER search results.
    It performs the following tasks:
    1. Reads and parses HMMER table outputs from nhmmer searches
    2. Filters spidroin hits based on E-value thresholds and HMM profile coverage
    3. Joins N-terminal (NTD) and C-terminal (CTD) domains to identify full spidroin genes
    4. Generates various plots including heatmaps, bar plots, and length distributions
    5. Extracts spidroin sequences from genome assemblies

Input:
    - HMMER search result directories (containing .tbl files)
    - HMM profile directory (containing .hmm files)
    - Genome assembly directory (containing genomic FASTA files)

Output:
    - Summary tables (TSV format) of spidroin hits
    - Filtered spidroin tables
    - Joined domain tables
    - Various plots (PNG and SVG format)
    - Extracted spidroin sequences (FASTA format)

Dependencies:
    - pandas
    - seaborn
    - matplotlib
    - biopython

Author: Adapted from Schöneberg et al. (2025)
Date: 2024
"""

import argparse
import glob
import os
import pathlib
import re

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from matplotlib.ticker import AutoMinorLocator


def read_hmmer_tbl(infile: str, header_names: list, species: str) -> pd.DataFrame:
    """
    Read and parse HMMER table output file.

    Args:
        infile: Path to HMMER table file
        header_names: Column names for the output DataFrame
        species: Species name to add as a column

    Returns:
        DataFrame with parsed HMMER results
    """
    cols_to_use = range(0, 15)

    try:
        table = pd.read_csv(
            infile, sep=r'\s+', skiprows=2,
            skipfooter=10, header=None, usecols=cols_to_use, engine="python"
        )
    except pd.errors.ParserError:
        print(f"Warning: Pandas parser Error, No HMMER hits in output file {infile}?")
        table = pd.DataFrame(columns=cols_to_use)
    table[len(table.columns)] = species
    table.columns = header_names

    table["query_name"] = [
        "_".join(spid_type.split("_")[:2]) for spid_type in table["query_name"].to_list()
    ]
    table["spidroin_type"] = [
        spid_type.split("_")[0] for spid_type in table["query_name"].to_list()
    ]
    table["domain"] = [
        spid_type.split("_")[1] for spid_type in table["query_name"].to_list()
    ]
    return table


def get_hmm_length(hmm_profile_dir: str) -> dict:
    """
    Get the length of each HMM profile.

    Args:
        hmm_profile_dir: Directory containing HMM profile files

    Returns:
        Dictionary mapping profile names to their lengths
    """
    hmmer_profile_files = glob.glob(hmm_profile_dir + "/*.hmm")
    hmmer_length_dict = {}
    for file in hmmer_profile_files:
        with open(file) as hmmer_file:
            hmmer_length = int(hmmer_file.readlines()[2].strip().split("  ")[-1])
        spidroin_name = pathlib.PurePath(file).stem
        hmmer_length_dict[spidroin_name] = hmmer_length
    return hmmer_length_dict


def get_spidroin_intervals(spidroin_table: pd.DataFrame) -> pd.IntervalIndex:
    """
    Get alignment intervals from spidroin table.

    Args:
        spidroin_table: DataFrame with spidroin hits

    Returns:
        IntervalIndex of alignment positions (empty if no data)
    """
    if len(spidroin_table) == 0:
        return pd.IntervalIndex.from_arrays([], [], closed="neither")

    intervals = []
    for _, row in spidroin_table.iterrows():
        interval = row[["ali_from", "ali_to"]].values.tolist()
        interval = [min(interval), max(interval)]
        intervals.append(interval)

    intervals = pd.DataFrame(intervals)
    intervals = pd.IntervalIndex.from_arrays(intervals[0], intervals[1], closed="neither")
    return intervals


def get_overlapping_spidroins(
    spidroin_table: pd.DataFrame,
    spidroin_hit: pd.Series,
    intervals: pd.IntervalIndex,
    offset: int = 0
) -> pd.DataFrame:
    """
    Find overlapping spidroin hits.

    Args:
        spidroin_table: DataFrame with all spidroin hits
        spidroin_hit: Single spidroin hit to check
        intervals: IntervalIndex of all alignments
        offset: Extension to apply to the interval

    Returns:
        DataFrame with overlapping hits
    """
    boundaries = spidroin_hit[["ali_from", "ali_to"]].values.tolist()
    aln_range = pd.Interval(min(boundaries) - offset, max(boundaries) + offset)
    scaffold = spidroin_hit["target_name"]
    species = spidroin_hit["Species"]
    duplicates = spidroin_table.loc[
        (spidroin_table["Species"] == species) &
        (spidroin_table["target_name"] == scaffold) &
        (intervals.overlaps(aln_range))
    ]
    return duplicates


def filter_spidroins(
    spidroin_hits: pd.DataFrame,
    e_value_threshold: float,
    hmm_lengths: dict,
    len_factor: float = 0.9,
    hit_num: int = 0
) -> pd.DataFrame:
    """
    Filter spidroin hits by E-value and HMM profile coverage.

    Args:
        spidroin_hits: DataFrame with spidroin hits
        e_value_threshold: Maximum E-value threshold
        hmm_lengths: Dictionary of HMM profile lengths
        len_factor: Minimum fraction of profile to be covered
        hit_num: Which hit to select when duplicates exist (0 = best)

    Returns:
        Filtered DataFrame
    """
    if len(spidroin_hits) == 0:
        return spidroin_hits.copy()

    filtered_table = spidroin_hits.loc[spidroin_hits["E-value"] < e_value_threshold]
    filtered_table = filtered_table.loc[filtered_table["score"] > filtered_table["bias"]]

    if len(filtered_table) == 0:
        return filtered_table.copy()

    filtered_table["hmm_length"] = filtered_table["hmm_to"] - filtered_table["hmm_from"]

    spidroin_profile_lens = []
    for _, hit in filtered_table.iterrows():
        spidroin_profile_lens.append(hmm_lengths.get(hit["query_name"]))

    filtered_table["profile_len"] = spidroin_profile_lens
    filtered_table = filtered_table.loc[
        filtered_table["hmm_length"] > (filtered_table["profile_len"] * len_factor)
    ]

    if len(filtered_table) == 0:
        return filtered_table.copy()

    intervals = get_spidroin_intervals(filtered_table)

    deduplicated_table = []
    for _, hit in filtered_table.iterrows():
        duplicates = get_overlapping_spidroins(filtered_table, hit, intervals)
        duplicates.sort_values("E-value", inplace=True)
        deduplicated_table.append(duplicates.iloc[hit_num].tolist())

    deduplicated_table = pd.DataFrame(
        deduplicated_table, columns=filtered_table.columns.values.tolist()
    )
    deduplicated_table.drop_duplicates(inplace=True)

    return deduplicated_table


def plot_spidroins_heatmap(
    spidroin_table: pd.DataFrame,
    outprefix: str,
    column_order: list = None,
    all_species: list = None
) -> None:
    """
    Generate heatmap of spidroin counts by species and type.

    Args:
        spidroin_table: DataFrame with spidroin data
        outprefix: Output file prefix
        column_order: Optional list to specify column order (spidroin types)
        all_species: Optional list of all species to include (even if all zeros)
    """
    sns.set_theme()
    pivot_table = spidroin_table.pivot_table(
        index='Species', columns='spidroin_type', aggfunc='size', fill_value=0
    )

    # Add missing species with all zeros
    if all_species is not None:
        for species in all_species:
            if species not in pivot_table.index:
                pivot_table.loc[species] = 0
        pivot_table = pivot_table.loc[all_species]

    # Reorder columns if specified
    if column_order is not None:
        existing_cols = [c for c in column_order if c in pivot_table.columns]
        pivot_table = pivot_table[existing_cols]

    pivot_table.to_csv(f"{outprefix}_pivot_table.tsv", sep="\t")

    plt.figure(figsize=(10, 7))
    sns.heatmap(pivot_table, cmap='YlGnBu', linewidths=.5, annot=True, fmt='d', cbar=True)
    plt.tight_layout()
    plt.xlabel('Spidroin Type')
    plt.ylabel('Species')
    plt.title('Spidroin Terminal Hits by Species')
    plt.savefig(f"{outprefix}_heatmap.png", bbox_inches="tight")
    plt.savefig(f"{outprefix}_heatmap.svg", bbox_inches="tight")
    plt.clf()


def plot_domain_ratio(spidroin_table: pd.DataFrame, outprefix: str) -> None:
    """
    Plot domain ratio bar charts.

    Args:
        spidroin_table: DataFrame with spidroin data
        outprefix: Output file prefix
    """
    pivot_table = spidroin_table.pivot_table(
        index='Species', columns='domain', aggfunc='size', fill_value=0
    )

    pivot_table.plot(kind="bar", stacked=True, colormap="Paired")
    plt.xlabel("Species")
    plt.ylabel("Occurences")
    plt.title("Detected Terminal Domains")
    plt.savefig(f"{outprefix}_barplot.png", bbox_inches="tight")
    plt.savefig(f"{outprefix}_barplot.svg", bbox_inches="tight")
    plt.clf()

    pivot_table.loc["all"] = pivot_table.sum(axis=0, numeric_only=True)
    normalized_table = pivot_table.div(pivot_table.sum(axis=1), axis=0)
    normalized_table.plot(kind="bar", stacked=True, colormap="Paired")
    plt.xlabel("Species")
    plt.ylabel("Occurences")
    plt.title("Detected Terminal Domains")
    plt.savefig(f"{outprefix}_normalized_barplot.png", bbox_inches="tight")
    plt.clf()


def plot_domain_evals(spidroin_table: pd.DataFrame, outprefix: str) -> None:
    """
    Plot E-value distribution by domain.

    Args:
        spidroin_table: DataFrame with spidroin data
        outprefix: Output file prefix
    """
    sns.boxplot(x="domain", y="E-value", data=spidroin_table)
    plt.xlabel("Terminal Domain")
    plt.ylabel("E-value")
    plt.yscale("log")
    plt.title("E-values of Terminal Domains")
    plt.savefig(f"{outprefix}_e-values_boxplot.png")
    plt.savefig(f"{outprefix}_e-values_boxplot.svg")
    plt.clf()


def plot_second_best_classified(
    spidroin_table: pd.DataFrame,
    all_spidroins: pd.DataFrame,
    outprefix: str
) -> None:
    """
    Plot the second best hit after unclassified spidroin.

    Args:
        spidroin_table: DataFrame with filtered spidroin data
        all_spidroins: DataFrame with all spidroin hits
        outprefix: Output file prefix
    """
    intervals = get_spidroin_intervals(all_spidroins)
    unclassified = spidroin_table.loc[spidroin_table["spidroin_type"] == "Spidroin"]

    second_hit = []
    no_secondbest = 0

    for _, hit in unclassified.iterrows():
        duplicates = get_overlapping_spidroins(all_spidroins, hit, intervals)
        duplicates.sort_values("E-value", inplace=True)
        duplicates = duplicates.loc[duplicates["spidroin_type"] != "Spidroin"]

        if len(duplicates) > 0:
            second_hit.append([
                hit["E-value"],
                duplicates.iloc[0]["spidroin_type"],
                duplicates.iloc[0]["E-value"]
            ])
        else:
            no_secondbest += 1
            second_hit.append(["None", "None", "None"])

    second_hit = pd.DataFrame(second_hit, columns=["uncl_E-value", "spidroin_type", "E-value"])
    second_hit["difference"] = second_hit["E-value"] - second_hit["uncl_E-value"]

    second_hit["difference"].loc[second_hit["difference"] != "None"].plot(
        kind="box", colormap="Paired", logy=True
    )
    plt.xlabel("E-value")
    plt.ylabel("Counts")
    plt.title("E-value difference between unclassified and best hit")
    plt.savefig(f"{outprefix}_evalue_barplot.png", bbox_inches="tight")
    plt.savefig(f"{outprefix}_evalue_barplot.svg", bbox_inches="tight")
    plt.clf()

    second_hit["spidroin_type"].value_counts().plot(kind="bar", colormap="Paired")
    plt.xlabel("Spidroin Type")
    plt.ylabel("Counts")
    plt.title("Second best fitting spidroin type for unclassified spidrions")
    plt.savefig(f"{outprefix}_spidroin_type_barplot.png", bbox_inches="tight")
    plt.savefig(f"{outprefix}_spidroin_type_barplot.svg", bbox_inches="tight")
    plt.clf()
    print(f"Number of hits without second best: {no_secondbest}")


def join_terminals(spidroin_table: pd.DataFrame) -> tuple:
    """
    Join N-terminal and C-terminal domains to identify full spidroin genes.

    Args:
        spidroin_table: DataFrame with spidroin hits

    Returns:
        Tuple of (full_spidroins DataFrame, singletons DataFrame)
    """
    n_terminals = spidroin_table.loc[spidroin_table["domain"] == "NTD"]
    c_terminals = spidroin_table.loc[spidroin_table["domain"] == "CTD"]
    full_spidroins = []
    singletons = []

    for _, ntd in n_terminals.iterrows():
        ctds = c_terminals.loc[
            (c_terminals["Species"] == ntd["Species"]) &
            (c_terminals["target_name"] == ntd["target_name"]) &
            (c_terminals["strand"] == ntd["strand"])
        ]
        if len(ctds) > 0:
            ctds["distance"] = abs(ctds["ali_from"] - ntd["ali_from"])
            ctds = ctds.sort_values("distance", ascending=True)
            ctd = ctds.iloc[0]
            if (
                (ctd["ali_from"] > ntd["ali_from"]) & (ntd["strand"] == "+")
            ) or (
                (ctd["ali_from"] < ntd["ali_from"]) & (ntd["strand"] == "-")
            ):
                gene_borders = [
                    ntd["ali_from"], ntd["ali_to"], ctd["ali_from"], ctd["ali_to"]
                ]
                gene_length = max(gene_borders) - min(gene_borders)
                full_spidroins.append([
                    ntd["Species"], ntd["target_name"], ntd["strand"],
                    ntd["ali_from"], ntd["ali_to"], ctd["ali_from"], ctd["ali_to"],
                    gene_length, ntd["E-value"], ctd["E-value"],
                    ntd["spidroin_type"], ctd["spidroin_type"]
                ])
            else:
                singletons.append(ntd)
        else:
            singletons.append(ntd)

    full_spidroins = pd.DataFrame(full_spidroins, columns=[
        "Species", "Scaffold", "Strand", "NTD_from", "NTD_to",
        "CTD_from", "CTD_to", "Length", "NTD_e-value", "CTD_e-value",
        "NTD_type", "CTD_type"
    ])

    for _, ctd in c_terminals.iterrows():
        hits = full_spidroins.loc[
            (full_spidroins["Scaffold"] == ctd["target_name"]) &
            (full_spidroins["Strand"] == ctd["strand"]) &
            (full_spidroins["CTD_from"] == ctd["ali_from"]) &
            (full_spidroins["CTD_to"] == ctd["ali_to"])
        ]
        if len(hits) == 0:
            singletons.append(ctd)

    singletons = pd.DataFrame(singletons, columns=spidroin_table.columns.to_list())
    return full_spidroins, singletons


def assign_spidroin_type(spidroin_table: pd.DataFrame) -> list:
    """
    Assign spidroin type based on NTD and CTD classifications.

    Args:
        spidroin_table: DataFrame with joined spidroin data

    Returns:
        List of assigned spidroin types
    """
    maj_ampullates = ["MaSp", "MaSp1", "MaSp2", "MaSp2B", "MaSp3", "MaSp3B"]
    assignments = []

    for _, spid in spidroin_table.iterrows():
        if spid["NTD_type"] == spid["CTD_type"]:
            assignments.append(spid["NTD_type"])
        elif (
            (spid["NTD_type"] == "Spidroin" and spid["CTD_type"] == "hypo") or
            (spid["CTD_type"] == "Spidroin" and spid["NTD_type"] == "hypo")
        ):
            assignments.append("Unclassified")
        elif spid["NTD_type"] == "Spidroin" or spid["CTD_type"] == "Spidroin":
            terminals = spid[["CTD_type", "NTD_type"]].values.tolist()
            terminals.remove("Spidroin")
            assignments.append(terminals[0])
        elif spid["NTD_type"] == "hypo" or spid["CTD_type"] == "hypo":
            terminals = spid[["CTD_type", "NTD_type"]].values.tolist()
            terminals.remove("hypo")
            assignments.append(terminals[0])
        elif spid["NTD_type"] == "MaSp" and spid["CTD_type"] in maj_ampullates:
            assignments.append(spid["CTD_type"])
        elif spid["CTD_type"] == "MaSp" and spid["NTD_type"] in maj_ampullates:
            assignments.append(spid["CTD_type"])
        else:
            assignments.append("Discordant")

    assignments = [
        "Unclassified" if spid_type == "Spidroin" else spid_type
        for spid_type in assignments
    ]

    return assignments


def plot_type_discordance(spidroin_table: pd.DataFrame, outprefix: str) -> None:
    """
    Plot spidroin type distribution pie chart.

    Args:
        spidroin_table: DataFrame with spidroin data
        outprefix: Output file prefix
    """
    types = spidroin_table.value_counts("spidroin_type")
    types.plot(kind="pie", title="Spidroin Type assignments", autopct='%1.1f%%')
    plt.savefig(f"{outprefix}spidroin_types.png", bbox_inches="tight")
    plt.savefig(f"{outprefix}spidroin_types.svg", bbox_inches="tight")
    plt.clf()


def plot_spidroin_length(
    spidroin_table: pd.DataFrame,
    outprefix: str,
    order: list
) -> None:
    """
    Plot spidroin length distribution by type.

    Args:
        spidroin_table: DataFrame with spidroin data
        outprefix: Output file prefix
        order: Order of spidroin types for plotting
    """
    sns.set_theme(style="ticks")

    f, ax = plt.subplots(figsize=(7, 6))

    sns.boxplot(
        x="Length", y="spidroin_type", data=spidroin_table,
        hue="spidroin_type", legend=False, order=order, palette="vlag"
    )

    ax.xaxis.grid(True)
    ax.set(ylabel="")
    sns.despine(trim=False, left=True)

    ax.set_xlim(0, 100000)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    plt.ylabel("Spidroin Type")
    plt.xlabel("Length [bp]")
    plt.title("Gene Lengths of different Spidroin Classes")

    plt.savefig(f"{outprefix}spidroin_lengths.png", bbox_inches="tight")
    plt.savefig(f"{outprefix}spidroin_lengths.svg", bbox_inches="tight")

    plt.clf()


def write_gff_files(spidroin_data: pd.DataFrame, output_folder: str) -> None:
    """
    Write GFF files for spidroin annotations, one file per species.

    Args:
        spidroin_data: DataFrame with joined spidroin data
        output_folder: Output directory for GFF files
    """
    gff_dir = os.path.join(output_folder, "GFF")
    os.makedirs(gff_dir, exist_ok=True)

    for species in spidroin_data["Species"].unique():
        species_data = spidroin_data[spidroin_data["Species"] == species]
        gff_file = os.path.join(gff_dir, f"{species}.gff")

        type_counters = {}
        with open(gff_file, "w") as f:
            f.write("##gff-version 3\n")

            for _, entry in species_data.iterrows():
                scaffold = entry["Scaffold"]
                strand = entry["Strand"]
                spidroin_type = entry["spidroin_type"]

                # Track type counter for unique gene IDs
                if spidroin_type not in type_counters:
                    type_counters[spidroin_type] = 1
                else:
                    type_counters[spidroin_type] += 1

                gene_id = f"{species}_{spidroin_type}_{type_counters[spidroin_type]}"

                # Gene boundaries
                gene_start = int(min(entry["NTD_from"], entry["NTD_to"],
                                     entry["CTD_from"], entry["CTD_to"]))
                gene_end = int(max(entry["NTD_from"], entry["NTD_to"],
                                   entry["CTD_from"], entry["CTD_to"]))

                # NTD boundaries
                ntd_start = int(min(entry["NTD_from"], entry["NTD_to"]))
                ntd_end = int(max(entry["NTD_from"], entry["NTD_to"]))

                # CTD boundaries
                ctd_start = int(min(entry["CTD_from"], entry["CTD_to"]))
                ctd_end = int(max(entry["CTD_from"], entry["CTD_to"]))

                # Write gene feature
                attrs = f"ID={gene_id};Name={spidroin_type};spidroin_type={spidroin_type}"
                f.write(f"{scaffold}\tspidroin_analysis\tgene\t{gene_start}\t{gene_end}\t.\t{strand}\t.\t{attrs}\n")

                # Write NTD feature
                ntd_attrs = f"ID={gene_id}_NTD;Parent={gene_id};Name=NTD"
                f.write(f"{scaffold}\tspidroin_analysis\tCDS\t{ntd_start}\t{ntd_end}\t.\t{strand}\t.\t{ntd_attrs}\n")

                # Write CTD feature
                ctd_attrs = f"ID={gene_id}_CTD;Parent={gene_id};Name=CTD"
                f.write(f"{scaffold}\tspidroin_analysis\tCDS\t{ctd_start}\t{ctd_end}\t.\t{strand}\t.\t{ctd_attrs}\n")

        print(f"Wrote GFF file: {gff_file}")


def extract_spidroin_sequences(
    spidroin_data: pd.DataFrame,
    genome_fasta: str,
    species: str,
    output_folder: str = "spidroin_analysis",
    nostop_file: str = "no_stop_spidroins.txt",
    nostart_file: str = "no_start_spidroins.txt"
) -> tuple:
    """
    Extract spidroin sequences from genome assembly.

    Args:
        spidroin_data: DataFrame with spidroin coordinates
        genome_fasta: Path to genome FASTA file
        species: Species name for output naming
        output_folder: Output directory
        nostop_file: File to log sequences without stop codons
        nostart_file: File to log sequences without start codons

    Returns:
        Tuple of (ntd_sequences, ctd_sequences, ntd_proteins, ctd_proteins, spidroin_sequences)
    """
    print(f"Assembly: {genome_fasta}")
    os.makedirs(output_folder, exist_ok=True)

    genome_sequences = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))

    ntd_sequences = []
    ctd_sequences = []
    ntd_proteins = []
    ctd_proteins = []
    spidroin_sequences = []
    type_counters = {}

    for _, entry in spidroin_data.iterrows():
        start_extension = 45
        stop_extension = 45
        scaffold = str(entry["Scaffold"])
        strand = entry["Strand"]
        scaff_seq = genome_sequences[scaffold]
        scaff_len = len(scaff_seq)
        print(f"{scaffold}: {scaff_len} bp")

        ntd_borders = entry[["NTD_from", "NTD_to"]].values.tolist()
        ntd_from, ntd_to = [min(ntd_borders), max(ntd_borders)]

        ctd_borders = entry[["CTD_from", "CTD_to"]].values.tolist()
        ctd_from, ctd_to = [min(ctd_borders), max(ctd_borders)]

        spidroin_borders = entry[["NTD_from", "NTD_to", "CTD_from", "CTD_to"]].values.tolist()
        spidroin_from, spidroin_to = [min(spidroin_borders), max(spidroin_borders)]

        spidroin_type = entry["spidroin_type"]

        if spidroin_type not in type_counters:
            type_counters[spidroin_type] = 1
        else:
            type_counters[spidroin_type] += 1

        if strand == "+":
            if ntd_from < start_extension:
                start_extension = ntd_from
            if (scaff_len - ctd_to) < stop_extension:
                stop_extension = scaff_len - ctd_to

            ntd_sequence = scaff_seq.seq[ntd_from - start_extension:ntd_to].upper()
            ctd_sequence = scaff_seq.seq[ctd_from:ctd_to + stop_extension].upper()

            try:
                start_index = min(re.search("ATG", str(ntd_sequence)).span())
            except AttributeError:
                print(f"Warning! Did not find a start codon for {species}_{spidroin_type}_{type_counters[spidroin_type]}")
                start_index = start_extension
                with open(os.path.join(output_folder, nostart_file), "a") as file:
                    file.write(f"{species}_{spidroin_type}_{type_counters[spidroin_type]}\n")
            start_offset = start_extension - start_index

            inv_ctd = ctd_sequence[::-1]
            ctd_len = ctd_to - ctd_from
            try:
                inv_index = min(re.search(r"GAT|AGT|AAT", str(inv_ctd)).span())
                stop_index = len(ctd_sequence) - inv_index
            except AttributeError:
                print(f"Warning! Did not find a stop codon for {species}_{spidroin_type}_{type_counters[spidroin_type]}")
                with open(os.path.join(output_folder, nostop_file), "a") as file:
                    file.write(f"{species}_{spidroin_type}_{type_counters[spidroin_type]}\n")
                stop_index = ctd_len

            stop_offset = (stop_index - ctd_len) + start_offset - 1

            ntd_sequence = ntd_sequence[start_index:]
            ctd_sequence = ctd_sequence[:stop_index]

            spidroin_sequence = scaff_seq.seq[
                spidroin_from - start_offset:spidroin_to + stop_offset
            ].upper()

        else:
            if ctd_from < start_extension:
                start_extension = ctd_from
            if (scaff_len - ntd_to) < stop_extension:
                stop_extension = scaff_len - ntd_to

            ntd_sequence = scaff_seq.seq[ntd_from:ntd_to + stop_extension].reverse_complement().upper()
            ctd_sequence = scaff_seq.seq[ctd_from - start_extension:ctd_to].reverse_complement().upper()

            try:
                start_index = min(re.search("ATG", str(ntd_sequence)).span())
            except AttributeError:
                print(f"Warning! Did not find a start codon for {species}_{spidroin_type}_{type_counters[spidroin_type]}")
                start_index = start_extension
                with open(os.path.join(output_folder, nostart_file), "a") as file:
                    file.write(f"{species}_{spidroin_type}_{type_counters[spidroin_type]}\n")

            inv_ctd = ctd_sequence[::-1]
            ctd_len = ctd_to - ctd_from
            try:
                inv_index = min(re.search(r"GAT|AGT|AAT", str(inv_ctd)).span())
                stop_index = len(ctd_sequence) - inv_index
            except AttributeError:
                print(f"Warning! Did not find a stop codon for {species}_{spidroin_type}_{type_counters[spidroin_type]}")
                with open(os.path.join(output_folder, nostop_file), "a") as file:
                    file.write(f"{species}_{spidroin_type}_{type_counters[spidroin_type]}\n")
                stop_index = ctd_len

            start_offset = start_extension - start_index
            stop_offset = (stop_index - ctd_len) + start_offset

            ntd_sequence = ntd_sequence[start_index:]
            ctd_sequence = ctd_sequence[:stop_index]

            spidroin_sequence = scaff_seq.seq[
                spidroin_from - stop_offset:spidroin_to + start_offset
            ].reverse_complement().upper()

        seq_id_suffix = f"{species}_{spidroin_type}_{type_counters[spidroin_type]}"

        ntd_sequences.append(SeqRecord(
            ntd_sequence,
            id=f"{seq_id_suffix}_NTD",
            description=f"{scaffold}_{strand}_{ntd_from}_{ntd_to}"
        ))
        ctd_sequences.append(SeqRecord(
            ctd_sequence,
            id=f"{seq_id_suffix}_CTD",
            description=f"{scaffold}_{strand}_{ctd_from}_{ctd_to}"
        ))

        # Translate to protein sequences
        # NTD: trim from end to make length multiple of 3 (keep start codon)
        # CTD: trim from start to make length multiple of 3 (keep stop codon)
        # Note: use to_stop=False to translate full region (internal stops shown as *)
        try:
            ntd_trim_len = len(ntd_sequence) - (len(ntd_sequence) % 3)
            ntd_for_translate = ntd_sequence[:ntd_trim_len]
            ntd_protein = ntd_for_translate.translate(to_stop=False)
            ntd_proteins.append(SeqRecord(
                ntd_protein,
                id=f"{seq_id_suffix}_NTD",
                description=f"{scaffold}_{strand}_{ntd_from}_{ntd_to}"
            ))
        except Exception as e:
            print(f"Warning: Failed to translate NTD for {seq_id_suffix}: {e}")

        try:
            ctd_remainder = len(ctd_sequence) % 3
            ctd_for_translate = ctd_sequence[ctd_remainder:] if ctd_remainder else ctd_sequence
            ctd_protein = ctd_for_translate.translate(to_stop=False)
            ctd_proteins.append(SeqRecord(
                ctd_protein,
                id=f"{seq_id_suffix}_CTD",
                description=f"{scaffold}_{strand}_{ctd_from}_{ctd_to}"
            ))
        except Exception as e:
            print(f"Warning: Failed to translate CTD for {seq_id_suffix}: {e}")

        spidroin_sequences.append(SeqRecord(
            spidroin_sequence,
            id=seq_id_suffix,
            description=f"{scaffold}_{strand}_{spidroin_from}_{spidroin_to}"
        ))

    return ntd_sequences, ctd_sequences, ntd_proteins, ctd_proteins, spidroin_sequences


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Spidroin Analysis Pipeline - Analyze spider silk protein data from HMMER search results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Single species mode (recommended for single genome analysis)
    python analyse_spidroins.py -i nhmmer_search/species -p hmmer_profiles -a genome.fa -s Species_name

    # Multi-species mode (batch analysis)
    python analyse_spidroins.py -i nhmmer_search_trimmed -p hmmer_profiles --assembly-dir genomes

    # Custom E-value threshold
    python analyse_spidroins.py -i nhmmer_search/species -p hmmer_profiles -a genome.fa -s Species -e 0.001

    # Skip sequence extraction (only generate tables and plots)
    python analyse_spidroins.py -i nhmmer_search/species -p hmmer_profiles -a genome.fa -s Species --skip-sequence-extraction
        """
    )

    # Required arguments
    parser.add_argument(
        "-i", "--input-dir",
        required=True,
        help="Directory containing HMMER search results (.tbl files)"
    )
    parser.add_argument(
        "-p", "--profile-dir",
        required=True,
        help="Directory containing HMM profile files (*.hmm)"
    )

    # Genome input (mutually exclusive: single file or directory)
    genome_group = parser.add_mutually_exclusive_group(required=True)
    genome_group.add_argument(
        "-a", "--assembly-file",
        help="Single genome assembly file (FASTA format) - use with -s/--species"
    )
    genome_group.add_argument(
        "--assembly-dir",
        help="Directory containing genome assemblies (pattern: assembly_dir/species/*_genomic.fna) - for batch mode"
    )

    parser.add_argument(
        "-s", "--species",
        help="Species name (required when using -a/--assembly-file)"
    )

    # Optional arguments
    parser.add_argument(
        "-o", "--output-dir",
        default="spidroin_analysis",
        help="Output directory (default: spidroin_analysis)"
    )
    parser.add_argument(
        "-e", "--e-value",
        type=float,
        default=0.0001,
        help="E-value threshold for filtering (default: 0.0001)"
    )
    parser.add_argument(
        "--hmm-length-factor",
        type=float,
        default=0.9,
        help="Minimum HMM profile coverage factor (default: 0.9)"
    )
    parser.add_argument(
        "--max-spidroin-length",
        type=int,
        default=100000,
        help="Maximum spidroin gene length in bp (default: 100000)"
    )

    # Boolean flags
    parser.add_argument(
        "--no-summary",
        action="store_true",
        help="Skip summary generation, use existing summary table"
    )
    parser.add_argument(
        "--no-filter",
        action="store_true",
        help="Skip filtering, use existing filtered table"
    )
    parser.add_argument(
        "--no-join",
        action="store_true",
        help="Skip domain joining, use existing joined table"
    )
    parser.add_argument(
        "--skip-plots",
        action="store_true",
        help="Skip plot generation"
    )
    parser.add_argument(
        "--skip-sequence-extraction",
        action="store_true",
        help="Skip spidroin sequence extraction"
    )

    return parser.parse_args()


def main():
    """Main function to run the spidroin analysis pipeline."""
    args = parse_arguments()

    # Validate arguments
    if args.assembly_file and not args.species:
        print("Error: --species (-s) is required when using --assembly-file (-a)")
        return 1

    # Determine mode: single species or batch
    single_species_mode = args.assembly_file is not None

    # Setup
    outdir = args.output_dir.rstrip("/") + "/"
    os.makedirs(outdir, exist_ok=True)

    spidroin_plotting_order = [
        "MiSp", "MaSp", "CrSp", "Pflag", "Flag",
        "AgSp", "AcSp", "CySp", "PySp", "Discordant", "Unclassified"
    ]

    column_names = [
        "target_name", "target_accession", "query_name", "query_accession",
        "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to",
        "sq_len", "strand", "E-value", "score", "bias", "Species"
    ]
    new_order = [
        "Species", "target_name", "strand", "ali_from", "ali_to",
        "spidroin_type", "domain", "E-value", "score", "bias",
        "hmm_from", "hmm_to", "env_from", "env_to", "sq_len", "query_name"
    ]

    generate_summary = not args.no_summary
    generate_filtered = not args.no_filter
    join_spidroins = not args.no_join

    # Generate or load summary table
    if generate_summary:
        print("Generating summary table from HMMER results...")
        summary_table = pd.DataFrame(columns=new_order)

        if single_species_mode:
            # Single species mode: read tbl files directly from input directory
            species = args.species
            hmmer_tbls = glob.glob(f"{args.input_dir}/*.tbl")
            print(f"Reading HMMER-Tables for {species}")

            for table in hmmer_tbls:
                summary_table = pd.concat([
                    summary_table,
                    read_hmmer_tbl(table, column_names, species)[new_order]
                ], ignore_index=True)
        else:
            # Batch mode: read from subdirectories
            hmmer_dirs = glob.glob(f"{args.input_dir}/*")
            for dir_path in hmmer_dirs:
                species = dir_path.split("/")[-1]
                hmmer_tbls = glob.glob(dir_path + "/*.tbl")
                print(f"Reading HMMER-Tables for {species}")

                for table in hmmer_tbls:
                    summary_table = pd.concat([
                        summary_table,
                        read_hmmer_tbl(table, column_names, species)[new_order]
                    ], ignore_index=True)

        print("Finished reading HMMER-Tables")
        print("Writing summary table")
        summary_table.to_csv(outdir + "spidroins_total.tsv", sep="\t", index=False)
    else:
        print("Using pre-computed summary table")
        summary_table = pd.read_csv(outdir + "spidroins_total.tsv", sep="\t")

    # Get HMM profile lengths
    print("Getting HMM profile lengths")
    profile_lengths = get_hmm_length(args.profile_dir)

    # Filter spidroins
    if generate_filtered:
        print("Applying filters...")
        filtered_spidroins = filter_spidroins(
            summary_table, args.e_value, profile_lengths, args.hmm_length_factor
        )

        if len(filtered_spidroins) > 0:
            filtered_spidroins.sort_values(by=["target_name", "query_name"], inplace=True)
            filtered_spidroins["spidroin_type"] = (
                filtered_spidroins["spidroin_type"]
                .str.rstrip("B")
                .str.replace(r"\d+", "", regex=True)
                .str.rstrip("B")
            )
        filtered_spidroins.to_csv(outdir + "spidroins_filtered.tsv", sep="\t", index=False)
        print(f"Found {len(filtered_spidroins)} spidroin hits after filtering")
    else:
        print("Using pre-filtered spidroin table")
        filtered_spidroins = pd.read_csv(outdir + "spidroins_filtered.tsv", sep="\t")

    if len(filtered_spidroins) == 0:
        print("Warning: No spidroin hits passed the filters!")
        print("Consider adjusting --e-value or --hmm-length-factor parameters.")
        print("Analysis complete (no results).")
        return

    # Filter excluding unclassified spidroins
    print("Removing Spidroin hits originating from unclassified spidroins")
    no_unclassified = summary_table.loc[summary_table["spidroin_type"] != "Spidroin"]
    filtered_no_unclassified = filter_spidroins(
        no_unclassified, args.e_value, profile_lengths, args.hmm_length_factor
    )

    if len(filtered_no_unclassified) > 0:
        filtered_no_unclassified.sort_values(by=["target_name", "query_name"], inplace=True)
    filtered_no_unclassified.to_csv(
        outdir + "spidroins_filtered_no_unclassified.tsv", sep="\t", index=False
    )

    # Generate plots
    if not args.skip_plots and len(filtered_spidroins) > 0:
        print("Generating plots...")
        all_species = sorted(summary_table["Species"].unique().tolist())
        plot_spidroins_heatmap(filtered_spidroins, outdir + "all_spidroins", spidroin_plotting_order, all_species)
        if len(filtered_no_unclassified) > 0:
            plot_spidroins_heatmap(filtered_no_unclassified, outdir + "no_unclassified", spidroin_plotting_order)
        plot_domain_ratio(filtered_spidroins, outdir + "domains")
        plot_domain_evals(filtered_spidroins, outdir + "domain")

        only_NTD = filtered_spidroins.loc[filtered_spidroins["domain"] == "NTD"]
        if len(only_NTD) > 0:
            plot_spidroins_heatmap(only_NTD, outdir + "only_NTD", spidroin_plotting_order)

        only_CTD = filtered_spidroins.loc[filtered_spidroins["domain"] == "CTD"]
        if len(only_CTD) > 0:
            plot_spidroins_heatmap(only_CTD, outdir + "only_CTD", spidroin_plotting_order)

        plot_second_best_classified(
            filtered_spidroins, summary_table, outdir + "second_best_classified"
        )

    # Join terminal domains
    if join_spidroins:
        print("Joining Terminal Domains...")
        full_spidroin_genes, singletons = join_terminals(filtered_spidroins)

        if len(full_spidroin_genes) > 0:
            full_spidroin_genes = full_spidroin_genes.loc[
                full_spidroin_genes["Length"] <= args.max_spidroin_length
            ]
            full_spidroin_genes["spidroin_type"] = assign_spidroin_type(full_spidroin_genes)
            full_spidroin_genes["spidroin_type"] = (
                full_spidroin_genes["spidroin_type"]
                .str.rstrip("B")
                .str.replace(r"\d+", "", regex=True)
                .str.rstrip("B")
            )
        full_spidroin_genes.to_csv(outdir + "joined_domains.tsv", index=False, sep="\t")
        singletons.to_csv(outdir + "singletons.tsv", sep="\t", index=False)
        print(f"Found {len(full_spidroin_genes)} joined spidroin genes, {len(singletons)} singletons")
    else:
        print("Using pre-joined spidroin table")
        full_spidroin_genes = pd.read_csv(outdir + "joined_domains.tsv", sep="\t")
        singletons = pd.read_csv(outdir + "singletons.tsv", sep="\t")

    # Generate additional plots
    all_species = sorted(summary_table["Species"].unique().tolist())
    if not args.skip_plots:
        if len(singletons) > 0:
            plot_spidroins_heatmap(singletons, outdir + "singletons", spidroin_plotting_order, all_species)
        if len(full_spidroin_genes) > 0:
            plot_type_discordance(full_spidroin_genes, outdir)
            plot_spidroin_length(full_spidroin_genes, outdir, spidroin_plotting_order)
            plot_spidroins_heatmap(full_spidroin_genes, outdir + "joined", spidroin_plotting_order, all_species)

    # Generate GFF files
    if len(full_spidroin_genes) > 0:
        print("Generating GFF files...")
        write_gff_files(full_spidroin_genes, outdir)

    # Extract spidroin sequences
    if not args.skip_sequence_extraction and len(full_spidroin_genes) > 0:
        print("Extracting spidroin sequences...")
        ntd_seqs = []
        ctd_seqs = []
        ntd_prots = []
        ctd_prots = []
        spid_seqs = []

        if single_species_mode:
            # Single species mode: use the provided assembly file
            species = args.species
            print(f"Extracting spidroin sequences for {species}")
            species_spidroins = full_spidroin_genes.loc[
                full_spidroin_genes["Species"] == species
            ]
            if len(species_spidroins) > 0:
                spec_ntd_seqs, spec_ctd_seqs, spec_ntd_prots, spec_ctd_prots, spec_spid_seqs = extract_spidroin_sequences(
                    species_spidroins, args.assembly_file, species, outdir
                )
                ntd_seqs.extend(spec_ntd_seqs)
                ctd_seqs.extend(spec_ctd_seqs)
                ntd_prots.extend(spec_ntd_prots)
                ctd_prots.extend(spec_ctd_prots)
                spid_seqs.extend(spec_spid_seqs)
            else:
                print(f"Warning: No joined spidroin genes found for {species}")
        else:
            # Batch mode: find assembly files by species name
            for species in full_spidroin_genes["Species"].unique():
                print(f"Extracting spidroin sequences for {species}")
                species_spidroins = full_spidroin_genes.loc[
                    full_spidroin_genes["Species"] == species
                ]
                # Try multiple patterns to find assembly file
                assembly_files = glob.glob(f"{args.assembly_dir}/{species}/*_genomic.fna")
                if not assembly_files:
                    assembly_files = glob.glob(f"{args.assembly_dir}/{species}/*.fa")
                if not assembly_files:
                    assembly_files = glob.glob(f"{args.assembly_dir}/{species}/*.fasta")
                if not assembly_files:
                    assembly_files = glob.glob(f"{args.assembly_dir}/{species}.fa")
                if not assembly_files:
                    assembly_files = glob.glob(f"{args.assembly_dir}/{species}.fasta")
                if not assembly_files:
                    print(f"Warning: No assembly file found for {species}")
                    continue

                assembly_file = assembly_files[0]
                spec_ntd_seqs, spec_ctd_seqs, spec_ntd_prots, spec_ctd_prots, spec_spid_seqs = extract_spidroin_sequences(
                    species_spidroins, assembly_file, species, outdir
                )

                ntd_seqs.extend(spec_ntd_seqs)
                ctd_seqs.extend(spec_ctd_seqs)
                ntd_prots.extend(spec_ntd_prots)
                ctd_prots.extend(spec_ctd_prots)
                spid_seqs.extend(spec_spid_seqs)

        # Write nucleotide sequences
        SeqIO.write(ntd_seqs, f"{outdir}/ntd_sequences.fasta", "fasta")
        SeqIO.write(ctd_seqs, f"{outdir}/ctd_sequences.fasta", "fasta")
        SeqIO.write(spid_seqs, f"{outdir}/spidroin_sequences.fasta", "fasta")

        # Write protein sequences
        SeqIO.write(ntd_prots, f"{outdir}/ntd_proteins.fasta", "fasta")
        SeqIO.write(ctd_prots, f"{outdir}/ctd_proteins.fasta", "fasta")

    print("Analysis complete!")


if __name__ == "__main__":
    main()
