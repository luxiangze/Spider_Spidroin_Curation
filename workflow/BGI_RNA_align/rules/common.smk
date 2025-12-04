# Common functions and variables for the workflow

import os
import glob
from pathlib import Path

# Get workflow base directory for resolving relative paths
WORKFLOW_DIR = Path(workflow.basedir)


def resolve_path(relative_path):
    """Resolve relative path to absolute path based on workflow directory."""
    path = Path(relative_path)
    if path.is_absolute():
        return str(path)
    return str((WORKFLOW_DIR / path).resolve())


def get_genomes():
    """Get list of genomes to process."""
    return list(config["samples"].keys())


def get_genome_fasta(wildcards):
    """Get genome fasta file path for a given genome."""
    # Genome files are organized in subfolders: genome_dir/genome_name/genome_name.fa
    return resolve_path(os.path.join(
        config["genome_dir"], wildcards.genome, f"{wildcards.genome}.fa"
    ))


def get_genome_gff(wildcards):
    """Get genome annotation gff file path for a given genome."""
    # GFF files are organized in subfolders: genome_dir/genome_name/genome_name.gff
    return resolve_path(os.path.join(
        config["genome_dir"], wildcards.genome, f"{wildcards.genome}.gff"
    ))


def get_fastq_folder(genome):
    """Get fastq folder path for a genome."""
    folder = config["samples"][genome]["fastq_folder"]
    return resolve_path(os.path.join(config["fastq_dir"], folder))


def get_fastq_files(genome):
    """Get R1 and R2 fastq files for a genome."""
    folder = get_fastq_folder(genome)
    r1_pattern = os.path.join(folder, f"*{config['fastq_suffix'][0]}")
    r2_pattern = os.path.join(folder, f"*{config['fastq_suffix'][1]}")
    r1_files = glob.glob(r1_pattern)
    r2_files = glob.glob(r2_pattern)
    if r1_files and r2_files:
        return sorted(r1_files)[0], sorted(r2_files)[0]
    return None, None


def get_fastq_r1(wildcards):
    """Get R1 fastq file for a genome."""
    r1, _ = get_fastq_files(wildcards.genome)
    return r1


def get_fastq_r2(wildcards):
    """Get R2 fastq file for a genome."""
    _, r2 = get_fastq_files(wildcards.genome)
    return r2
