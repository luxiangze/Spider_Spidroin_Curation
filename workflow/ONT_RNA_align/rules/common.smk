# Common functions and variables for the ONT RNA workflow

import os
import glob
from pathlib import Path

# Get workflow base directory for resolving relative paths
# workflow.basedir points to the directory containing the Snakefile
WORKFLOW_DIR = Path(workflow.basedir).parent if Path(workflow.basedir).name == "rules" else Path(workflow.basedir)


# Wildcard constraints to avoid ambiguous rules
wildcard_constraints:
    sample="[^.]+"  # Sample name cannot contain dots


def resolve_path(relative_path):
    """Resolve relative path to absolute path based on workflow directory."""
    path = Path(relative_path)
    if path.is_absolute():
        return str(path)
    return str((WORKFLOW_DIR / path).resolve())


def get_samples():
    """Get list of samples to process."""
    return list(config["samples"].keys())


def get_genome_for_sample(sample):
    """Get the genome name for a given sample."""
    return config["samples"][sample]["genome"]


def get_genome_fasta(wildcards):
    """Get genome fasta file path for a given sample."""
    genome = get_genome_for_sample(wildcards.sample)
    return os.path.join(config["genome_dir"], genome, genome + ".fa")


def get_genome_gff(wildcards):
    """Get genome annotation gff file path for a given sample."""
    genome = get_genome_for_sample(wildcards.sample)
    return os.path.join(config["genome_dir"], genome, genome + ".gff")


def get_fastq(wildcards):
    """Get fastq file path for a given sample."""
    return os.path.join(config["fastq_dir"], wildcards.sample, "pass.fq.gz")


def get_star_sj_tab(wildcards):
    """Get STAR SJ.out.tab file from BGI RNA-seq alignment."""
    genome = get_genome_for_sample(wildcards.sample)
    return os.path.join(config["bgi_star_dir"], genome, "SJ.out.tab")
