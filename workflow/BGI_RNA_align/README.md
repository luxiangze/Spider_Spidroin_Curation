# BGI RNA Alignment Workflow

Multi-genome RNA-seq alignment pipeline using Snakemake and STAR.

## Overview

This workflow aligns RNA-seq reads from BGI sequencing to 10 spider genomes. Each genome has its corresponding RNA-seq sample.

### Pipeline Steps

1. **Quality Control**: FastQC (raw) → fastp trimming → FastQC (trimmed)
2. **Index Building**: STAR genomeGenerate
3. **Alignment**: STAR alignReads → sorted BAM
4. **Post-processing**: samtools index → BigWig generation
5. **Statistics**: samtools flagstat/stats + MultiQC reports

## Data Structure

```
data/raw/
├── spider_genome/                    # Genome files
│   ├── Araneus_ventricosus/
│   │   ├── Araneus_ventricosus.fa
│   │   └── Araneus_ventricosus.gff
│   └── ...
└── BGI_RNA_10samples/                # RNA-seq data
    ├── Araneus_ventricosus/
    │   ├── *_R1.fq.gz
    │   └── *_R2.fq.gz
    └── ...
```

## Workflow Structure

```
BGI_RNA_align/
├── Snakefile           # Main workflow
├── config/
│   └── config.yaml     # Configuration (genomes, threads, memory)
├── rules/
│   ├── common.smk      # Helper functions
│   ├── qc.smk          # FastQC + fastp
│   ├── index.smk       # STAR index
│   ├── align.smk       # STAR alignment
│   ├── stats.smk       # Statistics
│   └── bigwig.smk      # BigWig generation
├── work/               # Intermediate files
└── results/            # Final outputs
```

## Requirements

All dependencies are managed by pixi:

- Snakemake >= 7.0
- STAR >= 2.7
- samtools >= 1.10
- FastQC >= 0.12
- fastp >= 1.0
- deeptools >= 3.5
- MultiQC >= 1.10

## Usage

### 1. Dry run

```bash
cd workflow/BGI_RNA_align
pixi run snakemake -n
```

### 2. Build indices only

```bash
pixi run snakemake build_all_indices --cores 80 --resources mem_gb=240
```

### 3. Run full pipeline

```bash
# Recommended: with memory limit (for 252GB server)
pixi run snakemake --cores 80 --resources mem_gb=240

# Without memory limit (risky if STAR jobs overlap)
pixi run snakemake --cores 80
```

### 4. Process specific genome

```bash
pixi run snakemake results/.done.Araneus_ventricosus --cores 80 --resources mem_gb=240
```

### 5. Generate workflow DAG

```bash
pixi run snakemake --dag | dot -Tpng > dag.png
```

## Resource Configuration

Default settings optimized for 80-core / 252GB server:

| Task | Threads | Memory | Max Parallel |
|------|---------|--------|--------------|
| star_index | 24 | 100GB | 2 |
| star_align | 12 | 50GB | 5 |
| fastp | 4 | - | 20 |
| bigwig | 8 | - | 10 |

Edit `config/config.yaml` to adjust for your server.

## Output

```
results/
├── bam/
│   ├── Araneus_ventricosus.bam
│   ├── Araneus_ventricosus.bam.bai
│   └── ...
├── bigwig/
│   ├── Araneus_ventricosus.bw
│   └── ...
├── stats/
│   ├── Araneus_ventricosus.flagstat.txt
│   ├── Araneus_ventricosus.stats.txt
│   ├── Araneus_ventricosus.star_log.txt
│   └── ...
└── multiqc/
    ├── raw_reads_multiqc.html      # Pre-trimming QC
    ├── trimmed_reads_multiqc.html  # Post-trimming QC
    └── multiqc_report.html         # Alignment QC
```

## Clean up

```bash
# Remove all generated files
pixi run snakemake clean

# Remove only work directory (keep results)
pixi run snakemake clean_work
```
