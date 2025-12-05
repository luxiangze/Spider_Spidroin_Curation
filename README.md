# Spider Spidroin Curation

<a target="_blank" href="https://cookiecutter-data-science.drivendata.org/">
    <img src="https://img.shields.io/badge/CCDS-Project%20template-328F97?logo=cookiecutter" />
</a>

**[中文版](README_CN.md)**

Automated annotation and expression analysis pipeline for spider silk proteins (Spidroins).

## Quick Start

```bash
# Install dependencies
pixi install

# Activate environment
pixi shell

# Run Spidroin annotation workflow (Jupyter Notebook, configure parameters first)
jupyter lab notebooks/Automated_spidroin_annotation.ipynb

# Run BGI short-read RNA-seq alignment
cd workflow/BGI_RNA_align
pixi run snakemake --cores 80 --resources mem_gb=240

# Run ONT long-read transcriptome assembly (requires BGI workflow to complete first)
cd workflow/ONT_RNA_align
pixi run snakemake --cores 80 --resources mem_gb=240
```

## Features

### Spidroin Automated Annotation
- ✅ nhmmer search for Spidroin N/C terminal sequences
- ✅ Multi-species batch analysis (10 spider genomes)
- ✅ Augustus gene structure prediction
- ✅ Protein sequence extraction

### RNA-seq Analysis
- ✅ BGI short-read STAR alignment (`workflow/BGI_RNA_align`)
- ✅ ONT long-read FLAIR transcriptome assembly (`workflow/ONT_RNA_align`)
- ✅ FastQC + fastp quality control
- ✅ BigWig visualization files
- ✅ MultiQC summary reports

## Project Structure

```
spider_silkome/
├── data/
│   ├── raw/                       # Raw data
│   │   ├── spider_genome/         # Spider genomes (fa + gff)
│   │   ├── BGI_RNA_10samples/     # BGI short-read RNA-seq
│   │   └── ONT_DRS_10samples/     # ONT long-read RNA-seq
│   ├── interim/                   # Intermediate results
│   └── processed/                 # Final outputs
│
├── workflow/
│   ├── BGI_RNA_align/             # BGI short-read STAR alignment
│   │   ├── Snakefile
│   │   ├── config/
│   │   └── rules/
│   └── ONT_RNA_align/             # ONT long-read FLAIR assembly
│       ├── Snakefile
│       ├── config/
│       └── rules/
│
├── spider_silkome_module/         # Python module
│   ├── config.py                  # Path configuration
│   ├── features.py                # Utility functions
│   └── plots.py                   # Visualization
│
├── notebooks/
│   └── Automated_spidroin_annotation.ipynb  # Spidroin annotation workflow
├── scripts/
│   └── analyse_spidroins.py       # Spidroin analysis script
└── Makefile
```

## Dependencies

This project uses [pixi](https://pixi.sh/) for dependency management. Key tools include:

- **Spidroin annotation**: nhmmer (HMMER), Augustus, BioPython
- **Short-read alignment**: STAR, samtools
- **Long-read alignment**: minimap2, FLAIR
- **Quality control**: FastQC, fastp, MultiQC
- **Visualization**: deeptools (bamCoverage)
- **Workflow**: Snakemake, Jupyter

## License

- **License**: See [LICENSE](LICENSE) file
- **Maintainer**: Yongkang Guo
- **Template**: Based on [Cookiecutter Data Science](https://cookiecutter-data-science.drivendata.org/)

