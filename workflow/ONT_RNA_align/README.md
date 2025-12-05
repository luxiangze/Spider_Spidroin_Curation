# ONT RNA Alignment Workflow

Oxford Nanopore long-read RNA-seq alignment pipeline using FLAIR for transcriptome assembly.

## Overview

This workflow processes ONT Direct RNA Sequencing (DRS) data for 10 spider species.

### Pipeline Steps

1. **minimap2** - Align ONT reads to genome
2. **flair transcriptome** - Transcriptome assembly (combines correct + collapse)
   - Uses BGI short-read STAR SJ.out.tab as splice junction reference
   - Uses `--stringent` mode to filter incomplete transcripts

## Data Structure

```
data/raw/
├── spider_genome/                    # Genome files
│   ├── Araneus_ventricosus/
│   │   ├── Araneus_ventricosus.fa
│   │   └── Araneus_ventricosus.gff
│   └── ...
└── ONT_DRS_10samples/                # ONT RNA-seq data
    ├── Araneus_ventricosus/
    │   └── pass.fq.gz
    └── ...
```

## Usage

### 1. Dry run (validate workflow)

```bash
cd workflow/ONT_RNA_align
pixi run snakemake -n
```

### 2. Run alignment only (can run before BGI workflow completes)

```bash
pixi run snakemake all_aligned --cores 80 --resources mem_gb=240
```

### 3. Run full pipeline (requires BGI RNA-seq to complete first)

```bash
# Ensure BGI_RNA_align has completed to get SJ.out.tab files
pixi run snakemake --cores 80 --resources mem_gb=240
```

### 4. Process single sample

```bash
pixi run snakemake results/.done.Araneus_ventricosus --cores 80 --resources mem_gb=240
```

## Resource Configuration

Default settings for 80-core / 252GB server:

| Task | Threads | Memory |
|------|---------|--------|
| minimap2 | 16 | 50GB |
| flair_transcriptome | 16 | 100GB |
| bigwig | 8 | - |

Edit `config/config.yaml` to adjust parameters.

## Output

```
results/                         # Final results
├── isoforms/
│   ├── {sample}.isoforms.gtf    # Transcript annotation (GTF)
│   ├── {sample}.isoforms.fa     # Transcript sequences (FASTA)
│   └── {sample}.isoforms.bed    # Transcript coordinates (BED)
├── bam/
│   ├── {sample}.bam             # Raw alignment
│   ├── {sample}.sorted.bam      # Sorted BAM (for visualization)
│   └── {sample}.sorted.bam.csi  # CSI index
├── bigwig/
│   └── {sample}.bw              # Coverage track
├── stats/
│   └── {sample}.flagstat.txt    # Alignment statistics
└── multiqc/
    └── multiqc_report.html      # Summary report

work/                            # Intermediate files (gitignored)
├── aligned/                     # minimap2 alignments
└── transcriptome/               # FLAIR transcriptome output
```

## Important Notes for Spidroin Analysis

⚠️ **FLAIR collapse may underestimate repetitive exon lengths!**

For spidroin genes with large repetitive regions:

1. **Don't trust collapse output blindly** - Assembled transcripts may be shorter than reality

2. **Verify in Geneious**:
   - Import `{sample}.isoforms.gtf` as reference
   - Import `{sample}.sorted.bam` to view raw reads
   - Zoom into repetitive regions, check actual read lengths

3. **`--stringent` mode** (enabled by default):
   - Requires reads to cover 80% of transcript
   - Requires first and last exon support
   - Helps filter fragmented transcripts

4. **Manual verification workflow**:
   ```
   1. Load genome FASTA
   2. Load BAM (raw read alignments)
   3. Compare read lengths vs FLAIR consensus in repetitive regions
   ```

## Clean up

```bash
# Remove all generated files
pixi run snakemake clean

# Remove only work directory (keep results)
pixi run snakemake clean_work
```
