# `scripts/` — Standalone Utility Scripts

Self-contained helper scripts that complement `spider_silkome_module/`.
Each script is runnable on its own (no need to import as a Python package)
and is grouped by purpose below.

---

## Analysis Pipeline

Scripts that process HMMER / miniprot outputs and produce intermediate or
final analysis artifacts.

| Script | Purpose |
|---|---|
| `analyse_spidroins.py` | End-to-end HMMER-based spidroin analysis pipeline |
| `extract_terminal_domains.py` | Extract non-repetitive NTD/CTD regions from spidroin proteins |
| `miniprot_dedup.py` | Deduplicate overlapping miniprot hits, keeping best score |
| `miniprot_pair.py` | Pair N/C-terminal miniprot hits into candidate gene loci |
| `merge_gff.py` | Merge curated spidroin GFF with Augustus gene predictions |
| `spidroin_repeat_visual.py` | Visualize tandem repeat motifs in spidroin proteins |

### `analyse_spidroins.py`

Parses nhmmer table outputs across multiple species, filters by E-value /
coverage, joins NTD/CTD hits into full-length candidates, and emits summary
TSVs plus heatmap / bar-plot / length-distribution figures.

```bash
pixi run python scripts/analyse_spidroins.py --help
```

### `extract_terminal_domains.py`

Auto-detects the repeat-unit length in each input protein (autocorrelation),
trims the repeat region, and writes the non-repetitive terminal domain
sequences grouped by `{SpidroinType}_{NTD|CTD}.faa`. `TuSp` and `CySp` are
merged into `TuSp_CySp`.

```bash
pixi run python scripts/extract_terminal_domains.py seeds.faa -o output/ -s 19
```

### `miniprot_dedup.py`

Reads a miniprot mRNA GFF, filters by `--min-positive` (homology fraction),
collapses overlapping hits sharing `(seqid, gene_type, terminal)` keeping the
best-scoring one, and writes a deduplicated TSV.

```bash
pixi run python scripts/miniprot_dedup.py \
  --miniprot-gff miniprot_out.mRNA.gff \
  --output dedup_hits.tsv --min-positive 0.85
```

### `miniprot_pair.py`

Consumes the deduplicated TSV from `miniprot_dedup.py`, pairs N- and
C-terminal hits with same gene type that are 15–90 kb apart on the same
scaffold, and emits paired loci as BED. Auto-generates `.chrom.sizes` from a
reference FASTA when missing.

```bash
pixi run python scripts/miniprot_pair.py \
  --input dedup_hits.tsv \
  --out-paired paired_genes.bed --out-unpaired unpaired_hits.bed \
  --pad 1000 --min-distance 15000 --max-distance 90000 \
  --ref-fasta /path/to/genome.fa.gz
```

### `merge_gff.py`

Combines the curated spidroin GFF (manual / pipeline output) with Augustus's
re-prediction GFF, used during locus refinement.

### `spidroin_repeat_visual.py`

Detects common spidroin repeat patterns (e.g. poly-A, GPGGX, GA-rich) in
each input protein and renders a per-sequence diagram of where each motif
occurs along the chain.

```bash
pixi run python scripts/spidroin_repeat_visual.py -i input.fasta -o output_dir/
```

---

## Reports & Publishing

| Script | Purpose |
|---|---|
| `publish_report.py` | Deploy `reports/*.html` to Cloudflare Pages with an auto-generated landing index |

### `publish_report.py`

Idempotent Cloudflare Pages deployer. Auto-creates the project on first run,
stages a clean deploy directory, renders a card-style index page (titles
extracted from each report's `<title>` tag), then uploads. Subsequent runs
reuse the same project so the public URL stays stable for sharing in
Yuque/Feishu/etc.

First-time setup (one shot, opens a browser for OAuth):

```bash
npx wrangler@latest login
```

```bash
# Default: publish every reports/*.html
pixi run python scripts/publish_report.py

# Selective: pass specific files
pixi run python scripts/publish_report.py reports/agent_vs_manual_report.html

# Custom project (gives a different .pages.dev URL)
pixi run python scripts/publish_report.py --project-name=my-name reports/*.html
```

Output URL: `https://<project-name>.pages.dev/` — the root shows a card
index, and `/<report_name>.html` serves each individual report.
