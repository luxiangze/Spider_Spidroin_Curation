# CLAUDE.md — Spider Silkome Project

## Project Goal

Build an automated spidroin identification agent to replace the current manual review workflow.

**Current workflow**: Run nHMMER domain search + Miniprot reference protein alignment → manual review in JBrowse2 → record into Feishu multi-dimensional table.

**Agent goal**: Automate the classification and completeness assessment steps, and output structured TSV results.

---

## Project Structure

```
spider_silkome/
├── agents/                         # Agent code (to be created)
│   └── spidroin_id_agent.py        # Main agent entry point
├── data/
│   ├── external/                   # External reference datasets (read-only)
│   ├── interim/                    # Intermediate processing data (read-only)
│   └── processed/                  # Processed results (read-only)
│       ├── miniprot_output/        # Miniprot GFF outputs (one .gff per species)
│       └── hmmer_search_20260211/
│           └── nhmmer_search_parsed/  # Parsed nhmmer results (one GFF per species: <species>/<species>.gff)
├── docs/
│   └── typing_rules.md             # Domain knowledge: spidroin typing rules (to be created)
├── models/                         # Trained model files
├── notebooks/                      # Jupyter notebooks for exploration
├── references/
│   └── 2025_Schoneberg_data/
│       └── hmmer_nucl_profile_trimmed/  # HMM profile files (.hmm) for NTD/CTD domains
├── reports/figures/                # Output figures
├── results/                        # Agent output results (to be created)
├── scripts/                        # Standalone analysis scripts
├── spider_silkome_module/          # Existing Python package
│   ├── config.py                   # Path constants (PROJ_ROOT, DATA_DIR, etc.)
│   ├── parse_nhmmer_results.py     # nhmmer tblout parser (reusable)
│   ├── run_nhmmer.py
│   ├── run_miniprot.py
│   └── utils/
├── tests/
├── workflow/                       # Snakemake/Nextflow workflows
├── pixi.toml                       # Package management (pixi)
└── pyproject.toml
```

**IMPORTANT**: `data/` directory is read-only. Never modify files under `data/`.

---

## Data Locations

| Data Type | Path |
|-----------|------|
| Miniprot GFF outputs | `data/processed/miniprot_output/*.gff` |
| nHMMER parsed results | `data/processed/hmmer_search_20260211/nhmmer_search_parsed/<species>/<species>.gff` |
| HMM profile files | `references/2025_Schoneberg_data/hmmer_nucl_profile_trimmed/*.hmm` |
| Reference spidroin proteins | `data/external/spider-silkome-database.v1.prot.fixed.renamed.fasta` |
| NTD/CTD seed sequences | `data/interim/spidroin_proteins/` |

---

## Agent Architecture

Agent code lives in `agents/`. Results are written to `results/`.

### Agent Pipeline Steps

1. **Parse nHMMER tblout** — extract domain hits: profile name, e-value, coordinates (env_from/env_to), strand
2. **Parse Miniprot GFF** — extract alignment features: reference protein, query scaffold, coordinates, identity
3. **Classify spidroin type** — based on NTD/CTD domain hit combination, determine type (MaSp/MiSp/PySp/AgSp/Flag/TuSp/AcSp/CySp/CrSp, etc.)
4. **Assess completeness** — categorize as: `Full_length` / `N-terminal` / `C-terminal` / `Repeat_only`
5. **Output structured TSV** — high-confidence entries written directly; low-confidence entries flagged as `needs_review`
6. **Feishu API interface** — reserved but not implemented; design as a pluggable output module

### Confidence Thresholds (default, adjustable)

- nHMMER: e-value < 1e-10, HMM coverage ≥ 90%, score > bias
- Miniprot: identity ≥ 70% (TODO: confirm threshold from domain expert)

### TSV Output Schema (`results/spidroin_candidates.tsv`)

| Column | Description |
|--------|-------------|
| species | Species name |
| scaffold | Genome scaffold/contig ID |
| start | Start coordinate |
| end | End coordinate |
| strand | Strand (+/-) |
| spidroin_type | Predicted type (MaSp1, MiSp, etc.) |
| completeness | Full_length / N-terminal / C-terminal / Repeat_only |
| ntd_hit | NTD domain profile name (or empty) |
| ctd_hit | CTD domain profile name (or empty) |
| ntd_evalue | NTD hit e-value |
| ctd_evalue | CTD hit e-value |
| miniprot_ref | Best matching reference protein |
| miniprot_identity | Miniprot alignment identity |
| confidence | high / low |
| review_flag | needs_review / auto_accepted |
| notes | Free-text notes for manual review |

---

## Existing Reusable Code

Before writing new parsers, check the existing module:

- `spider_silkome_module/parse_nhmmer_results.py` — `parse_tbl_file()`, `filter_hits()`, `process_species()` are already implemented and tested
- `spider_silkome_module/config.py` — defines `PROJ_ROOT`, `DATA_DIR`, `PROCESSED_DATA_DIR`, etc.
- `spider_silkome_module/spidroin_analysis_aa.py` — amino acid level analysis (check before reimplementing)

**Reuse existing parsers rather than rewriting them.**

---

## HMM Profile Naming Convention

Profile filenames follow the pattern `{SpidroinType}_{Domain}.hmm`, e.g.:
- `MaSp1_NTD.hmm`, `MaSp1_CTD.hmm`
- `Flag_NTD.hmm`, `Flag_CTD.hmm`
- `Spidroin_NTD.hmm` (generic spidroin NTD, fallback)

Available types inferred from `data/interim/spidroin_proteins/`:
AcSp, AgSp1, AgSp2, CrSp, CySp, Flag, MaSp, MaSp1, MaSp2, MaSp2B, MaSp3, MaSp3B, MiSp, Pflag, PySp, Spidroin (generic)

---

## Coding Conventions

- **Language**: Python 3.10+ (project uses 3.12, see `pixi.toml`)
- **Package management**: `pixi` — do not use pip/conda directly
- **Data frames**: use `polars` (already a project dependency), not pandas, for new agent code
- **Logging**: use `loguru` (already a project dependency)
- **CLI**: use `typer` for command-line interfaces
- **Path handling**: use `pathlib.Path`; import path constants from `spider_silkome_module.config`
- **Server**: Ubuntu, paths are case-sensitive
- **Code comments**: English only
- **Formatting**: `ruff` is configured (run `pixi run ruff check` to lint)

---

## Domain Knowledge

Typing rules are maintained separately in `docs/typing_rules.md`.

**Do not hardcode typing logic in agent code.** Load rules from `docs/typing_rules.md` or a structured config (YAML/JSON) so rules can be updated without code changes.

---

## Development Rules

1. **`data/` is read-only** — never write to `data/`
2. **`results/` for outputs** — all agent-generated files go here
3. **`agents/` for agent code** — keep agent logic separate from the existing `spider_silkome_module/`
4. **Reuse `spider_silkome_module`** — import existing parsers and utilities; don't duplicate
5. **Write TODO comments** when logic is uncertain — never silently assume domain behavior
6. **Validate inputs** — check file existence and format at system boundaries; log counts at key steps
7. **No Feishu API yet** — design the output module with a clear interface so Feishu integration can be plugged in later

---

## Running the Agent

```bash
# Activate environment
pixi shell

# Run agent (example, exact CLI TBD)
python agents/spidroin_id_agent.py \
  --nhmmer-dir data/processed/hmmer_search_20260211/nhmmer_search_parsed \
  --miniprot-dir data/processed/miniprot_output \
  --output results/spidroin_candidates.tsv
```

---

## Feishu Integration (Reserved Interface)

The agent should expose a function:

```python
def export_to_feishu(df: pl.DataFrame, table_id: str) -> None:
    # TODO: implement Feishu multi-dimensional table API write
    raise NotImplementedError("Feishu export not yet implemented")
```

Do not implement this until explicitly requested.
