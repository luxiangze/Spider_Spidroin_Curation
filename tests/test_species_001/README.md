# Species 001 Integration Tests — *Allagelena difficilis*

## Overview

This test suite validates the automated spidroin typing agent against manually
curated results for species 001 (*Allagelena difficilis*, abbreviated `Aldi`).

All tests run `_process_species()` directly (bypassing the CLI) with annotation,
codon-check, and RNA-seq evidence disabled (`skip_anno`, `skip_codon`, `skip_rna`)
to avoid external path dependencies and keep runtime short (~30 s).

---

## Ground Truth

`ground_truth.csv` — 23 spidroin loci manually reviewed gene-by-gene in JBrowse2
on 2026-03-22. Each entry has a **confidence score** (1–5):

| Score | Meaning |
|-------|---------|
| 5 | Unambiguous — high-quality NTD and CTD hits, strong miniprot support |
| 4 | Confident — clear NTD+CTD pair |
| 3 | Likely correct — weaker evidence |
| 2 | Partial — only a CTD detected (~330 bp) |
| 1 | Very uncertain |

Evidence screenshots are stored in `attachments/`.

**Summary:**

| Type | Completeness | Count | Score range |
|------|-------------|-------|-------------|
| AcSp | Full_length | 5 | 3–4 |
| AcSp | C-terminal  | 7 | 1–2 |
| PySp | Full_length | 3 | 5 |
| CySp | Full_length | 2 | 4 |
| MiSp | Full_length | 5 | 5 |

---

## Key Test Case: chr07 Ambiguous Window

Region `chr07:71,243,836–71,465,633` contains **1 NTD + 6 CTD** hits after
suppression — the classic "multiple adjacent CTDs" scenario shown in `agents/samples/sample7.png`.

Human review decision:
- **spid00004–00008**: 5 × C-terminal (each ~330 bp, score 2, not full-length)
- **spid00009**: Full_length AcSp at `71,442,768–71,465,633` (score 4)
  — NTD at 71,465,162 paired with the **nearest CTD** at 71,442,774 (22,860 bp span),
  supported by 11 non-overlapping miniprot alignments spanning the interval.

This window specifically tests whether the pairing algorithm avoids incorrectly
pairing the NTD with more distant CTDs.

---

## How to Run

```bash
# Fast (no LLM calls, ~30 s)
pixi run pytest tests/test_species_001/ -v

# With LLM-based pairing (requires ANTHROPIC_API_KEY in .env, ~2 min)
pixi run pytest tests/test_species_001/ -v --run-llm
```

---

## Output Files

After running, `output/` contains:

| File | Description |
|------|-------------|
| `001.no_llm.tsv` | Results using greedy pairing (no LLM) |
| `001.llm.tsv` | Results using LangGraph + Claude pairing (--run-llm only) |

Columns include standard Feishu-aligned fields plus `llm_confidence` and `llm_reasoning`
for inspecting LLM decisions.

---

## Test Descriptions

| Test | Description |
|------|-------------|
| `test_pipeline_runs` | Smoke test — agent completes without errors, all loci have valid structure |
| `test_high_confidence_full_length` | All 13 score-≥4 Full_length loci from human review must be found |
| `test_ambiguous_window_chr07` | Correct NTD-CTD pairing in the complex window; ≥4 unpaired CTDs remain |
| `test_no_false_full_length_in_ambiguous_window` | The NTD must not be paired with a distant wrong CTD |
| `test_llm_pairing` *(--run-llm)* | LLM path gives `high` confidence with spanning miniprot evidence |
| `test_output_tsv_written` | TSV output file exists and has ≥15 data rows |
