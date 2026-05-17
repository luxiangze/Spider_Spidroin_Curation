"""
Shared fixtures for the species 001 integration test suite.

Runs _process_species() once per session (module scope) and writes results to
tests/test_species_001/output/ for visual inspection.
"""

from __future__ import annotations

import csv
from pathlib import Path

import polars as pl
import pytest
import yaml

from agents.typing_agent import _process_species
from spider_silkome_module.config import PROCESSED_DATA_DIR, PROJ_ROOT

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

HERE = Path(__file__).parent
SPECIES_ID = "001.Allagelena_difficilis"
NHMMER_DIR = PROCESSED_DATA_DIR / "hmmer_search_20260211" / "nhmmer_search_parsed"
MINIPROT_DIR = PROCESSED_DATA_DIR / "miniprot_output"
RULES_FILE = PROJ_ROOT / "docs" / "typing_rules.yaml"
OUTPUT_DIR = HERE / "output"


# ---------------------------------------------------------------------------
# CLI option
# ---------------------------------------------------------------------------

def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption(
        "--run-llm",
        action="store_true",
        default=False,
        help="Run @pytest.mark.llm tests that require ANTHROPIC_API_KEY (slow, ~2 min)",
    )


def pytest_configure(config: pytest.Config) -> None:
    config.addinivalue_line(
        "markers",
        "llm: tests that require ANTHROPIC_API_KEY and the --run-llm flag",
    )


def pytest_collection_modifyitems(
    config: pytest.Config, items: list[pytest.Item]
) -> None:
    if not config.getoption("--run-llm"):
        skip_llm = pytest.mark.skip(reason="Pass --run-llm to enable LLM tests")
        for item in items:
            if "llm" in item.keywords:
                item.add_marker(skip_llm)


# ---------------------------------------------------------------------------
# Helper: build _process_species kwargs from typing_rules.yaml
# (mirrors the parsing logic in typing_agent.main(), lines 1362–1382)
# ---------------------------------------------------------------------------

def _build_kwargs(rules: dict, use_llm: bool = False) -> dict:
    thresholds = rules.get("thresholds", {})
    nhmmer_cfg = thresholds.get("nhmmer", {})
    miniprot_cfg = thresholds.get("miniprot", {})
    type_to_family: dict[str, str] = {
        member: family
        for family, members in rules.get("type_families", {}).items()
        for member in members
        if member != family
    }
    return dict(
        miniprot_dir=MINIPROT_DIR,
        e_value_threshold=nhmmer_cfg.get("e_value", 1e-4),
        e_value_high=nhmmer_cfg.get("e_value_high", 1e-10),
        min_span_bp=int(thresholds.get("min_span_bp", 1_000)),
        max_span_bp=int(thresholds.get("max_span_bp", 150_000)),
        positive_high=miniprot_cfg.get("positive_high", 0.6),
        type_mapping=rules.get("type_mapping", {}),
        specificity_pairs={(a, b) for a, b in rules.get("type_specificity", [])},
        type_to_family=type_to_family,
        anno_new=None,
        anno_old=None,
        bw_cache=PROJ_ROOT / "cache" / "bw",
        bw_index=None,
        skip_anno=True,
        skip_codon=True,
        skip_rna=True,
        use_llm=use_llm,
    )


# ---------------------------------------------------------------------------
# Helper: write loci list to a TSV for visual inspection
# ---------------------------------------------------------------------------

def _write_tsv(loci_list: list[dict], path: Path) -> None:
    rows = []
    for i, locus in enumerate(
        sorted(loci_list, key=lambda x: (x.get("seqid", ""), x.get("start", 0))), 1
    ):
        start = locus.get("start") or 0
        end = locus.get("end") or 0
        rows.append(
            {
                "Spidroin_ID": f"Aldi_spid{i:05d}",
                "Chr": locus.get("seqid", ""),
                "Start": start,
                "End": end,
                "Stran": locus.get("strand", ""),
                "Spidroin_type": locus.get("spidroin_type", ""),
                "Full_length": locus.get("completeness") == "Full_length",
                "Length": end - start if end and start else 0,
                "Hint_type": locus.get("completeness", ""),
                "ntd_profile": locus.get("ntd_profile") or "",
                "ntd_e_value": locus.get("ntd_e_value") or "",
                "ctd_profile": locus.get("ctd_profile") or "",
                "ctd_e_value": locus.get("ctd_e_value") or "",
                "miniprot_ref": locus.get("miniprot_ref") or "",
                "miniprot_positive": locus.get("miniprot_positive") or "",
                "confidence": locus.get("confidence") or "",
                "needs_review": locus.get("needs_review") or "",
                "llm_confidence": locus.get("llm_confidence") or "",
                "llm_reasoning": locus.get("llm_reasoning") or "",
            }
        )
    df = pl.DataFrame(rows, infer_schema_length=None)
    df.write_csv(path, separator="\t")


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def rules() -> dict:
    with open(RULES_FILE) as fh:
        return yaml.safe_load(fh)


@pytest.fixture(scope="module")
def loci(rules: dict) -> list[dict]:
    """Run _process_species without LLM (fast, ~20–30 s). Writes 001.no_llm.tsv."""
    OUTPUT_DIR.mkdir(exist_ok=True)
    result = _process_species(
        SPECIES_ID,
        NHMMER_DIR / SPECIES_ID,
        **_build_kwargs(rules, use_llm=False),
    )
    _write_tsv(result, OUTPUT_DIR / "001.no_llm.tsv")
    return result


@pytest.fixture(scope="module")
def loci_llm(rules: dict) -> list[dict]:
    """Run _process_species with LLM pairing. Writes 001.llm.tsv. Requires --run-llm."""
    OUTPUT_DIR.mkdir(exist_ok=True)
    result = _process_species(
        SPECIES_ID,
        NHMMER_DIR / SPECIES_ID,
        **_build_kwargs(rules, use_llm=True),
    )
    _write_tsv(result, OUTPUT_DIR / "001.llm.tsv")
    return result


@pytest.fixture(scope="module")
def ground_truth() -> list[dict]:
    """Load and parse the human-curated loci from ground_truth.csv."""
    rows: list[dict] = []
    with open(HERE / "ground_truth.csv", newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            chr_val = row.get("Chr", "").strip()
            if not chr_val:
                continue
            score_raw = row.get("Scoring", "").strip()
            rows.append(
                {
                    "seqid": chr_val,
                    "start": int(float(row["Start"])),
                    "end": int(float(row["End"])),
                    "strand": row.get("Stran", "").strip(),
                    "spidroin_type": row.get("Spidroin_type", "").strip(),
                    "full_length": row.get("Full_length", "").strip().lower() == "true",
                    "score": int(score_raw) if score_raw else 0,
                    "spidroin_id": row.get("Spidroin_ID", "").strip(),
                }
            )
    return rows
