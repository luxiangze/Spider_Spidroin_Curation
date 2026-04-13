"""
Integration tests for the spidroin typing agent on species 001 (Allagelena difficilis).

Ground truth: 23 loci manually reviewed in JBrowse2 (ground_truth.csv).
Tests run _process_species() with skip_anno/codon/rna=True to avoid external dependencies.

Run:
    pixi run pytest tests/test_species_001/ -v              # fast (~30s)
    pixi run pytest tests/test_species_001/ -v --run-llm    # with LLM (~2 min)
"""

from __future__ import annotations

from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _type_compatible(agent_type: str | None, human_type: str) -> bool:
    """
    Return True if the agent's spidroin_type is compatible with the human-reviewed type.
    Handles the MiSp vs MaSp/MiSp conflict (agent may classify chr15 MiSp as MaSp/MiSp
    due to HMM cross-matching; both are acceptable for the same biological locus).
    """
    if not agent_type or not human_type:
        return False
    return human_type in agent_type or agent_type in human_type


def _find_matching_locus(
    loci: list[dict],
    human_row: dict,
    coord_tol: int = 1_000,
) -> dict | None:
    """
    Find the agent locus closest to a human-annotated entry.
    Matching criteria: same seqid, start within coord_tol, end within coord_tol.
    Type-compatible matches are preferred over type-mismatched ones.
    Returns None if no candidate is found.
    """
    seqid = human_row["seqid"]
    h_start = human_row["start"]
    h_end = human_row["end"]
    h_type = human_row["spidroin_type"]

    candidates = [
        l for l in loci
        if l.get("seqid") == seqid
        and abs((l.get("start") or 0) - h_start) < coord_tol
        and abs((l.get("end") or 0) - h_end) < coord_tol
    ]
    if not candidates:
        return None

    typed = [c for c in candidates if _type_compatible(c.get("spidroin_type"), h_type)]
    return typed[0] if typed else candidates[0]


# ---------------------------------------------------------------------------
# Smoke test
# ---------------------------------------------------------------------------

def test_pipeline_runs(loci: list[dict]) -> None:
    """Agent completes without errors; all loci have valid structure."""
    assert len(loci) > 0, "No loci returned"
    for locus in loci:
        start = locus.get("start")
        end = locus.get("end")
        assert start is not None and end is not None, f"Missing coordinates: {locus}"
        assert start < end, f"Invalid span (start >= end): {locus}"
        assert locus.get("confidence") in {"high", "medium", "low"}, (
            f"Unexpected confidence value: {locus.get('confidence')}"
        )


# ---------------------------------------------------------------------------
# Core regression: all high-confidence human loci must be found
# ---------------------------------------------------------------------------

def test_high_confidence_full_length(
    loci: list[dict], ground_truth: list[dict]
) -> None:
    """
    Every score->=4 Full_length entry from the human review must be detected
    as a Full_length locus (completeness=='Full_length').
    Type matching allows MiSp / MaSp/MiSp equivalence.
    """
    targets = [r for r in ground_truth if r["score"] >= 4 and r["full_length"]]
    assert len(targets) == 13, f"Expected 13 high-conf Full_length entries, got {len(targets)}"

    failures: list[str] = []
    for hr in targets:
        match = _find_matching_locus(loci, hr)
        if match is None:
            failures.append(
                f"  NOT FOUND: {hr['spidroin_id']} {hr['seqid']}:{hr['start']}-{hr['end']} "
                f"({hr['spidroin_type']})"
            )
        elif match.get("completeness") != "Full_length":
            failures.append(
                f"  WRONG COMPLETENESS: {hr['spidroin_id']} {hr['seqid']}:{hr['start']}-{hr['end']} "
                f"→ agent completeness={match.get('completeness')}"
            )

    assert not failures, (
        f"{len(failures)}/13 high-confidence Full_length loci failed:\n"
        + "\n".join(failures)
    )


# ---------------------------------------------------------------------------
# Key case: chr07 ambiguous window (1 NTD + 6 CTD)
# ---------------------------------------------------------------------------

def test_ambiguous_window_chr07(loci: list[dict]) -> None:
    """
    The chr07:71,243,836–71,465,633 window has 1 NTD + 6 CTD after suppression.
    Expected outcome (Aldi_spid00009 from human review):
      - Exactly 1 Full_length AcSp: NTD at ~71,465,162, CTD at ~71,442,774 (22,860 bp)
      - At least 4 remaining CTDs are CTD-only candidates
    """
    # --- Full_length check ---
    full_length_hits = [
        l for l in loci
        if l.get("seqid") == "chr07"
        and l.get("has_ntd") and l.get("has_ctd")
        and abs((l.get("ntd_start") or 0) - 71_465_162) < 500
        and "AcSp" in (l.get("spidroin_type") or "")
    ]
    assert len(full_length_hits) == 1, (
        f"Expected 1 Full_length AcSp in ambiguous window, got {len(full_length_hits)}: "
        f"{[(l.get('start'), l.get('end'), l.get('spidroin_type')) for l in full_length_hits]}"
    )
    locus = full_length_hits[0]

    assert locus.get("confidence") == "high", (
        f"Expected high confidence for spid00009 equivalent, got {locus.get('confidence')}"
    )
    assert locus.get("completeness") == "Full_length"

    ctd_start = locus.get("ctd_start") or 0
    assert abs(ctd_start - 71_442_774) < 500, (
        f"NTD paired to wrong CTD: ctd_start={ctd_start} (expected ~71,442,774). "
        "The NTD must pair with the nearest CTD, not a more distant one."
    )

    # --- CTD-only check ---
    ctd_only = [
        l for l in loci
        if l.get("seqid") == "chr07"
        and not l.get("has_ntd")
        and l.get("has_ctd")
        and 71_240_000 < (l.get("ctd_start") or 0) < 71_460_000
    ]
    assert len(ctd_only) >= 4, (
        f"Expected ≥4 CTD-only loci in chr07 window, got {len(ctd_only)}: "
        f"{[(l.get('ctd_start'), l.get('ctd_e_value')) for l in ctd_only]}"
    )


def test_no_false_full_length_in_ambiguous_window(loci: list[dict]) -> None:
    """
    The NTD at ~71,465,162 must NOT be paired with CTD[1] at ~71,288,648
    (176 kb away, well beyond the 150 kb max span). This would be a false full-length.
    """
    false_pairs = [
        l for l in loci
        if l.get("seqid") == "chr07"
        and l.get("has_ntd") and l.get("has_ctd")
        and abs((l.get("ntd_start") or 0) - 71_465_162) < 500
        and abs((l.get("ctd_start") or 0) - 71_288_648) < 500
    ]
    assert len(false_pairs) == 0, (
        f"Found incorrect NTD–CTD pair (distance >150 kb): {false_pairs}"
    )


# ---------------------------------------------------------------------------
# LLM-specific assertions (--run-llm required)
# ---------------------------------------------------------------------------

@pytest.mark.llm
def test_llm_pairing(loci_llm: list[dict]) -> None:
    """
    With LLM pairing enabled, the chr07 ambiguous window must still produce the
    correct Full_length AcSp, with llm_confidence='high' and spanning-miniprot evidence.
    """
    fl = [
        l for l in loci_llm
        if l.get("seqid") == "chr07"
        and l.get("has_ntd") and l.get("has_ctd")
        and abs((l.get("ntd_start") or 0) - 71_465_162) < 500
        and "AcSp" in (l.get("spidroin_type") or "")
    ]
    assert len(fl) == 1, f"Expected 1 Full_length AcSp in LLM run, got {len(fl)}"
    locus = fl[0]

    assert locus.get("llm_confidence") == "high", (
        f"LLM confidence should be 'high', got {locus.get('llm_confidence')}"
    )
    assert locus.get("llm_reasoning"), "LLM reasoning should not be empty"

    evidence = " ".join(locus.get("llm_evidence") or [])
    assert "spanning_miniprot" in evidence, (
        f"Expected 'spanning_miniprot' in LLM evidence tags, got: {evidence!r}"
    )


# ---------------------------------------------------------------------------
# Output file check
# ---------------------------------------------------------------------------

def test_output_tsv_written(loci: list[dict]) -> None:
    """TSV output file must exist and contain at least 15 data rows."""
    _ = loci  # ensure loci fixture (and thus TSV write) has run
    tsv = Path(__file__).parent / "output" / "001.no_llm.tsv"
    assert tsv.exists(), f"Output TSV not written: {tsv}"
    lines = tsv.read_text(encoding="utf-8").splitlines()
    data_rows = len(lines) - 1  # subtract header
    assert data_rows >= 15, f"Output TSV has only {data_rows} data rows (expected ≥15)"
