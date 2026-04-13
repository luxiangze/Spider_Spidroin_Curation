"""
LLM-based spidroin type conflict resolver.

Called after resolve_spidroin_type() returns a cross-family conflict
(e.g. "MaSp/MiSp") to synthesise all available evidence and determine
the most likely type.

Evidence used:
  1. All raw nHMMER hits at the locus (before profile suppression) — so the
     LLM sees the full ranking of profile scores, not just the winning one.
  2. Miniprot vote summary — counts and best Positive score per type.
  3. Locus coordinates and strand for context.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import re

from dotenv import load_dotenv
from langchain_core.messages import HumanMessage, SystemMessage
from langchain_openai import ChatOpenAI
from loguru import logger

load_dotenv(Path(__file__).parent.parent.parent / ".env")

_DASHSCOPE_BASE_URL = "https://dashscope.aliyuncs.com/compatible-mode/v1"

_SYSTEM_PROMPT = """\
You are an expert in spider silk protein (spidroin) genomics.

You are given evidence for a spidroin gene locus where the automated classifier
produced a TYPE CONFLICT — the NTD (N-terminal domain) and CTD (C-terminal domain)
matched different spidroin profile families. Your task is to determine the most
likely true spidroin type.

## Key evidence to evaluate

1. **nHMMER profile hits** — HMM profile matches ranked by e-value (lower = better).
   The classifier used only the single best hit per terminus; you can see all of them.
   A very small e-value difference between profiles (< 1 order of magnitude) means
   the profiles are nearly equally good matches and cross-matching is likely.

2. **Miniprot reference protein vote** — count and best Positive score (0–1, higher = better)
   for each protein type that aligned to this locus. This reflects which reference
   proteins actually align well to the genomic sequence, independent of HMM profiles.
   Miniprot evidence is particularly valuable when nHMMER gives ambiguous results.

3. **Strand and genomic position** — for context only.

## Decision rules

- If the CTD clearly matches one type and the NTD's best profile differs only slightly
  in e-value from the same type's profile, prefer the CTD type supported by multiple
  miniprot hits.
- If miniprot votes strongly support one family (e.g. >60% of high-confidence hits),
  favour that family.
- If evidence is genuinely ambiguous (e.g. 50/50 miniprot split, NTD and CTD equally
  matched to different profiles), return the conflict string (e.g. "MaSp/MiSp") and
  set needs_review=true.
- Never invent a type not present in the evidence.

## Output format (JSON only, no markdown fences)

{
  "resolved_type": "<type string, e.g. MiSp, MaSp, AcSp — or the conflict string>",
  "needs_review": <bool>,
  "confidence": "high" | "medium" | "low",
  "reasoning": "<2-3 sentences>"
}
"""


def _summarise_miniprot_votes(mp_hits: list[dict]) -> str:
    """Tally miniprot hits by type extracted from ref_protein name."""
    from collections import defaultdict
    votes: dict[str, list[float]] = defaultdict(list)
    for h in mp_hits:
        ref = h.get("ref_protein", "")
        pos = h.get("positive", 0.0)
        # Extract type token: "6079_MiSp_NTD_319aa" → "MiSp"
        parts = ref.replace(".", "_").split("_")
        # Known type tokens (order matters — check longer names first)
        type_tokens = ["MaSp3B", "MaSp2B", "MaSp3", "MaSp2", "MaSp1", "MaSp",
                       "MiSp", "AcSp", "CySp", "PySp", "Flag", "Pflag",
                       "TuSp", "AgSp1", "AgSp2", "AgSp", "CrSp"]
        matched = None
        for tok in type_tokens:
            if tok in parts:
                matched = tok
                break
        if matched is None:
            matched = "Spidroin/unknown"
        votes[matched].append(pos)

    if not votes:
        return "  (no miniprot hits)"

    lines = []
    for typ, scores in sorted(votes.items(), key=lambda x: -len(x[1])):
        best = max(scores)
        lines.append(f"  {typ}: {len(scores)} hits, best_positive={best:.3f}")
    return "\n".join(lines)


def _build_evidence_text(
    locus: dict,
    raw_hmmer_hits: list[dict],
    mp_hits: list[dict],
    conflict_type: str,
) -> str:
    seqid = locus.get("seqid", "?")
    start = locus.get("start", 0)
    end   = locus.get("end", 0)
    strand = locus.get("strand", "?")

    # All nHMMER hits overlapping the locus, sorted by e-value
    locus_hits = [
        h for h in raw_hmmer_hits
        if h.get("seqid") == seqid
        and h.get("start", 0) <= end
        and h.get("end", 0)   >= start
    ]
    ntd_hits = sorted([h for h in locus_hits if h.get("terminus") == "NTD"],
                      key=lambda h: h.get("e_value", 1.0))
    ctd_hits = sorted([h for h in locus_hits if h.get("terminus") == "CTD"],
                      key=lambda h: h.get("e_value", 1.0))

    ntd_lines = "\n".join(
        f"  {h['profile']:20s}  e={h['e_value']:.2e}"
        for h in ntd_hits[:10]  # top 10 is enough
    ) or "  (none)"
    ctd_lines = "\n".join(
        f"  {h['profile']:20s}  e={h['e_value']:.2e}"
        for h in ctd_hits[:10]
    ) or "  (none)"

    mp_overlap = [
        h for h in mp_hits
        if h.get("seqid") == seqid
        and h.get("start", 0) <= end
        and h.get("end", 0)   >= start
    ]
    mp_summary = _summarise_miniprot_votes(mp_overlap)

    return (
        f"Locus: {seqid}:{start:,}-{end:,}  strand={strand}\n"
        f"Classifier conflict: {conflict_type}\n\n"
        f"nHMMER NTD hits (all profiles, ranked by e-value):\n{ntd_lines}\n\n"
        f"nHMMER CTD hits (all profiles, ranked by e-value):\n{ctd_lines}\n\n"
        f"Miniprot reference protein vote ({len(mp_overlap)} hits):\n{mp_summary}"
    )


def resolve_type_conflict_llm(
    locus: dict,
    conflict_type: str,
    raw_hmmer_hits: list[dict],
    mp_hits: list[dict],
) -> tuple[str, bool, str]:
    """
    Use LLM to resolve a spidroin type conflict.

    Args:
        locus:         The locus dict (needs seqid, start, end, strand).
        conflict_type: The conflict string from resolve_spidroin_type (e.g. "MaSp/MiSp").
        raw_hmmer_hits: All nHMMER rows for this species (unfiltered).
        mp_hits:       All miniprot rows for this species.

    Returns:
        (resolved_type, needs_review, reasoning)
        Falls back to the original conflict_type if LLM fails.
    """
    if not os.environ.get("DASHSCOPE_API_KEY"):
        return conflict_type, True, "llm_unavailable:no_api_key"

    evidence = _build_evidence_text(locus, raw_hmmer_hits, mp_hits, conflict_type)

    model = ChatOpenAI(
        model=os.environ.get("TYPE_RESOLVER_LLM_MODEL", "qwen3.6-plus"),
        openai_api_key=os.environ.get("DASHSCOPE_API_KEY", ""),
        openai_api_base=_DASHSCOPE_BASE_URL,
        max_tokens=512,
        temperature=0,
        extra_body={"enable_thinking": False},
        timeout=60,
    )

    messages = [
        SystemMessage(content=_SYSTEM_PROMPT),
        HumanMessage(content=evidence),
    ]

    try:
        response = model.invoke(messages)
        raw = response.content if isinstance(response.content, str) else str(response.content)
    except Exception as exc:
        logger.warning(
            f"[type_resolver] LLM call failed for {locus.get('seqid')}:{locus.get('start')}: {exc}"
        )
        return conflict_type, True, f"llm_error:{exc}"

    # Extract JSON
    match = re.search(r"\{[\s\S]*?\}", raw)
    if not match:
        logger.warning(f"[type_resolver] No JSON in LLM response: {raw[:200]}")
        return conflict_type, True, "llm_parse_error:no_json"

    try:
        result = json.loads(match.group())
    except json.JSONDecodeError as exc:
        logger.warning(f"[type_resolver] JSON parse error: {exc}")
        return conflict_type, True, f"llm_parse_error:{exc}"

    resolved = result.get("resolved_type", conflict_type).strip()
    needs_review = bool(result.get("needs_review", True))
    reasoning = result.get("reasoning", "").strip()
    conf = result.get("confidence", "medium")

    logger.debug(
        f"[type_resolver] {conflict_type} → {resolved} "
        f"(conf={conf}, needs_review={needs_review}): {reasoning[:80]}"
    )
    return resolved, needs_review, reasoning
