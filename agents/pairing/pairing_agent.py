"""
LangGraph-based NTD/CTD pairing agent.

Graph:
  build_evidence → call_claude → parse_response
                                      ↓ (parse ok)   ↓ (parse failed)
                                  emit_loci ←── fallback_greedy

Entry point: run_pairing_agent()
"""

from __future__ import annotations

import json
import os
import re
from pathlib import Path
from typing import TypedDict

from dotenv import load_dotenv
from langchain_openai import ChatOpenAI

# Load .env from the project root (two levels up from this file: agents/pairing/ → project root)
load_dotenv(Path(__file__).parent.parent.parent / ".env")

_DASHSCOPE_BASE_URL = "https://dashscope.aliyuncs.com/compatible-mode/v1"
from langchain_core.messages import HumanMessage, SystemMessage
from langgraph.graph import END, START, StateGraph
from loguru import logger

from agents.pairing.evidence_builder import build_evidence_text
from agents.pairing.prompts import FEW_SHOT_EXAMPLE, SYSTEM_PROMPT

# ---------------------------------------------------------------------------
# State
# ---------------------------------------------------------------------------

class PairingState(TypedDict):
    # Inputs
    species: str
    seqid: str
    ntd_hits: list[dict]
    ctd_hits: list[dict]
    miniprot_hits: list[dict]
    min_span_bp: int
    max_span_bp: int

    # Intermediate
    evidence_text: str
    llm_raw_response: str

    # Outputs
    final_loci: list[dict]
    used_fallback: bool
    parse_error: str


# ---------------------------------------------------------------------------
# Import _make_locus lazily to avoid circular imports
# ---------------------------------------------------------------------------

def _make_locus_shim(species: str, seqid: str, ntd: dict | None, ctd: dict | None) -> dict:
    from agents.typing_agent import _make_locus
    return _make_locus(species, seqid, ntd, ctd)


# ---------------------------------------------------------------------------
# Greedy fallback (mirrors typing_agent logic, self-contained)
# ---------------------------------------------------------------------------

def _greedy_pair_local(
    species: str,
    seqid: str,
    ntd_hits: list[dict],
    ctd_hits: list[dict],
    min_span_bp: int,
    max_span_bp: int,
) -> list[dict]:
    """Fallback greedy pairing, identical to the original assemble_loci inner loop."""
    ntd_used = [False] * len(ntd_hits)
    ctd_used = [False] * len(ctd_hits)

    pairs: list[tuple[float, int, int]] = []
    for i, ntd in enumerate(ntd_hits):
        ntd_mid = (ntd["start"] + ntd["end"]) / 2
        for j, ctd in enumerate(ctd_hits):
            ctd_mid = (ctd["start"] + ctd["end"]) / 2
            dist = abs(ntd_mid - ctd_mid)
            if dist <= max_span_bp:
                pairs.append((dist, i, j))
    pairs.sort()

    loci: list[dict] = []
    for _dist, i, j in pairs:
        if ntd_used[i] or ctd_used[j]:
            continue
        ntd = ntd_hits[i]
        ctd = ctd_hits[j]
        span = max(ntd["end"], ctd["end"]) - min(ntd["start"], ctd["start"]) + 1
        if span < min_span_bp:
            continue
        loci.append(_make_locus_shim(species, seqid, ntd=ntd, ctd=ctd))
        ntd_used[i] = True
        ctd_used[j] = True

    for i, ntd in enumerate(ntd_hits):
        if not ntd_used[i]:
            loci.append(_make_locus_shim(species, seqid, ntd=ntd, ctd=None))
    for j, ctd in enumerate(ctd_hits):
        if not ctd_used[j]:
            loci.append(_make_locus_shim(species, seqid, ntd=None, ctd=ctd))

    return loci


# ---------------------------------------------------------------------------
# Graph nodes
# ---------------------------------------------------------------------------

def node_build_evidence(state: PairingState) -> dict:
    text = build_evidence_text(
        species=state["species"],
        seqid=state["seqid"],
        ntd_hits=state["ntd_hits"],
        ctd_hits=state["ctd_hits"],
        mp_hits=state["miniprot_hits"],
        min_span_bp=state["min_span_bp"],
        max_span_bp=state["max_span_bp"],
    )
    return {"evidence_text": text}


def node_call_claude(state: PairingState) -> dict:
    model = ChatOpenAI(
        model=os.environ.get("PAIRING_LLM_MODEL", "qwen3.6-plus"),
        openai_api_key=os.environ.get("DASHSCOPE_API_KEY", ""),
        openai_api_base=_DASHSCOPE_BASE_URL,
        max_tokens=2048,
        temperature=0,
        # Disable chain-of-thought for qwen3 thinking models to avoid slow responses.
        # extra_body is forwarded as-is to the API request body by langchain-openai.
        extra_body={"enable_thinking": False},
        timeout=120,
    )
    user_content = (
        f"{FEW_SHOT_EXAMPLE}\n\n"
        f"---\n\n"
        f"Now process the following window:\n\n"
        f"{state['evidence_text']}"
    )
    messages = [
        SystemMessage(content=SYSTEM_PROMPT),
        HumanMessage(content=user_content),
    ]
    response = model.invoke(messages)
    raw = response.content if isinstance(response.content, str) else str(response.content)
    logger.debug(
        f"[{state['species']}] LLM pairing response for {state['seqid']}:\n{raw[:400]}…"
    )
    return {"llm_raw_response": raw}


def node_parse_response(state: PairingState) -> dict:
    """
    Parse LLM JSON, validate indices, build locus candidates.
    Sets used_fallback=True and parse_error on failure.
    """
    raw = state.get("llm_raw_response", "")
    ntd_hits = state["ntd_hits"]
    ctd_hits = state["ctd_hits"]
    species = state["species"]
    seqid = state["seqid"]

    # Extract first JSON object from response
    match = re.search(r"\{[\s\S]*\}", raw)
    if not match:
        return {
            "used_fallback": True,
            "parse_error": "no_json_found",
            "final_loci": [],
        }

    try:
        result = json.loads(match.group())
    except json.JSONDecodeError as exc:
        return {
            "used_fallback": True,
            "parse_error": f"json_decode_error:{exc}",
            "final_loci": [],
        }

    loci: list[dict] = []
    seen_ntd: set[int] = set()
    seen_ctd: set[int] = set()

    for pair in result.get("pairs", []):
        ni = pair.get("ntd_index")
        ci = pair.get("ctd_index")

        # Validate indices
        if not isinstance(ni, int) or not isinstance(ci, int):
            continue
        if ni < 0 or ni >= len(ntd_hits) or ci < 0 or ci >= len(ctd_hits):
            logger.warning(
                f"[{species}] LLM returned out-of-range index ntd={ni} ctd={ci}, skipping"
            )
            continue
        if ni in seen_ntd or ci in seen_ctd:
            logger.warning(
                f"[{species}] LLM reused index ntd={ni} or ctd={ci}, skipping duplicate"
            )
            continue

        locus = _make_locus_shim(species, seqid, ntd=ntd_hits[ni], ctd=ctd_hits[ci])
        locus["llm_confidence"] = pair.get("confidence", "medium")
        locus["llm_needs_review"] = pair.get("needs_review", True)
        locus["llm_reasoning"] = pair.get("reasoning", "")
        locus["llm_evidence"] = pair.get("supporting_evidence", [])
        loci.append(locus)
        seen_ntd.add(ni)
        seen_ctd.add(ci)

    # Unpaired NTDs
    for ni in result.get("unpaired_ntds", []):
        if isinstance(ni, int) and 0 <= ni < len(ntd_hits) and ni not in seen_ntd:
            loci.append(_make_locus_shim(species, seqid, ntd=ntd_hits[ni], ctd=None))
            seen_ntd.add(ni)

    # Unpaired CTDs
    for ci in result.get("unpaired_ctds", []):
        if isinstance(ci, int) and 0 <= ci < len(ctd_hits) and ci not in seen_ctd:
            loci.append(_make_locus_shim(species, seqid, ntd=None, ctd=ctd_hits[ci]))
            seen_ctd.add(ci)

    # Safety net: any index not mentioned by LLM → single-terminus candidate
    for ni in range(len(ntd_hits)):
        if ni not in seen_ntd:
            logger.warning(f"[{species}] LLM omitted NTD[{ni}], emitting as N-terminal-only")
            loci.append(_make_locus_shim(species, seqid, ntd=ntd_hits[ni], ctd=None))
    for ci in range(len(ctd_hits)):
        if ci not in seen_ctd:
            logger.warning(f"[{species}] LLM omitted CTD[{ci}], emitting as C-terminal-only")
            loci.append(_make_locus_shim(species, seqid, ntd=None, ctd=ctd_hits[ci]))

    return {
        "final_loci": loci,
        "used_fallback": False,
        "parse_error": "",
    }


def node_fallback(state: PairingState) -> dict:
    """Greedy fallback — used when LLM call or parsing fails."""
    logger.warning(
        f"[{state['species']}] Falling back to greedy pairing for {state['seqid']} "
        f"(reason: {state.get('parse_error', 'unknown')})"
    )
    loci = _greedy_pair_local(
        species=state["species"],
        seqid=state["seqid"],
        ntd_hits=state["ntd_hits"],
        ctd_hits=state["ctd_hits"],
        min_span_bp=state["min_span_bp"],
        max_span_bp=state["max_span_bp"],
    )
    for locus in loci:
        locus["llm_confidence"] = None
        locus["llm_needs_review"] = None
        locus["llm_reasoning"] = "fallback_greedy"
        locus["llm_evidence"] = []
    return {"final_loci": loci, "used_fallback": True}


# ---------------------------------------------------------------------------
# Routing
# ---------------------------------------------------------------------------

def _route_after_parse(state: PairingState) -> str:
    if state.get("used_fallback"):
        return "fallback"
    return END


# ---------------------------------------------------------------------------
# Build graph (compiled once at module import)
# ---------------------------------------------------------------------------

def _build_graph() -> object:
    builder = StateGraph(PairingState)
    builder.add_node("build_evidence", node_build_evidence)
    builder.add_node("call_claude", node_call_claude)
    builder.add_node("parse_response", node_parse_response)
    builder.add_node("fallback", node_fallback)

    builder.add_edge(START, "build_evidence")
    builder.add_edge("build_evidence", "call_claude")
    builder.add_edge("call_claude", "parse_response")
    builder.add_conditional_edges(
        "parse_response",
        _route_after_parse,
        {"fallback": "fallback", END: END},
    )
    builder.add_edge("fallback", END)

    return builder.compile()


try:
    _GRAPH = _build_graph()
    _GRAPH_AVAILABLE = True
except Exception as _exc:
    logger.warning(f"LangGraph pairing graph could not be compiled: {_exc}")
    _GRAPH = None
    _GRAPH_AVAILABLE = False


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def run_pairing_agent(
    species: str,
    seqid: str,
    ntd_hits: list[dict],
    ctd_hits: list[dict],
    miniprot_hits: list[dict],
    min_span_bp: int,
    max_span_bp: int,
) -> list[dict]:
    """
    Run the LangGraph pairing agent for one ambiguous window.
    Falls back to greedy pairing if LLM is unavailable or fails.

    Returns a list of locus dicts (same schema as _make_locus output).
    """
    if not _GRAPH_AVAILABLE or not os.environ.get("DASHSCOPE_API_KEY"):
        logger.warning(
            f"[{species}] LLM pairing unavailable "
            "(DASHSCOPE_API_KEY not set or LangGraph not compiled). "
            "Using greedy fallback."
        )
        return _greedy_pair_local(
            species=species,
            seqid=seqid,
            ntd_hits=ntd_hits,
            ctd_hits=ctd_hits,
            min_span_bp=min_span_bp,
            max_span_bp=max_span_bp,
        )

    logger.info(
        f"[{species}] LLM pairing: {seqid} "
        f"({len(ntd_hits)} NTD, {len(ctd_hits)} CTD, {len(miniprot_hits)} miniprot hits)"
    )

    initial_state: PairingState = {
        "species": species,
        "seqid": seqid,
        "ntd_hits": ntd_hits,
        "ctd_hits": ctd_hits,
        "miniprot_hits": miniprot_hits,
        "min_span_bp": min_span_bp,
        "max_span_bp": max_span_bp,
        "evidence_text": "",
        "llm_raw_response": "",
        "final_loci": [],
        "used_fallback": False,
        "parse_error": "",
    }

    try:
        final_state = _GRAPH.invoke(initial_state)
    except Exception as exc:
        logger.error(f"[{species}] LangGraph pairing agent crashed: {exc}. Using greedy fallback.")
        loci = _greedy_pair_local(species, seqid, ntd_hits, ctd_hits, min_span_bp, max_span_bp)
        for locus in loci:
            locus["llm_confidence"] = None
            locus["llm_needs_review"] = None
            locus["llm_reasoning"] = f"fallback_crash:{exc}"
            locus["llm_evidence"] = []
        return loci

    return final_state.get("final_loci", [])
