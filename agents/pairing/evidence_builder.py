"""
Build a structured evidence packet for the LLM pairing node.

For each candidate NTD-CTD pair the packet includes:
  - pre-computed topological feasibility
  - spanning miniprot count (non-overlapping hits between NTD and CTD)
  - best miniprot Positive score in the interval
  - whether any same-terminus hit lies between the two endpoints
"""

from __future__ import annotations


# ---------------------------------------------------------------------------
# Helper: non-overlapping miniprot count in an interval
# ---------------------------------------------------------------------------

def _count_spanning_miniprot(
    interval_start: int,
    interval_end: int,
    strand: str,
    mp_hits: list[dict],
) -> tuple[int, float]:
    """
    Count non-overlapping miniprot hits wholly contained in [interval_start, interval_end]
    on the given strand. Returns (count, best_positive_score).

    Uses a greedy sweep (sort by start) to count the maximum non-overlapping set.
    """
    candidates = [
        h for h in mp_hits
        if h["strand"] == strand
        and h["start"] >= interval_start
        and h["end"] <= interval_end
    ]
    if not candidates:
        return 0, 0.0

    candidates.sort(key=lambda h: h["start"])
    count = 0
    last_end = -1
    best_positive = 0.0
    for h in candidates:
        if h["start"] > last_end:
            count += 1
            last_end = h["end"]
            best_positive = max(best_positive, h.get("positive", 0.0))
    return count, best_positive


# ---------------------------------------------------------------------------
# Helper: check for intervening same-terminus hits
# ---------------------------------------------------------------------------

def _count_intervening_hits(
    interval_start: int,
    interval_end: int,
    strand: str,
    all_hits: list[dict],
) -> int:
    """
    Count NTD or CTD hits (from any profile) that lie within
    (interval_start, interval_end) exclusive on the given strand.
    These intervening hits indicate that the interval is likely occupied by
    another gene, making the proposed pair suspicious.
    """
    return sum(
        1 for h in all_hits
        if h["strand"] == strand
        and h["start"] > interval_start
        and h["end"] < interval_end
    )


# ---------------------------------------------------------------------------
# Topology check
# ---------------------------------------------------------------------------

def _is_valid_topology(ntd: dict, ctd: dict, min_span_bp: int, max_span_bp: int) -> str:
    """
    Returns "OK" or a short reason string if the pair is invalid.
    """
    if ntd["strand"] != ctd["strand"]:
        return "strand_mismatch"
    strand = ntd["strand"]
    # On + strand: NTD start < CTD start (NTD is 5' of CTD)
    # On - strand: NTD start > CTD start (NTD is 5' in reverse)
    if strand == "+" and ntd["start"] >= ctd["start"]:
        return "wrong_order_plus"
    if strand == "-" and ntd["start"] <= ctd["start"]:
        return "wrong_order_minus"
    span = max(ntd["end"], ctd["end"]) - min(ntd["start"], ctd["start"]) + 1
    if span < min_span_bp:
        return f"span_too_small_{span}bp"
    if span > max_span_bp:
        return f"span_too_large_{span}bp"
    return "OK"


# ---------------------------------------------------------------------------
# E-value confidence label
# ---------------------------------------------------------------------------

_HIGH_EVALUE = 1e-10

def _evalue_label(e: float) -> str:
    return "HIGH" if e < _HIGH_EVALUE else "medium"


# ---------------------------------------------------------------------------
# Main: build_evidence_packet
# ---------------------------------------------------------------------------

def build_evidence_text(
    species: str,
    seqid: str,
    ntd_hits: list[dict],
    ctd_hits: list[dict],
    mp_hits: list[dict],
    min_span_bp: int,
    max_span_bp: int,
) -> str:
    """
    Format all evidence into a structured text block for the LLM.
    """
    all_hits = ntd_hits + ctd_hits
    window_start = min(h["start"] for h in all_hits)
    window_end = max(h["end"] for h in all_hits)

    lines: list[str] = []
    lines.append(f"=== WINDOW SUMMARY ===")
    lines.append(f"Species: {species}")
    lines.append(f"Scaffold: {seqid}  Window: {window_start:,} – {window_end:,} bp")
    lines.append("")

    # NTD hits
    lines.append("=== nHMMER NTD HITS ===")
    if ntd_hits:
        for i, h in enumerate(ntd_hits):
            lines.append(
                f"  [{i}] profile={h['profile']}  "
                f"pos={h['start']:,}-{h['end']:,}  strand={h['strand']}  "
                f"e_value={h['e_value']:.2e}  ({_evalue_label(h['e_value'])})"
            )
    else:
        lines.append("  (none)")
    lines.append("")

    # CTD hits
    lines.append("=== nHMMER CTD HITS ===")
    if ctd_hits:
        for i, h in enumerate(ctd_hits):
            lines.append(
                f"  [{i}] profile={h['profile']}  "
                f"pos={h['start']:,}-{h['end']:,}  strand={h['strand']}  "
                f"e_value={h['e_value']:.2e}  ({_evalue_label(h['e_value'])})"
            )
    else:
        lines.append("  (none)")
    lines.append("")

    # Miniprot hits in window
    mp_in_window = [
        h for h in mp_hits
        if h["start"] <= window_end and h["end"] >= window_start
    ]
    lines.append("=== MINIPROT ALIGNMENTS IN WINDOW ===")
    if mp_in_window:
        for h in sorted(mp_in_window, key=lambda x: x["start"]):
            lines.append(
                f"  ref={h['ref_protein']}  "
                f"pos={h['start']:,}-{h['end']:,}  strand={h['strand']}  "
                f"positive={h.get('positive', 0):.3f}  "
                f"identity={h.get('identity', 0):.3f}  rank={h.get('rank', '?')}"
            )
    else:
        lines.append("  (none)")
    lines.append("")

    # Pre-computed pair features
    lines.append("=== PRE-COMPUTED PAIR FEATURES ===")
    if ntd_hits and ctd_hits:
        for i, ntd in enumerate(ntd_hits):
            for j, ctd in enumerate(ctd_hits):
                topo = _is_valid_topology(ntd, ctd, min_span_bp, max_span_bp)
                if topo != "OK":
                    lines.append(f"  NTD[{i}]-CTD[{j}]: INVALID ({topo})")
                    continue

                # Interval for spanning analysis (inner boundaries)
                iv_start = min(ntd["start"], ctd["start"])
                iv_end = max(ntd["end"], ctd["end"])
                strand = ntd["strand"]

                span_count, best_pos = _count_spanning_miniprot(iv_start, iv_end, strand, mp_hits)
                intervening = _count_intervening_hits(iv_start, iv_end, strand, all_hits)
                dist = abs(
                    (ntd["start"] + ntd["end"]) // 2 - (ctd["start"] + ctd["end"]) // 2
                )

                # Type consistency
                ntd_type = ntd.get("type_prefix", "?")
                ctd_type = ctd.get("type_prefix", "?")
                if ntd_type == ctd_type:
                    type_note = f"type={ntd_type}/{ctd_type} (consistent)"
                else:
                    type_note = f"type={ntd_type}/{ctd_type} (CONFLICT)"

                evalue_note = (
                    "both_high_evalue"
                    if ntd["e_value"] < _HIGH_EVALUE and ctd["e_value"] < _HIGH_EVALUE
                    else "one_or_both_medium_evalue"
                )

                intervening_note = (
                    f"intervening_hits={intervening}"
                    if intervening > 0
                    else "no_intervening_hits"
                )

                lines.append(
                    f"  NTD[{i}]-CTD[{j}]: {type_note}, "
                    f"dist={dist:,}bp, span_bp={(iv_end - iv_start + 1):,}, "
                    f"spanning_miniprot={span_count} (best_positive={best_pos:.2f}), "
                    f"{intervening_note}, {evalue_note}"
                )
    else:
        lines.append(
            "  (no pairs possible — window has only "
            f"{'NTD' if ntd_hits else 'CTD'} hits)"
        )

    return "\n".join(lines)
