"""
typing_agent.py — automated spidroin identification agent.

Pipeline:
  1. Parse nHMMER GFF files (domain hits: NTD/CTD per species)
  2. Parse miniprot GFF files (reference protein alignments)
  3. Assemble spidroin locus candidates (pair NTD+CTD within max_span_bp)
  4. Classify each candidate (spidroin type + completeness)
  5. Join with miniprot evidence (best overlapping alignment)
  6. Add genome annotation overlap evidence
  7. Add RNA-seq transcription support evidence
  8. Compute confidence and needs_review flags
  9. Write results TSV to results/

Usage:
  pixi run python agents/typing_agent.py --help
  pixi run python agents/typing_agent.py --dry-run
"""

from pathlib import Path
from typing import Optional

from loguru import logger
import polars as pl
import typer
import yaml

from agents.evidence.annotation_overlap import annotate_loci as anno_annotate
from agents.evidence.annotation_overlap import upgrade_loci_from_annotation as anno_upgrade
from agents.evidence.codon_check import annotate_loci as codon_annotate
from agents.evidence.rna_support import annotate_loci as rna_annotate
from agents.evidence.rna_support import build_bw_index
from spider_silkome_module.config import PROCESSED_DATA_DIR, PROJ_ROOT, RAW_DATA_DIR

app = typer.Typer(add_completion=False)

# ---------------------------------------------------------------------------
# Default paths
# ---------------------------------------------------------------------------
_DEFAULT_NHMMER_DIR = PROCESSED_DATA_DIR / "nhmmer_search_parsed"
_DEFAULT_MINIPROT_DIR = PROCESSED_DATA_DIR / "miniprot_output"
_DEFAULT_ANNO_NEW = RAW_DATA_DIR / "spider_anno2"
_DEFAULT_ANNO_OLD = RAW_DATA_DIR / "01.ref_gff"
_DEFAULT_BW_CACHE = PROJ_ROOT / "cache" / "bw"
_DEFAULT_OUTPUT = PROCESSED_DATA_DIR / "typing_results"
_DEFAULT_RULES = PROJ_ROOT / "docs" / "typing_rules.yaml"


# ---------------------------------------------------------------------------
# GFF parsers
# ---------------------------------------------------------------------------

def _parse_attrs(attr_str: str) -> dict[str, str]:
    """Parse GFF3 attribute field into a dict."""
    attrs: dict[str, str] = {}
    for part in attr_str.strip().split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            attrs[k.strip()] = v.strip()
    return attrs


def parse_nhmmer_gff(gff_path: Path) -> list[dict]:
    """
    Parse a single species nHMMER GFF file.
    Only NTD and CTD feature lines are returned.
    E-value (lower = better) is read from the attribute field.
    """
    rows: list[dict] = []
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            seqid, _src, feature, start, end, score, strand, _phase, attr_str = parts[:9]
            if feature not in ("NTD", "CTD"):
                continue
            attrs = _parse_attrs(attr_str)
            profile = attrs.get("Name", "")
            species = attrs.get("Species", "")
            try:
                e_value = float(attrs["E-value"])
                hmmer_score = float(score)
            except (KeyError, ValueError):
                continue
            # Type prefix: "MaSp3_NTD" → "MaSp3"; "Spidroin_CTD" → "Spidroin"
            type_prefix = profile[: -(len(feature) + 1)]  # strip "_NTD" or "_CTD"
            rows.append(
                {
                    "seqid": seqid,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "terminus": feature,
                    "profile": profile,
                    "type_prefix": type_prefix,
                    "e_value": e_value,
                    "score": hmmer_score,
                    "species": species,
                }
            )
    return rows


def parse_miniprot_gff(gff_path: Path, species: str) -> list[dict]:
    """
    Parse a miniprot GFF file, extracting only mRNA-level records.
    Returns one row per alignment hit (Rank=1 preferred but all ranks returned).
    """
    rows: list[dict] = []
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            seqid, _src, feature, start, end, _score, strand, _phase, attr_str = parts[:9]
            if feature != "mRNA":
                continue
            attrs = _parse_attrs(attr_str)
            target = attrs.get("Target", "")
            ref_protein = target.split()[0] if target else ""
            try:
                rows.append(
                    {
                        "seqid": seqid,
                        "start": int(start),
                        "end": int(end),
                        "strand": strand,
                        "ref_protein": ref_protein,
                        "identity": float(attrs.get("Identity", 0)),
                        "positive": float(attrs.get("Positive", 0)),
                        "rank": int(attrs.get("Rank", 1)),
                        "species": species,
                    }
                )
            except ValueError:
                continue
    return rows


# ---------------------------------------------------------------------------
# Locus assembly
# ---------------------------------------------------------------------------

def _suppress_generic_hits(hits: list[dict]) -> list[dict]:
    """
    Within a single (species, seqid, terminus) pool, remove any hit with
    type_prefix='Spidroin' that overlaps a more specific (non-Spidroin) hit.
    Spidroin hits in regions with no specific competition are kept.
    """
    generic = [h for h in hits if h["type_prefix"] == "Spidroin"]
    specific = [h for h in hits if h["type_prefix"] != "Spidroin"]

    if not generic or not specific:
        return hits  # nothing to filter

    kept_generic = [
        g for g in generic
        if not any(g["start"] <= s["end"] and g["end"] >= s["start"] for s in specific)
    ]
    return specific + kept_generic


def _suppress_overlapping_hits(hits: list[dict]) -> list[dict]:
    """
    Within a single (species, seqid, terminus) pool, keep only the hit with
    the lowest e-value per overlapping group.

    Greedy algorithm (sorted by e-value ascending):
      - Accept a hit if it does not overlap any already-accepted hit.
      - Reject (suppress) any hit that overlaps an already-accepted hit.

    This prevents multiple profiles matching the same domain region from
    producing redundant single-terminus locus candidates.
    """
    sorted_hits = sorted(hits, key=lambda h: h["e_value"])
    kept: list[dict] = []
    for hit in sorted_hits:
        if not any(hit["start"] <= k["end"] and hit["end"] >= k["start"] for k in kept):
            kept.append(hit)
    return kept


def _cluster_hits_by_proximity(
    ntd_hits: list[dict],
    ctd_hits: list[dict],
    max_span_bp: int,
) -> list[tuple[list[dict], list[dict]]]:
    """
    Group NTD and CTD hits into local clusters by genomic proximity.

    Hits are sorted by start position. Two adjacent hits belong to the same
    cluster if the gap between them (start_next - end_current) is <= max_span_bp.
    Hits separated by more than max_span_bp can never be in the same valid pair,
    so they are placed in different clusters.

    Returns a list of (cluster_ntd_hits, cluster_ctd_hits) tuples.
    This ensures that hits from distinct genes spread across a chromosome are
    processed in independent local windows rather than as one giant ambiguous group.
    """
    if not ntd_hits and not ctd_hits:
        return []

    # Tag and sort all hits by genomic start
    tagged: list[tuple[int, str, dict]] = [
        (h["start"], "NTD", h) for h in ntd_hits
    ] + [
        (h["start"], "CTD", h) for h in ctd_hits
    ]
    tagged.sort(key=lambda x: x[0])

    clusters: list[list[tuple[str, dict]]] = []
    current: list[tuple[str, dict]] = []
    cluster_end = -1

    for start, terminus, hit in tagged:
        if not current:
            current.append((terminus, hit))
            cluster_end = hit["end"]
        elif start - cluster_end <= max_span_bp:
            current.append((terminus, hit))
            cluster_end = max(cluster_end, hit["end"])
        else:
            clusters.append(current)
            current = [(terminus, hit)]
            cluster_end = hit["end"]
    if current:
        clusters.append(current)

    return [
        (
            [h for t, h in c if t == "NTD"],
            [h for t, h in c if t == "CTD"],
        )
        for c in clusters
    ]


def _is_ambiguous(ntd_hits: list[dict], ctd_hits: list[dict]) -> bool:
    """
    Return True if this local cluster has multiple NTD or CTD hits after suppression,
    requiring LLM-based pairing instead of the simple greedy algorithm.
    """
    return len(ntd_hits) > 1 or len(ctd_hits) > 1


def _greedy_pair(
    species: str,
    seqid: str,
    ntd_hits: list[dict],
    ctd_hits: list[dict],
    min_span_bp: int,
    max_span_bp: int,
) -> list[dict]:
    """
    Simple greedy NTD-CTD pairing by minimum midpoint distance.
    Used for unambiguous windows (exactly 1 NTD or 1 CTD).
    """
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
        loci.append(_make_locus(species, seqid, ntd=ntd, ctd=ctd))
        ntd_used[i] = True
        ctd_used[j] = True

    for i, ntd in enumerate(ntd_hits):
        if not ntd_used[i]:
            loci.append(_make_locus(species, seqid, ntd=ntd, ctd=None))
    for j, ctd in enumerate(ctd_hits):
        if not ctd_used[j]:
            loci.append(_make_locus(species, seqid, ntd=None, ctd=ctd))

    return loci


def assemble_loci(
    hmmer_rows: list[dict],
    e_value_threshold: float,
    min_span_bp: int,
    max_span_bp: int,
    miniprot_rows: list[dict] | None = None,
    use_llm: bool = True,
) -> list[dict]:
    """
    Pair NTD and CTD hits from filtered nHMMER results into spidroin locus candidates.

    Algorithm per (species, seqid):
      1. Filter hits by E-value threshold (lower = better; keep e_value < threshold).
      2. Build all NTD–CTD pairs whose midpoints are within max_span_bp.
      3. Span filter: only pair if actual locus span >= min_span_bp (default 1 kb).
         Pairs below min_span_bp remain as single-terminal candidates.
      4. Greedy pairing: starting from the closest pair, assign each hit to at most one partner.
      5. Unpaired NTD → N-terminal-only candidate; unpaired CTD → C-terminal-only candidate.

    miniprot_rows: Optional list of miniprot alignment dicts. When provided and
        use_llm=True, ambiguous windows (multiple NTD or CTD hits after
        suppression) are resolved by the LangGraph pairing agent instead of
        the greedy algorithm.
    use_llm: Set to False to disable LLM pairing and always use greedy.

    Returns a list of locus dicts.
    """
    if not hmmer_rows:
        return []

    # Filter by E-value
    filtered = [r for r in hmmer_rows if r["e_value"] < e_value_threshold]
    if not filtered:
        return []

    # Group by (species, seqid)
    groups: dict[tuple[str, str], list[dict]] = {}
    for row in filtered:
        key = (row["species"], row["seqid"])
        groups.setdefault(key, []).append(row)

    loci: list[dict] = []

    for (species, seqid), hits in groups.items():
        ntd_hits = sorted(
            _suppress_overlapping_hits(
                _suppress_generic_hits([h for h in hits if h["terminus"] == "NTD"])
            ),
            key=lambda h: h["e_value"],
        )
        ctd_hits = sorted(
            _suppress_overlapping_hits(
                _suppress_generic_hits([h for h in hits if h["terminus"] == "CTD"])
            ),
            key=lambda h: h["e_value"],
        )

        # Split into per-gene local clusters before routing to LLM or greedy.
        # Hits separated by > max_span_bp can never be paired, so they go to
        # separate clusters and are processed independently.
        clusters = _cluster_hits_by_proximity(ntd_hits, ctd_hits, max_span_bp)

        # Pre-filter miniprot hits to this (species, seqid) once
        mp_seqid_hits = [
            r for r in (miniprot_rows or [])
            if r.get("species") == species and r.get("seqid") == seqid
        ]

        for cluster_ntd, cluster_ctd in clusters:
            if use_llm and _is_ambiguous(cluster_ntd, cluster_ctd):
                # Restrict miniprot to the cluster's genomic extent
                c_start = min(h["start"] for h in cluster_ntd + cluster_ctd)
                c_end = max(h["end"] for h in cluster_ntd + cluster_ctd)
                mp_cluster = [
                    r for r in mp_seqid_hits
                    if r["start"] <= c_end and r["end"] >= c_start
                ]
                from agents.pairing.pairing_agent import run_pairing_agent
                window_loci = run_pairing_agent(
                    species=species,
                    seqid=seqid,
                    ntd_hits=cluster_ntd,
                    ctd_hits=cluster_ctd,
                    miniprot_hits=mp_cluster,
                    min_span_bp=min_span_bp,
                    max_span_bp=max_span_bp,
                )
            else:
                window_loci = _greedy_pair(
                    species, seqid, cluster_ntd, cluster_ctd, min_span_bp, max_span_bp
                )
            # Mark loci from ambiguous clusters so that upgrade_loci_from_miniprot
            # Pass 2 does not second-guess the pairing decision already made here.
            if _is_ambiguous(cluster_ntd, cluster_ctd):
                for loc in window_loci:
                    loc["_from_ambiguous_cluster"] = True
            loci.extend(window_loci)

    return loci


def _make_locus(
    species: str,
    seqid: str,
    ntd: dict | None,
    ctd: dict | None,
) -> dict:
    """Construct a locus candidate dict from NTD and/or CTD hit records."""
    if ntd and ctd:
        strand = ntd["strand"]
    elif ntd:
        strand = ntd["strand"]
    else:
        assert ctd is not None
        strand = ctd["strand"]

    # Extend nHMMER CTD boundary by 3 bp outward to cover the stop codon.
    # nHMMER domain models end just before the stop codon, so the raw coordinates
    # do not include it.  +strand: push ctd_end forward; -strand: push ctd_start back.
    if ctd is not None:
        ctd_start_adj = ctd["start"] - 3 if strand == "-" else ctd["start"]
        ctd_end_adj = ctd["end"] + 3 if strand == "+" else ctd["end"]
    else:
        ctd_start_adj = None
        ctd_end_adj = None

    if ntd and ctd:
        locus_start = min(ntd["start"], ctd_start_adj)
        locus_end = max(ntd["end"], ctd_end_adj)
    elif ntd:
        locus_start, locus_end = ntd["start"], ntd["end"]
    else:
        locus_start, locus_end = ctd_start_adj, ctd_end_adj

    return {
        "species": species,
        "seqid": seqid,
        "start": locus_start,
        "end": locus_end,
        "strand": strand,
        "ntd_profile": ntd["profile"] if ntd else None,
        "ntd_type_prefix": ntd["type_prefix"] if ntd else None,
        "ntd_e_value": ntd["e_value"] if ntd else None,
        "ntd_start": ntd["start"] if ntd else None,
        "ntd_end": ntd["end"] if ntd else None,
        "ctd_profile": ctd["profile"] if ctd else None,
        "ctd_type_prefix": ctd["type_prefix"] if ctd else None,
        "ctd_e_value": ctd["e_value"] if ctd else None,
        "ctd_start": ctd_start_adj,
        "ctd_end": ctd_end_adj,
        "has_ntd": ntd is not None,
        "has_ctd": ctd is not None,
    }


# ---------------------------------------------------------------------------
# Type classification
# ---------------------------------------------------------------------------

def resolve_spidroin_type(
    ntd_prefix: str | None,
    ctd_prefix: str | None,
    type_mapping: dict[str, str],
    specificity_pairs: set[tuple[str, str]],
    type_to_family: dict[str, str],
) -> tuple[str, bool, str]:
    """
    Determine spidroin type from NTD and CTD type prefixes.
    Returns (spidroin_type, needs_review, reason).

    Resolution priority:
      1. Both absent → unknown, needs_review.
      2. One absent → use the present side.
      3. Both equal after mapping → use that type.
      4. One is a more specific sub-type (specificity_pairs) → use specific, no conflict.
      5. One side is "unknown" (generic Spidroin/hypo profile) → use the specific side,
         needs_review=True (generic NTD/CTD reduces confidence).
      6. Same type family (e.g. MaSp1 vs MaSp2) → resolve to family level, needs_review=True.
      7. Genuinely different families → record both, needs_review=True.
    """
    mapped_ntd = type_mapping.get(ntd_prefix, "unknown") if ntd_prefix else None
    mapped_ctd = type_mapping.get(ctd_prefix, "unknown") if ctd_prefix else None

    # Case 1: no hits at all
    if mapped_ntd is None and mapped_ctd is None:
        return "unknown", True, "no_terminus"

    # Case 2: only one side present
    if mapped_ntd is None:
        return mapped_ctd, mapped_ctd == "unknown", ""  # type: ignore[return-value]
    if mapped_ctd is None:
        return mapped_ntd, mapped_ntd == "unknown", ""

    # Case 3: both sides map to the same type
    if mapped_ntd == mapped_ctd:
        return mapped_ntd, mapped_ntd == "unknown", ""

    # Case 4: specificity hierarchy (e.g. MaSp3_NTD + MaSp_CTD → MaSp3)
    if (mapped_ntd, mapped_ctd) in specificity_pairs:
        return mapped_ntd, False, ""  # NTD is more specific
    if (mapped_ctd, mapped_ntd) in specificity_pairs:
        return mapped_ctd, False, ""  # CTD is more specific

    # Case 5: one side is "unknown" (generic profile), use the specific side
    if mapped_ntd == "unknown":
        return mapped_ctd, True, "NTD_generic_profile"
    if mapped_ctd == "unknown":
        return mapped_ntd, True, "CTD_generic_profile"

    # Case 6: same family but different subtypes (e.g. MaSp1 vs MaSp2)
    family_ntd = type_to_family.get(mapped_ntd)
    family_ctd = type_to_family.get(mapped_ctd)
    if family_ntd and family_ntd == family_ctd:
        return family_ntd, True, f"subtype_conflict:{mapped_ntd}_vs_{mapped_ctd}"

    # Case 7: genuinely different spidroin families → conflict
    return f"{mapped_ntd}/{mapped_ctd}", True, f"type_conflict:{mapped_ntd}_vs_{mapped_ctd}"


def resolve_type_by_context(
    loci: list[dict],
    all_hmmer_rows: list[dict],
    type_mapping: dict[str, str],
    window_bp: int = 50_000,
) -> list[dict]:
    """
    For loci with type conflicts (reason contains "conflict"), use neighboring
    nHMMER hits on the same scaffold as a voting signal.

    Looks at all nHMMER hits within ±window_bp of the locus boundaries, maps
    their type_prefix through type_mapping, and tallies votes.  If a clear
    majority type is found, it is recorded in the locus notes as
    "context_vote:<type>(N)" without overriding the existing classification —
    the annotation is advisory only.

    Modifies loci in-place and returns the list.
    """
    # Build per-(species, seqid) hit list once for all lookups
    hits_by_key: dict[tuple[str, str], list[dict]] = {}
    for row in all_hmmer_rows:
        key = (row["species"], row["seqid"])
        hits_by_key.setdefault(key, []).append(row)

    for locus in loci:
        # Only apply to conflict cases
        notes = locus.get("notes", "") or ""
        if "conflict" not in notes:
            continue

        key = (locus["species"], locus["seqid"])
        nearby = hits_by_key.get(key, [])
        locus_start = locus["start"]
        locus_end = locus["end"]

        # Collect types of neighboring hits (excluding the locus's own NTD/CTD positions)
        vote_counts: dict[str, int] = {}
        for hit in nearby:
            # Skip hits that overlap the locus itself (they're already used in classification)
            if hit["start"] <= locus_end and hit["end"] >= locus_start:
                continue
            # Only count hits within the window
            if (locus_start - window_bp) <= hit["start"] <= (locus_end + window_bp):
                hit_type = type_mapping.get(hit["type_prefix"], "unknown")
                if hit_type != "unknown":
                    vote_counts[hit_type] = vote_counts.get(hit_type, 0) + 1

        if not vote_counts:
            continue

        # Determine majority vote
        top_type = max(vote_counts, key=lambda t: vote_counts[t])
        top_count = vote_counts[top_type]
        total = sum(vote_counts.values())

        # Only record if majority is clear (>50% of votes)
        if top_count > total / 2:
            vote_note = f"context_vote:{top_type}({top_count}/{total})"
            sep = "; " if notes else ""
            locus["notes"] = notes + sep + vote_note

    return loci


def assess_completeness(
    has_ntd: bool,
    has_ctd: bool,
    ntd_e_value: float | None,
    ctd_e_value: float | None,
    high_e_threshold: float,
) -> str:
    """
    Classify spidroin completeness.
    Full_length requires both termini with high-confidence E-values.
    """
    ntd_confident = has_ntd and ntd_e_value is not None and ntd_e_value < high_e_threshold
    ctd_confident = has_ctd and ctd_e_value is not None and ctd_e_value < high_e_threshold

    if ntd_confident and ctd_confident:
        return "Full_length"
    if ntd_confident:
        return "N-terminal"
    if ctd_confident:
        return "C-terminal"
    # Both present but low confidence, or neither present
    # TODO: detect repeat-only once repeat-region evidence is integrated
    return "Repeat_only"


def compute_confidence(
    locus: dict,
    high_e_threshold: float,
    type_needs_review: bool,
    miniprot_positive_high: float,
) -> tuple[str, bool]:
    """
    Compute confidence level following evidence priority:
      nhmmer (primary) > miniprot > genome annotation > RNA-seq

    Returns (confidence, needs_review).
    RNA-seq and annotation evidence do not upgrade confidence but their
    absence is not penalized.
    """
    ntd_e = locus.get("ntd_e_value")
    ctd_e = locus.get("ctd_e_value")
    has_ntd = locus.get("has_ntd", False)
    has_ctd = locus.get("has_ctd", False)
    miniprot_positive = locus.get("miniprot_positive")

    ntd_high = has_ntd and ntd_e is not None and ntd_e < high_e_threshold
    ctd_high = has_ctd and ctd_e is not None and ctd_e < high_e_threshold

    # Start with type-resolution review flag; may only be set True, never cleared here
    needs_review = type_needs_review

    if ntd_high and ctd_high:
        # Both termini confirmed at high threshold: only high if no type issues
        confidence = "high" if not needs_review else "medium"
    elif ntd_high or ctd_high:
        # One terminus high-confidence
        confidence = "medium"
        needs_review = True
    elif has_ntd or has_ctd:
        # nhmmer hit present but below high threshold (between e_value and e_value_high)
        confidence = "medium"
        needs_review = True
    elif miniprot_positive is not None and miniprot_positive >= miniprot_positive_high:
        # No nhmmer; strong miniprot alignment is secondary evidence
        confidence = "medium"
        needs_review = True
    else:
        # Only annotation/RNA evidence or nothing
        confidence = "low"
        needs_review = True

    return confidence, needs_review


# ---------------------------------------------------------------------------
# Miniprot join
# ---------------------------------------------------------------------------

def join_miniprot(
    loci: list[dict],
    miniprot_rows: list[dict],
    positive_threshold: float,
) -> list[dict]:
    """
    For each locus, find the best overlapping miniprot mRNA alignment.
    Overlap criterion: miniprot seqid matches locus seqid AND regions overlap.
    Best = highest Positive score among rank-1 hits; falls back to any rank if no rank-1.

    Adds to each locus:
      miniprot_ref, miniprot_positive, miniprot_identity, miniprot_rank,
      miniprot_start, miniprot_end
    """
    # Group miniprot rows by species+seqid for fast lookup
    mp_by_key: dict[tuple[str, str], list[dict]] = {}
    for row in miniprot_rows:
        key = (row["species"], row["seqid"])
        mp_by_key.setdefault(key, []).append(row)

    for locus in loci:
        key = (locus["species"], locus["seqid"])
        candidates = mp_by_key.get(key, [])

        # Keep only hits that overlap the locus region
        overlapping = [
            r for r in candidates
            if r["start"] <= locus["end"] and r["end"] >= locus["start"]
        ]

        # Prefer rank-1 hits; fall back to all ranks if none found
        rank1 = [r for r in overlapping if r["rank"] == 1]
        pool = rank1 if rank1 else overlapping

        if pool:
            best = max(pool, key=lambda r: r["positive"])
            locus["miniprot_ref"] = best["ref_protein"]
            locus["miniprot_positive"] = best["positive"]
            locus["miniprot_identity"] = best["identity"]
            locus["miniprot_rank"] = best["rank"]
            locus["miniprot_start"] = best["start"]
            locus["miniprot_end"] = best["end"]
        else:
            locus["miniprot_ref"] = None
            locus["miniprot_positive"] = None
            locus["miniprot_identity"] = None
            locus["miniprot_rank"] = None
            locus["miniprot_start"] = None
            locus["miniprot_end"] = None

        # If no nHMMER-based type but miniprot has a strong hit, infer type from ref_protein name
        # Example: "21_MaSp1_CTD_873aa" → contains "MaSp1"
        if locus.get("spidroin_type") == "unknown" and locus.get("miniprot_ref"):
            # TODO: improve type inference from miniprot ref protein name
            #       using a regex or the type_mapping lookup on the ref protein tokens
            pass

    return loci


# ---------------------------------------------------------------------------
# Species abbreviation map (for Spidroin_ID generation)
# ---------------------------------------------------------------------------

def build_abbrev_map(species_ids: list[str]) -> dict[str, str]:
    """
    Build a collision-free abbreviation for each species_id.

    species_id format: "{NNN}.{Genus_species}" (e.g. "064.Araneus_ventricosus")
    Base abbreviation: Genus[:2] + species_epithet[:2] (e.g. "Arve")

    Collision resolution: progressively extend genus prefix, then species prefix,
    until all abbreviations are unique.  If still not unique after exhausting all
    extensions, append a numeric suffix (e.g. "Arve2").
    """
    # Extract (genus, species_epithet) from each species_id
    parsed: dict[str, tuple[str, str]] = {}
    for sid in species_ids:
        name_part = sid.split(".", 1)[-1]  # strip NNN prefix
        parts = name_part.split("_", 1)
        genus = parts[0] if parts else name_part
        epithet = parts[1] if len(parts) > 1 else ""
        parsed[sid] = (genus, epithet)

    abbrev_map: dict[str, str] = {}

    # Try progressively longer combinations until no collisions remain
    # combination format: genus[:g_len] + epithet[:e_len]
    remaining = list(species_ids)
    g_len, e_len = 2, 2

    while remaining:
        candidate: dict[str, str] = {}
        for sid in remaining:
            genus, epithet = parsed[sid]
            abbrev = (genus[:g_len] + epithet[:e_len]).replace("_", "")
            candidate[sid] = abbrev

        # Detect collisions
        used: dict[str, list[str]] = {}
        for sid, abbrev in candidate.items():
            used.setdefault(abbrev, []).append(sid)

        resolved: list[str] = []
        for sid in remaining:
            if len(used[candidate[sid]]) == 1:
                abbrev_map[sid] = candidate[sid]
                resolved.append(sid)

        remaining = [sid for sid in remaining if sid not in resolved]

        # Extend prefix lengths for the next round
        if remaining:
            genus_max = max(len(parsed[sid][0]) for sid in remaining)
            epithet_max = max(len(parsed[sid][1]) for sid in remaining)
            if g_len < genus_max:
                g_len += 1
            elif e_len < epithet_max:
                e_len += 1
            else:
                # Cannot disambiguate further; add numeric suffix
                suffix_counter: dict[str, int] = {}
                for sid in remaining:
                    genus, epithet = parsed[sid]
                    base = (genus[:g_len] + epithet[:e_len]).replace("_", "")
                    suffix_counter[base] = suffix_counter.get(base, 1) + 1
                    abbrev_map[sid] = f"{base}{suffix_counter[base]}"
                remaining = []

    return abbrev_map


# ---------------------------------------------------------------------------
# Per-species worker
# ---------------------------------------------------------------------------

def _process_species(
    species_id: str,
    sp_dir: Path,
    miniprot_dir: Path,
    e_value_threshold: float,
    e_value_high: float,
    min_span_bp: int,
    max_span_bp: int,
    positive_high: float,
    type_mapping: dict[str, str],
    specificity_pairs: set[tuple[str, str]],
    type_to_family: dict[str, str],
    anno_new: Optional[Path],
    anno_old: Optional[Path],
    bw_cache: Path,
    bw_index: dict | None,
    skip_anno: bool,
    skip_codon: bool,
    skip_rna: bool,
    use_llm: bool = True,
) -> list[dict]:
    """Run all pipeline steps for one species. Returns annotated locus dicts."""
    logger.info(f"[{species_id}] Starting")

    # Step 1: Parse nHMMER GFF
    gff = sp_dir / f"{sp_dir.name}.gff"
    if not gff.exists():
        logger.warning(f"[{species_id}] nHMMER GFF not found: {gff}")
        return []
    hmmer_rows = parse_nhmmer_gff(gff)

    # Step 2: Parse miniprot GFF
    miniprot_rows: list[dict] = []
    mp_file = miniprot_dir / f"{species_id}.gff"
    if mp_file.exists():
        miniprot_rows = parse_miniprot_gff(mp_file, species_id)

    # Step 3: Assemble loci
    loci = assemble_loci(
        hmmer_rows,
        e_value_threshold,
        min_span_bp,
        max_span_bp,
        miniprot_rows=miniprot_rows,
        use_llm=use_llm,
    )
    if not loci:
        logger.info(f"[{species_id}] No loci assembled")
        return []

    # Step 4: Classify (type + completeness)
    for locus in loci:
        spidroin_type, type_needs_review, type_reason = resolve_spidroin_type(
            locus.get("ntd_type_prefix"),
            locus.get("ctd_type_prefix"),
            type_mapping,
            specificity_pairs,
            type_to_family,
        )
        locus["spidroin_type"] = spidroin_type
        locus["_type_needs_review"] = type_needs_review
        locus["_type_reason"] = type_reason
        locus["completeness"] = assess_completeness(
            locus["has_ntd"],
            locus["has_ctd"],
            locus.get("ntd_e_value"),
            locus.get("ctd_e_value"),
            e_value_high,
        )
        locus["notes"] = type_reason if type_reason else ""
        locus["confidence"] = None
        locus["needs_review"] = type_needs_review

    # Step 5: Context vote (neighbor hits, same species only)
    loci = resolve_type_by_context(loci, hmmer_rows, type_mapping)

    # Step 5b: LLM-based type conflict resolution
    # Only runs when use_llm=True and a cross-family conflict is detected.
    if use_llm:
        from agents.pairing.type_resolver import resolve_type_conflict_llm
        for locus in loci:
            reason = locus.get("_type_reason", "")
            if "type_conflict" not in reason:
                continue
            conflict_type = locus["spidroin_type"]
            resolved, nr, llm_reason = resolve_type_conflict_llm(
                locus, conflict_type, hmmer_rows, miniprot_rows
            )
            if resolved != conflict_type:
                logger.info(
                    f"[{species_id}] Type resolved {conflict_type} → {resolved} "
                    f"at {locus.get('seqid')}:{locus.get('start')}"
                )
            locus["spidroin_type"] = resolved
            locus["_type_needs_review"] = nr
            existing = locus.get("notes", "") or ""
            sep = "; " if existing else ""
            locus["notes"] = existing + sep + f"llm_type_resolve:{llm_reason[:120]}"

    # Step 6: Join miniprot
    loci = join_miniprot(loci, miniprot_rows, positive_high)

    # Step 6b: Type rescue + completeness upgrade from Miniprot evidence
    loci = upgrade_loci_from_miniprot(loci, miniprot_rows, max_span_bp, type_mapping)

    # Step 7: Finalize confidence
    for locus in loci:
        confidence, needs_review = compute_confidence(
            locus,
            e_value_high,
            locus["_type_needs_review"],
            positive_high,
        )
        locus["confidence"] = confidence
        locus["needs_review"] = needs_review

    # Step 8: Annotation overlap
    if not skip_anno:
        effective_anno_new = anno_new if (anno_new and anno_new.exists()) else None
        effective_anno_old = anno_old if (anno_old and anno_old.exists()) else None
        if effective_anno_new or effective_anno_old:
            loci = anno_annotate(loci, effective_anno_new, effective_anno_old)
        else:
            for locus in loci:
                locus.setdefault("anno_gene_id", None)
                locus.setdefault("anno_overlap_frac", None)
                locus.setdefault("anno_start", None)
                locus.setdefault("anno_end", None)
    else:
        for locus in loci:
            locus["anno_gene_id"] = None
            locus["anno_overlap_frac"] = None
            locus["anno_start"] = None
            locus["anno_end"] = None

    # Step 8b: Annotation-guided completeness upgrade for N-terminal-only loci
    # Runs after anno_annotate (which populates anno_start/anno_end) so that
    # the proxy ctd_start/ctd_end set here are available for Step 9 codon check.
    if not skip_anno:
        loci = anno_upgrade(loci)

    # Step 9: Codon check (genome FASTA loaded and freed inside codon_annotate)
    if not skip_codon:
        fasta_root_old = anno_old if (anno_old and anno_old.exists()) else None
        fasta_root_new = anno_new if (anno_new and anno_new.exists()) else None
        if fasta_root_old or fasta_root_new:
            loci = codon_annotate(loci, fasta_root_old, fasta_root_new=fasta_root_new)
        else:
            for locus in loci:
                locus["ntd_start_codon"] = None
                locus["ctd_stop_codon"] = None
    else:
        for locus in loci:
            locus["ntd_start_codon"] = None
            locus["ctd_stop_codon"] = None

    # Step 10: RNA support
    if not skip_rna and bw_index is not None:
        loci = rna_annotate(loci, bw_cache, bw_index)
    else:
        for locus in loci:
            locus["bgi_rna_support"] = None
            locus["ont_rna_support"] = None
            locus["rna_mean_signal"] = None

    logger.info(f"[{species_id}] Done: {len(loci)} loci")
    return loci


# ---------------------------------------------------------------------------
# Feishu integration placeholder
# ---------------------------------------------------------------------------

def export_to_feishu(df: pl.DataFrame, table_id: str) -> None:
    # TODO: implement Feishu multi-dimensional table API write
    raise NotImplementedError("Feishu export not yet implemented")


# ---------------------------------------------------------------------------
# Miniprot-guided locus upgrade
# ---------------------------------------------------------------------------

def _parse_miniprot_ref_type(ref: str) -> tuple[str, str] | None:
    """
    Parse (type_prefix, terminus) from a reference protein name.
    Naming convention: '{id}_{Type}_{CTD|NTD}_{length}aa[.rank]'
    Example: '2666_MiSp_CTD_426aa.1' → ('MiSp', 'CTD')
    Returns None if the name does not match the expected pattern.
    """
    ref_clean = ref.split(".")[0]  # strip rank suffix
    parts = ref_clean.split("_")
    for i, part in enumerate(parts):
        if part in ("CTD", "NTD") and i >= 2:
            return "_".join(parts[1:i]), part
    return None


def upgrade_loci_from_miniprot(
    loci: list[dict],
    miniprot_rows: list[dict],
    max_span_bp: int,
    type_mapping: dict[str, str],
) -> list[dict]:
    """
    Two upgrade passes using Miniprot evidence:

    Pass 1 — Type rescue for loci where nHMMER type resolved to 'unknown':
      Parse the type from miniprot_ref (best Miniprot hit already attached by
      join_miniprot). If a valid specific type is found, override spidroin_type.

    Pass 2 — Completeness upgrade for single-terminus loci:
      Type is NOT modified here; nHMMER already determined it.
      Searches for Miniprot hits of the opposing terminus type within max_span_bp
      downstream (N-terminal loci) or upstream (C-terminal loci).
      If found: sets ctd_start/ctd_end (or ntd_start/ntd_end), extends locus
      coordinates, upgrades completeness to Full_length, adds a note.
    """
    mp_index: dict[tuple[str, str], list[dict]] = {}
    for row in miniprot_rows:
        key = (row["species"], row["seqid"])
        mp_index.setdefault(key, []).append(row)

    # Build interval index of all assembled loci (pre-upgrade) to avoid stealing
    # miniprot hits that already belong to a neighbouring locus on the same scaffold.
    locus_intervals: dict[tuple[str, str], list[tuple[int, int]]] = {}
    for loc in loci:
        ikey = (loc["species"], loc["seqid"])
        locus_intervals.setdefault(ikey, []).append((loc["start"], loc["end"]))

    for locus in loci:
        # ── Pass 1: type rescue for 'unknown' loci ──────────────────────
        if locus.get("spidroin_type") == "unknown":
            ref = locus.get("miniprot_ref")
            if ref:
                parsed = _parse_miniprot_ref_type(ref)
                if parsed:
                    inferred_prefix, _ = parsed
                    mapped = type_mapping.get(inferred_prefix, inferred_prefix)
                    if mapped and mapped != "unknown":
                        locus["spidroin_type"] = mapped
                        locus["needs_review"] = True
                        existing = locus.get("notes") or ""
                        sep = "; " if existing else ""
                        locus["notes"] = existing + sep + "type_from_miniprot"

        # ── Pass 2: completeness upgrade for single-terminus loci ────────
        has_ntd = locus.get("has_ntd", False)
        has_ctd = locus.get("has_ctd", False)
        if has_ntd == has_ctd:
            continue  # Full_length or Repeat_only: nothing to upgrade

        # Skip loci from ambiguous clusters — their pairing was already decided
        # by LLM or greedy in assemble_loci; do not second-guess that decision.
        if locus.get("_from_ambiguous_cluster"):
            continue

        key = (locus["species"], locus["seqid"])
        candidates = mp_index.get(key, [])
        strand = locus.get("strand", "+")

        # Intervals of OTHER loci on this scaffold (excluding self) used to
        # prevent a miniprot hit that belongs to a neighbouring locus from
        # being used to upgrade completeness of this locus.
        locus_own = (locus["start"], locus["end"])
        other_intervals = [
            iv for iv in locus_intervals.get(key, [])
            if iv != locus_own
        ]

        def _in_other_locus(r: dict) -> bool:
            for ivs, ive in other_intervals:
                if r["start"] <= ive and r["end"] >= ivs:
                    return True
            return False

        if has_ntd and not has_ctd:
            # Search downstream for CTD-type Miniprot hits
            if strand == "+":
                window = (locus["end"], locus["end"] + max_span_bp)
            else:
                window = (locus["start"] - max_span_bp, locus["start"])
            hits = [
                r for r in candidates
                if r.get("strand") == strand
                and "_CTD_" in (r.get("ref_protein") or "")
                and r["start"] >= window[0]
                and r["end"] <= window[1]
                and not _in_other_locus(r)
            ]
            if not hits:
                continue
            inferred_start = min(r["start"] for r in hits)
            inferred_end = max(r["end"] for r in hits)
            locus["ctd_start"] = inferred_start
            locus["ctd_end"] = inferred_end
            locus["end"] = max(locus["end"], inferred_end)
            locus["completeness"] = "Full_length"
            locus["needs_review"] = True
            existing = locus.get("notes") or ""
            sep = "; " if existing else ""
            locus["notes"] = existing + sep + f"completeness_from_miniprot_ctd({len(hits)}hits)"

        else:  # has_ctd and not has_ntd
            # Search upstream for NTD-type Miniprot hits
            if strand == "+":
                window = (locus["start"] - max_span_bp, locus["start"])
            else:
                window = (locus["end"], locus["end"] + max_span_bp)
            hits = [
                r for r in candidates
                if r.get("strand") == strand
                and "_NTD_" in (r.get("ref_protein") or "")
                and r["start"] >= window[0]
                and r["end"] <= window[1]
                and not _in_other_locus(r)
            ]
            if not hits:
                continue
            inferred_start = min(r["start"] for r in hits)
            inferred_end = max(r["end"] for r in hits)
            locus["ntd_start"] = inferred_start
            locus["ntd_end"] = inferred_end
            locus["start"] = min(locus["start"], inferred_start)
            locus["completeness"] = "Full_length"
            locus["needs_review"] = True
            existing = locus.get("notes") or ""
            sep = "; " if existing else ""
            locus["notes"] = existing + sep + f"completeness_from_miniprot_ntd({len(hits)}hits)"

    return loci


# ---------------------------------------------------------------------------
# GFF output
# ---------------------------------------------------------------------------

_GFF_ESCAPE = str.maketrans({";": "%3B", "=": "%3D", "&": "%26", ",": "%2C"})


def _gff_escape(value: str) -> str:
    return str(value).translate(_GFF_ESCAPE)


def _gff_lines_from_rows(rows: list[dict], species_id: str) -> list[str]:
    """
    Build GFF3 content lines from a list of row dicts (locus or TSV-derived).
    Expected keys: seqid, start, end, strand, Spidroin_ID, spidroin_type,
    Hint_type, confidence, ntd_e_value, ctd_e_value, Note,
    ntd_start, ntd_end, ntd_profile, ctd_start, ctd_end, ctd_profile.
    NTD/CTD sub-feature lines are emitted when ntd_start / ctd_start are non-null.
    """
    lines: list[str] = [
        "##gff-version 3",
        f"# Generated by spidroin_agent | species: {species_id}",
    ]
    sorted_rows = sorted(
        rows,
        key=lambda r: (str(r.get("seqid", "")), int(r.get("start") or 0)),
    )
    for row in sorted_rows:
        seqid = row.get("seqid", ".")
        start = row.get("start", ".")
        end = row.get("end", ".")
        strand = row.get("strand") or "."
        spidroin_id = row.get("Spidroin_ID", ".")
        spidroin_type = row.get("spidroin_type") or "Unknown"
        hint_type = row.get("Hint_type") or ""
        confidence = row.get("confidence") or ""
        ntd_ev = row.get("ntd_e_value") or ""
        ctd_ev = row.get("ctd_e_value") or ""
        note = _gff_escape(row.get("Note") or "")

        attrs = (
            f"ID={spidroin_id}"
            f";Name={_gff_escape(spidroin_type)}"
            f";hint_type={hint_type}"
            f";confidence={confidence}"
            f";ntd_evalue={ntd_ev}"
            f";ctd_evalue={ctd_ev}"
        )
        if note:
            attrs += f";Note={note}"

        lines.append(f"{seqid}\tspidroin_agent\tspidroin_gene\t{start}\t{end}\t.\t{strand}\t.\t{attrs}")

        ntd_start = row.get("ntd_start")
        if ntd_start is not None and ntd_start != "":
            ntd_profile = row.get("ntd_profile") or ""
            lines.append(
                f"{seqid}\tspidroin_agent\tNTD\t{ntd_start}\t{row.get('ntd_end', '.')}\t.\t{strand}\t."
                f"\tParent={spidroin_id};Name={_gff_escape(ntd_profile)}"
            )

        ctd_start = row.get("ctd_start")
        if ctd_start is not None and ctd_start != "":
            ctd_profile = row.get("ctd_profile") or ""
            lines.append(
                f"{seqid}\tspidroin_agent\tCTD\t{ctd_start}\t{row.get('ctd_end', '.')}\t.\t{strand}\t."
                f"\tParent={spidroin_id};Name={_gff_escape(ctd_profile)}"
            )

    return lines


def write_loci_gff(sp_loci: list[dict], species_id: str, out_path: Path) -> None:
    """Write a GFF3 file for one species from enriched locus dicts."""
    out_path.write_text("\n".join(_gff_lines_from_rows(sp_loci, species_id)) + "\n")


def write_fasta_and_hints(
    sp_loci: list[dict],
    species_id: str,
    out_dir: Path,
    fasta_root: Path | None,
    fasta_root_new: Path | None,
) -> None:
    """
    Write spidroin_full_length.fasta and hints.gff for full-length loci with
    verified start and stop codons.

    Only loci satisfying all three conditions are included:
      • completeness == "Full_length"
      • ntd_start_codon is True  (ATG confirmed by codon_check)
      • ctd_stop_codon  is True  (stop codon confirmed by codon_check)

    FASTA format
    ============
    Sequences are extracted from the genome and oriented 5'→3' (reverse-complement
    for '−' strand loci), so every entry starts with ATG and ends with a stop codon.
    Header: >{Spidroin_ID} {seqid}:{start}-{end}({strand}) type={spidroin_type}
    Lines wrapped at 60 characters.

    hints.gff format (Augustus)
    ===========================
    Coordinates are relative to the extracted FASTA (1-based, always '+' strand).
    Start codon: positions 1–3.
    Stop  codon: positions (seq_len − 2)–seq_len.
    Attributes:  grp=1;pri=4;src=M  (one gene per sequence → grp always 1)

    Files written only when at least one qualifying locus is found.
    """
    import pysam

    from agents.evidence.codon_check import _find_genome_fasta, _revcomp

    candidates = [
        loc for loc in sp_loci
        if loc.get("completeness") == "Full_length"
        and loc.get("ntd_start_codon") is True
        and loc.get("ctd_stop_codon") is True
    ]
    if not candidates:
        return

    fasta_path = _find_genome_fasta(fasta_root, species_id, fasta_root_new)
    if fasta_path is None:
        logger.warning(f"[{species_id}] No genome FASTA found; skipping FASTA/hints output")
        return

    try:
        fasta = pysam.FastaFile(str(fasta_path))
    except Exception as exc:
        logger.warning(f"[{species_id}] Cannot open genome FASTA for sequence extraction: {exc}")
        return

    fasta_records: list[str] = []
    hints_lines:   list[str] = []

    with fasta:
        for locus in candidates:
            seqid        = locus["seqid"]
            strand       = locus.get("strand", "+")
            spidroin_id  = locus.get("Spidroin_ID", ".")
            spidroin_type = locus.get("spidroin_type", "Unknown")

            # Use confirmed codon positions as extraction boundaries so the
            # sequence always starts with ATG and ends with the stop codon.
            ntd_start = locus.get("ntd_start")
            ntd_end   = locus.get("ntd_end")
            ctd_start = locus.get("ctd_start")
            ctd_end   = locus.get("ctd_end")

            if strand == "+":
                if ntd_start is None or ctd_end is None:
                    continue
                if ntd_start > ctd_end:
                    # CTD is upstream of NTD — discordant domain order; skip
                    logger.debug(
                        f"[{species_id}] {spidroin_id}: discordant domain order on + strand; skipping"
                    )
                    continue
                seq_gstart = ntd_start   # 1-based, first base of ATG
                seq_gend   = ctd_end     # 1-based, last base of stop codon
            else:
                if ctd_start is None or ntd_end is None:
                    continue
                if ntd_end < ctd_start:
                    # NTD is downstream of CTD — discordant domain order; skip
                    logger.debug(
                        f"[{species_id}] {spidroin_id}: discordant domain order on - strand; skipping"
                    )
                    continue
                seq_gstart = ctd_start   # 1-based, lowest coord (3' of mRNA)
                seq_gend   = ntd_end     # 1-based, highest coord (5' of mRNA = ATG)

            if seqid not in fasta.references:
                logger.debug(f"[{species_id}] {seqid} not in FASTA; skipping {spidroin_id}")
                continue

            chrom_len   = fasta.get_reference_length(seqid)
            fetch_start = max(0, seq_gstart - 1)    # 0-based
            fetch_end   = min(chrom_len, seq_gend)  # 0-based exclusive

            try:
                seq = fasta.fetch(seqid, fetch_start, fetch_end).upper()
            except Exception as exc:
                logger.debug(f"[{species_id}] Fetch error {seqid}:{seq_gstart}-{seq_gend}: {exc}")
                continue

            if strand == "-":
                seq = _revcomp(seq)

            seq_len = len(seq)
            if seq_len < 6:
                # Too short to contain both start and stop codon; skip silently
                continue

            # FASTA record
            header = (
                f">{spidroin_id} {seqid}:{seq_gstart}-{seq_gend}({strand})"
                f" type={spidroin_type}"
            )
            wrapped = "\n".join(seq[i : i + 60] for i in range(0, seq_len, 60))
            fasta_records.append(f"{header}\n{wrapped}")

            # Augustus hints (coordinates relative to extracted FASTA, always '+')
            stop_start = seq_len - 2   # 1-based
            stop_end   = seq_len       # 1-based
            attrs = "grp=1;pri=4;src=M"
            hints_lines.append(
                f"{spidroin_id}\tspidroin_agent\tstart\t1\t3\t.\t+\t0\t{attrs}"
            )
            hints_lines.append(
                f"{spidroin_id}\tspidroin_agent\tstop\t{stop_start}\t{stop_end}\t.\t+\t0\t{attrs}"
            )

    if fasta_records:
        (out_dir / "spidroin_full_length.fasta").write_text(
            "\n".join(fasta_records) + "\n"
        )
    if hints_lines:
        (out_dir / "hints.gff").write_text(
            "\n".join(hints_lines) + "\n"
        )

    logger.info(
        f"[{species_id}] {len(fasta_records)} full-length sequences → "
        f"spidroin_full_length.fasta + hints.gff"
    )


def write_gff_from_tsv(tsv_path: Path, out_path: Path) -> None:
    """
    Regenerate a GFF3 file from an existing species TSV without re-running the pipeline.
    Requires the TSV to include ntd_start/ntd_end/ctd_start/ctd_end columns
    (written by this agent since those columns were added to diag_cols).
    """
    df = pl.read_csv(tsv_path, separator="\t", missing_utf8_is_empty_string=True, infer_schema_length=0)
    species_id = tsv_path.stem
    rename = {"Chr": "seqid", "Start": "start", "End": "end", "Stran": "strand"}
    df = df.rename({k: v for k, v in rename.items() if k in df.columns})
    lines = _gff_lines_from_rows(df.to_dicts(), species_id)
    out_path.write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# Main CLI
# ---------------------------------------------------------------------------

@app.command()
def main(
    nhmmer_dir: Path = typer.Option(
        _DEFAULT_NHMMER_DIR,
        help="Directory of nHMMER parsed GFF files (<species>/<species>.gff)",
    ),
    miniprot_dir: Path = typer.Option(
        _DEFAULT_MINIPROT_DIR,
        help="Directory of miniprot GFF files (<species>.gff)",
    ),
    anno_new: Optional[Path] = typer.Option(
        _DEFAULT_ANNO_NEW,
        help="spider_anno2 root directory (<NNNspd>/gene.gff). Pass empty string to disable.",
    ),
    anno_old: Optional[Path] = typer.Option(
        _DEFAULT_ANNO_OLD,
        help="01.ref_gff root directory (fallback annotation + genome FASTA).",
    ),
    bw_cache: Path = typer.Option(
        _DEFAULT_BW_CACHE,
        help="Local directory for caching downloaded BigWig files.",
    ),
    output: Path = typer.Option(
        _DEFAULT_OUTPUT,
        help="Output directory. One TSV file per species is written here.",
    ),
    rules_file: Path = typer.Option(
        _DEFAULT_RULES,
        help="YAML file with typing rules and thresholds.",
    ),
    min_span: Optional[int] = typer.Option(
        None,
        help="Min locus span (bp) to form an NTD+CTD pair. Overrides rules file.",
    ),
    max_span: Optional[int] = typer.Option(
        None,
        help="Max NTD–CTD midpoint distance (bp) to form a locus. Overrides rules file.",
    ),
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        help="Parse and classify without writing output files.",
    ),
    resume: bool = typer.Option(
        False,
        "--resume",
        help="Skip species whose output TSV already exists in the output directory.",
    ),
    skip_rna: bool = typer.Option(
        False,
        "--skip-rna",
        help="Skip RNA-seq BigWig evidence (faster, no MinIO access needed).",
    ),
    skip_anno: bool = typer.Option(
        False,
        "--skip-anno",
        help="Skip genome annotation overlap evidence.",
    ),
    skip_codon: bool = typer.Option(
        False,
        "--skip-codon",
        help="Skip start/stop codon boundary check (no genome FASTA access).",
    ),
    workers: int = typer.Option(
        1,
        "--workers",
        "-j",
        help="Number of parallel worker threads. 1 = sequential (default).",
        min=1,
    ),
    use_llm: bool = typer.Option(True, "--use-llm/--no-llm", help="Use LangGraph LLM pairing for ambiguous windows (requires ANTHROPIC_API_KEY)."),
) -> None:
    """
    Automated spidroin identification agent.
    Parses nHMMER + miniprot results, assembles locus candidates, and writes a TSV summary
    aligned with the Feishu multi-dimensional table schema (docs/feishu_column.md).
    """
    # ------------------------------------------------------------------
    # Load rules
    # ------------------------------------------------------------------
    if not rules_file.exists():
        logger.error(f"Rules file not found: {rules_file}")
        raise typer.Exit(1)

    with open(rules_file) as fh:
        rules = yaml.safe_load(fh)

    thresholds = rules.get("thresholds", {})
    nhmmer_cfg = thresholds.get("nhmmer", {})
    miniprot_cfg = thresholds.get("miniprot", {})

    e_value_threshold: float = nhmmer_cfg.get("e_value", 1e-4)
    e_value_high: float = nhmmer_cfg.get("e_value_high", 1e-10)
    min_span_bp: int = min_span if min_span is not None else int(thresholds.get("min_span_bp", 1_000))
    max_span_bp: int = max_span if max_span is not None else int(thresholds.get("max_span_bp", 150_000))
    positive_high: float = miniprot_cfg.get("positive_high", 0.6)

    type_mapping: dict[str, str] = rules.get("type_mapping", {})
    specificity_list: list[list[str]] = rules.get("type_specificity", [])
    specificity_pairs: set[tuple[str, str]] = {(a, b) for a, b in specificity_list}

    # Build reverse lookup: member type → family name (e.g. "MaSp1" → "MaSp")
    type_to_family: dict[str, str] = {
        member: family
        for family, members in rules.get("type_families", {}).items()
        for member in members
        if member != family  # exclude the family name itself (no self-mapping needed)
    }

    logger.info(
        f"Thresholds: nhmmer e_value<{e_value_threshold}, high<{e_value_high}, "
        f"span=[{min_span_bp}, {max_span_bp}] bp, miniprot positive_high>={positive_high}"
    )

    # ------------------------------------------------------------------
    # Build bw_index once (RNA enabled path) — must not be called per-thread
    # ------------------------------------------------------------------
    bw_index: dict | None = None
    if not skip_rna:
        bw_index = build_bw_index(bw_cache)
        if not bw_index:
            logger.warning("No BigWig index entries found; RNA support will be skipped per locus")

    # ------------------------------------------------------------------
    # Discover species directories
    # ------------------------------------------------------------------
    if not nhmmer_dir.exists():
        logger.error(f"nHMMER directory not found: {nhmmer_dir}")
        raise typer.Exit(1)

    species_dirs = sorted([d for d in nhmmer_dir.iterdir() if d.is_dir()])
    logger.info(f"Found {len(species_dirs)} nHMMER species directories in {nhmmer_dir}")

    if resume and not dry_run:
        skip_all: list[Path] = []
        gff_only: list[Path] = []
        to_process: list[Path] = []
        for d in species_dirs:
            tsv = output / d.name / f"{d.name}.tsv"
            gff = output / d.name / f"{d.name}.gff"
            if tsv.exists() and gff.exists():
                skip_all.append(d)
            elif tsv.exists():
                gff_only.append(d)
            else:
                to_process.append(d)
        if skip_all:
            logger.info(f"Resume mode: skipping {len(skip_all)} already-completed species")
        if gff_only:
            logger.info(
                f"Resume mode: regenerating GFF for {len(gff_only)} species (TSV exists, GFF missing)"
            )
            for d in gff_only:
                tsv_path = output / d.name / f"{d.name}.tsv"
                gff_path = output / d.name / f"{d.name}.gff"
                write_gff_from_tsv(tsv_path, gff_path)
                logger.info(f"[{d.name}] Wrote {gff_path.name} from existing TSV")
        species_dirs = to_process

    # ------------------------------------------------------------------
    # Dispatch per-species workers (sequential or threaded)
    # ------------------------------------------------------------------
    worker_kwargs: dict = dict(
        miniprot_dir=miniprot_dir,
        e_value_threshold=e_value_threshold,
        e_value_high=e_value_high,
        min_span_bp=min_span_bp,
        max_span_bp=max_span_bp,
        positive_high=positive_high,
        type_mapping=type_mapping,
        specificity_pairs=specificity_pairs,
        type_to_family=type_to_family,
        anno_new=anno_new,
        anno_old=anno_old,
        bw_cache=bw_cache,
        bw_index=bw_index,
        skip_anno=skip_anno,
        skip_codon=skip_codon,
        skip_rna=skip_rna,
        use_llm=use_llm,
    )

    all_loci: list[dict] = []

    if workers == 1:
        for sp_dir in species_dirs:
            all_loci.extend(_process_species(sp_dir.name, sp_dir, **worker_kwargs))
    else:
        from concurrent.futures import ThreadPoolExecutor, as_completed

        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(_process_species, sp_dir.name, sp_dir, **worker_kwargs): sp_dir.name
                for sp_dir in species_dirs
            }
            for fut in as_completed(futures):
                species_id = futures[fut]
                try:
                    all_loci.extend(fut.result())
                except Exception as exc:
                    logger.error(f"Species {species_id} failed: {exc}")

    logger.info(
        f"Assembled {len(all_loci)} total locus candidates "
        f"({sum(1 for loc in all_loci if loc['has_ntd'] and loc['has_ctd'])} Full/potential, "
        f"{sum(1 for loc in all_loci if loc['has_ntd'] and not loc['has_ctd'])} NTD-only, "
        f"{sum(1 for loc in all_loci if not loc['has_ntd'] and loc['has_ctd'])} CTD-only)"
    )

    # ------------------------------------------------------------------
    # Post-processing: abbreviation map + Spidroin_ID assignment
    # (requires all species together; must run after workers complete)
    # ------------------------------------------------------------------
    all_species_ids = list({loc["species"] for loc in all_loci})
    abbrev_map = build_abbrev_map(all_species_ids)

    # Sort for deterministic Spidroin_ID assignment regardless of worker completion order
    all_loci.sort(key=lambda loc: (loc["species"], loc["seqid"], loc["start"]))

    species_counter: dict[str, int] = {}
    for locus in all_loci:
        sid = locus["species"]
        n = species_counter.get(sid, 0) + 1
        species_counter[sid] = n
        abbrev = abbrev_map.get(sid, sid[:4].replace(".", ""))
        locus["Spidroin_ID"] = f"{abbrev}_spid_{n:05d}"

    # ------------------------------------------------------------------
    # Build output DataFrame (Feishu-aligned columns first, then diagnostics)
    # ------------------------------------------------------------------

    # completeness → Hint_type mapping ("Repeat_only" → "Internal" per Feishu schema)
    _hint_map = {
        "Full_length": "Full_length",
        "N-terminal": "N-terminal",
        "C-terminal": "C-terminal",
        "Repeat_only": "Internal",
    }

    for locus in all_loci:
        hint = _hint_map.get(locus.get("completeness", ""), locus.get("completeness", ""))
        locus["Hint_type"] = hint
        locus["Full_length"] = hint == "Full_length"
        # Species display: replace underscores with spaces in the name part, keep NNN prefix
        sp = locus["species"]
        nnn, name = sp.split(".", 1) if "." in sp else ("", sp)
        locus["Species_display"] = f"{nnn}.{name.replace('_', ' ')}" if nnn else name.replace("_", " ")
        locus["Length"] = locus["end"] - locus["start"] + 1
        # Note column: trim to avoid empty trailing semicolons
        locus["Note"] = (locus.get("notes") or "").strip("; ")

    # Feishu-aligned columns (match docs/feishu_column.md)
    feishu_cols = [
        ("Spidroin_ID", "Spidroin_ID"),
        ("Species_display", "Species"),
        ("seqid", "Chr"),
        ("start", "Start"),
        ("end", "End"),
        ("strand", "Stran"),
        ("spidroin_type", "Spidroin_type"),
        ("Full_length", "Full_length"),
        ("Length", "Length"),
        ("Hint_type", "Hint_type"),
        ("Note", "Note"),
    ]
    # Diagnostic columns (for internal review; not imported to Feishu)
    diag_cols = [
        "ntd_profile", "ntd_e_value", "ntd_start", "ntd_end",
        "ctd_profile", "ctd_e_value", "ctd_start", "ctd_end",
        "miniprot_ref", "miniprot_positive", "miniprot_identity", "miniprot_rank",
        "anno_gene_id", "anno_overlap_frac",
        "bgi_rna_support", "ont_rna_support", "rna_mean_signal",
        "ntd_start_codon", "ctd_stop_codon",
        "confidence", "needs_review",
    ]

    # Ensure all expected keys exist
    for locus in all_loci:
        for _, col in feishu_cols:
            locus.setdefault(col, None)
        for col in diag_cols:
            locus.setdefault(col, None)

    # Group loci by species; all_loci is already sorted by (species, seqid, start)
    from itertools import groupby as _groupby

    confidence_order = {"high": 0, "medium": 1, "low": 2}

    total_written = 0
    total_high = total_medium = total_low = total_review = 0

    for species_id, sp_loci_iter in _groupby(all_loci, key=lambda loc: loc["species"]):
        sp_loci = list(sp_loci_iter)

        # Build rows for this species
        rows_out = []
        for locus in sp_loci:
            row: dict = {}
            for src_key, dst_col in feishu_cols:
                row[dst_col] = locus.get(src_key)
            for col in diag_cols:
                row[col] = locus.get(col)
            rows_out.append(row)

        df = pl.DataFrame(rows_out)

        # Sort within species: confidence (high first), Chr, Start
        df = df.with_columns(
            pl.col("confidence").replace(confidence_order, default=3).alias("_conf_order")
        ).sort(["_conf_order", "Chr", "Start"]).drop("_conf_order")

        high = (df["confidence"] == "high").sum()
        medium = (df["confidence"] == "medium").sum()
        low = (df["confidence"] == "low").sum()
        review = df["needs_review"].sum() if df["needs_review"].dtype == pl.Boolean else 0

        total_written += len(df)
        total_high += high
        total_medium += medium
        total_low += low
        total_review += review

        if dry_run:
            logger.info(
                f"[{species_id}] {len(df)} loci | high={high} medium={medium} low={low} | needs_review={review}"
            )
        else:
            sp_out = output / species_id
            sp_out.mkdir(parents=True, exist_ok=True)
            out_file = sp_out / f"{species_id}.tsv"
            df.write_csv(out_file, separator="\t", null_value="")
            gff_file = sp_out / f"{species_id}.gff"
            write_loci_gff(sp_loci, species_id, gff_file)
            fasta_root_out = anno_old if (anno_old and anno_old.exists()) else None
            fasta_root_new_out = anno_new if (anno_new and anno_new.exists()) else None
            write_fasta_and_hints(sp_loci, species_id, sp_out, fasta_root_out, fasta_root_new_out)
            logger.info(
                f"[{species_id}] Wrote {len(df)} rows → {sp_out.name}/{out_file.name} + {gff_file.name} "
                f"(high={high} medium={medium} low={low})"
            )

    logger.info(
        f"Results: {total_written} total loci | "
        f"high={total_high} medium={total_medium} low={total_low} | needs_review={total_review}"
    )
    if dry_run:
        logger.info("Dry-run mode: no files written")
    else:
        logger.success(f"Wrote {total_written} rows across {len(all_species_ids)} species to {output}")


@app.command("generate-gff")
def generate_gff_cmd(
    output: Path = typer.Option(
        _DEFAULT_OUTPUT,
        help="Output directory containing per-species subdirectories.",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        help="Overwrite existing GFF files.",
    ),
) -> None:
    """
    Generate (or regenerate) GFF files from existing TSV results without re-running the pipeline.
    Reads {output}/{species}/{species}.tsv and writes {output}/{species}/{species}.gff.
    """
    if not output.is_dir():
        logger.error(f"Output directory not found: {output}")
        raise typer.Exit(1)

    species_dirs = sorted([d for d in output.iterdir() if d.is_dir()])
    if not species_dirs:
        logger.warning(f"No species subdirectories found in {output}")
        raise typer.Exit(0)

    generated = skipped = 0
    for sp_dir in species_dirs:
        tsv_path = sp_dir / f"{sp_dir.name}.tsv"
        gff_path = sp_dir / f"{sp_dir.name}.gff"
        if not tsv_path.exists():
            logger.debug(f"[{sp_dir.name}] No TSV found, skipping")
            continue
        if gff_path.exists() and not force:
            logger.debug(f"[{sp_dir.name}] GFF already exists, skipping (use --force to overwrite)")
            skipped += 1
            continue
        write_gff_from_tsv(tsv_path, gff_path)
        logger.info(f"[{sp_dir.name}] Wrote {gff_path.name}")
        generated += 1

    logger.success(f"Generated {generated} GFF files, skipped {skipped} (already present)")


if __name__ == "__main__":
    app()
