#!/usr/bin/env python3
# filename: miniprot_dedup.py
# Purpose: Deduplicate miniprot GFF output:
# - Parse miniprot mRNA GFF output
# - Filter by positive (homology) threshold
# - Deduplicate overlapping hits by (seqid, gene_type, terminal), keeping best score
# - Output: deduplicated hits TSV
#
# Usage:
#   python miniprot_dedup.py \
#       --miniprot-gff miniprot_out.mRNA.gff \
#       --output dedup_hits.tsv \
#       --min-positive 0.85
#
# Notes:
# - Input: miniprot GFF output with mRNA records
# - Target field format: id|num|Family|Genus|species|GeneType|GeneType|Terminal start end
# - Terminal: NTD (N-terminal) or CTD (C-terminal)
# - Deduplication: overlapping hits with same (seqid, gene_type, terminal) keep best score

import argparse
import csv
from collections import defaultdict
from typing import Dict, List

from natsort import natsort_keygen

# Natural sort key for chromosome names (Chr1, Chr2, ..., Chr10, Chr11)
natural_sort_key = natsort_keygen()


def parse_gff_attributes(attr_str: str) -> Dict[str, str]:
    """Parse GFF9 attributes field"""
    attrs = {}
    for item in attr_str.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value
    return attrs


def parse_target_field(target_str: str) -> Dict[str, str]:
    """
    Parse Target field: id|num|Family|Genus|species|GeneType|GeneType|Terminal start end
    Example: 2600|5931|Linyphiidae|Lepthyphantes|minutus|PySp|PySp|CTD 1 225
    """
    parts = target_str.split("|")
    if len(parts) < 8:
        return {}

    # Last part contains "Terminal start end"
    last_part = parts[7].split()
    terminal = last_part[0] if last_part else ""

    return {
        "target_id": parts[0],
        "target_num": parts[1],
        "family": parts[2],
        "genus": parts[3],
        "species": parts[4],
        "gene_type": parts[5],  # e.g., MaSp, Flag, MiSp, PySp, AcSp
        "gene_type2": parts[6],
        "terminal": terminal,  # NTD or CTD
    }


def read_miniprot_gff(path: str) -> List[Dict]:
    """Read miniprot GFF output and extract mRNA records"""
    rows = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) < 9:
                continue

            # Only process mRNA records
            if fields[2] != "mRNA":
                continue

            seqid = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            score = float(fields[5])
            strand = fields[6]
            attrs = parse_gff_attributes(fields[8])

            target_info = {}
            if "Target" in attrs:
                target_info = parse_target_field(attrs["Target"])

            row = {
                "seqid": seqid,
                "start": min(start, end),
                "end": max(start, end),
                "score": score,
                "strand": strand,
                "mp_id": attrs.get("ID", ""),
                "rank": int(attrs.get("Rank", 1)),
                "identity": float(attrs.get("Identity", 0)),
                "positive": float(attrs.get("Positive", 0)),
                **target_info,
            }
            rows.append(row)
    return rows


def intervals_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    """Check if two intervals overlap"""
    return a_start <= b_end and b_start <= a_end


def select_best_hits(hits: List[Dict]) -> List[Dict]:
    """
    Deduplicate hits by overlap detection, keeping best hit per overlapping region.
    All overlapping hits on the same seqid are deduplicated regardless of gene_type/terminal.
    """
    # Sort by seqid, then score desc, identity desc
    hits_sorted = sorted(
        hits,
        key=lambda x: (
            x["seqid"],
            -x["score"],
            -x["identity"],
        ),
    )

    # Group by seqid only for overlap detection
    grouped: Dict[str, List[Dict]] = defaultdict(list)
    for h in hits_sorted:
        grouped[h["seqid"]].append(h)

    chosen = []
    for seqid, group in grouped.items():
        # Greedy: pick best score first, skip any that overlap with already chosen
        used = []
        for h in group:  # Already sorted by score desc
            overlaps = False
            for chosen_hit in used:
                if intervals_overlap(h["start"], h["end"], chosen_hit["start"], chosen_hit["end"]):
                    overlaps = True
                    break
            if not overlaps:
                chosen.append(h)
                used.append(h)

    return chosen


def write_dedup_report(hits: List[Dict], out_path: str):
    """Write deduplicated hits report"""
    cols = [
        "seqid",
        "start",
        "end",
        "score",
        "identity",
        "positive",
        "strand",
        "gene_type",
        "terminal",
        "family",
        "genus",
        "species",
        "mp_id",
    ]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for h in sorted(hits, key=lambda x: (natural_sort_key(x["seqid"]), x["start"])):
            w.writerow(h)


def main():
    ap = argparse.ArgumentParser(
        description="Deduplicate miniprot GFF output for spidroin annotation."
    )
    ap.add_argument("--miniprot-gff", required=True, help="Miniprot mRNA GFF output file.")
    ap.add_argument("--output", required=True, help="Output: deduplicated hits TSV.")
    ap.add_argument(
        "--min-positive",
        type=float,
        default=0.85,
        help="Minimum positive (homology) threshold for filtering (default 0.85).",
    )
    args = ap.parse_args()

    # Read miniprot GFF
    print(f"Reading miniprot GFF: {args.miniprot_gff}")
    hits = read_miniprot_gff(args.miniprot_gff)
    print(f"  Found {len(hits)} mRNA records")

    # Filter by positive threshold
    hits = [h for h in hits if h.get("positive", 0) >= args.min_positive]
    print(f"  {len(hits)} hits after positive >= {args.min_positive} filter")

    # Deduplicate
    print("Deduplicating overlapping hits...")
    dedup_hits = select_best_hits(hits)
    print(f"  {len(dedup_hits)} hits after deduplication")

    # Write output
    write_dedup_report(dedup_hits, args.output)
    print(f"Wrote deduplicated hits to: {args.output}")


if __name__ == "__main__":
    main()
