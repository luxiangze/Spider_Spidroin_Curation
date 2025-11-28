#!/usr/bin/env python3
# filename: miniprot_to_coordinates.py
# Purpose: Process miniprot GFF output to extract spidroin gene coordinates:
# - Parse miniprot mRNA GFF output
# - Round start pos to nearest 500 bp for binning
# - Sort by seqid -> rounded_start -> Score desc
# - Deduplicate by (seqid, rounded_start) keeping best hit
# - Pair N-terminal and C-terminal hits into genes (15-90kb apart, same gene type)
# - Output: 1) deduplicated hits, 2) paired genes with padding, 3) unpaired hits without padding
#
# Usage:
#   python miniprot_to_coordinates.py \
#       --miniprot-gff miniprot_out.mRNA.gff \
#       --out-dedup dedup_hits.tsv \
#       --out-paired paired_genes.bed \
#       --out-unpaired unpaired_hits.bed \
#       --bin-size 500 \
#       --pad 1000 \
#       --min-distance 15000 \
#       --max-distance 90000 \
#       --ref-fasta /path/to/reference.fa
#       --min-positive 0.85
#
# Notes:
# - Input: miniprot GFF output with mRNA records
# - Target field format: id|num|Family|Genus|species|GeneType|GeneType|Terminal start end
# - Terminal: NTD (N-terminal) or CTD (C-terminal)
# - Pairing criteria:
#   - Same chromosome/scaffold
#   - N/C terminals 15,000-90,000 bp apart
#   - Same gene type (MaSp, Flag, MiSp, etc.)
#   - Prefer same/closely related species

import argparse
import csv
import gzip
from collections import defaultdict
from typing import List, Dict, Tuple, Optional

try:
    from Bio import SeqIO
except ImportError:
    SeqIO = None


def mround(value: int, base: int) -> int:
    """Excel MROUND equivalent: round to nearest multiple of base"""
    return int(base * round(value / float(base)))


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
        for line_num, line in enumerate(f, start=1):
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


def build_contig_lengths(ref_fasta: Optional[str]) -> Dict[str, int]:
    """Build contig length dictionary from reference FASTA (supports .gz)"""
    if not ref_fasta:
        return {}
    if SeqIO is None:
        raise RuntimeError(
            "Biopython not installed; install biopython or omit --ref-fasta."
        )
    contig_len = {}
    # Support gzip compressed files
    open_func = gzip.open if ref_fasta.endswith(".gz") else open
    mode = "rt" if ref_fasta.endswith(".gz") else "r"
    with open_func(ref_fasta, mode) as f:
        for rec in SeqIO.parse(f, "fasta"):
            contig_len[rec.id] = len(rec.seq)
    return contig_len


def select_best_hits(
    hits: List[Dict], bin_size: int
) -> Tuple[List[Dict], List[Dict]]:
    """
    Deduplicate hits by binned position, keeping best hit per bin.
    Returns (chosen_hits, report_rows)
    """
    for h in hits:
        h["rounded_start"] = mround(h["start"], bin_size)

    # Sort: seqid -> rounded_start -> score desc -> identity desc
    hits_sorted = sorted(
        hits,
        key=lambda x: (x["seqid"], x["rounded_start"], -x["score"], -x["identity"]),
    )

    chosen = []
    report_rows = []
    grouped: Dict[Tuple[str, int], List[Dict]] = defaultdict(list)
    for h in hits_sorted:
        grouped[(h["seqid"], h["rounded_start"])].append(h)

    for key, group in grouped.items():
        best = group[0]
        chosen.append(best)
        # Report top 5 candidates for this bin
        for rank, cand in enumerate(group[:5], start=1):
            report_rows.append(
                {
                    "seqid": cand["seqid"],
                    "rounded_bin": cand["rounded_start"],
                    "rank": rank,
                    "start": cand["start"],
                    "end": cand["end"],
                    "score": cand["score"],
                    "identity": cand["identity"],
                    "gene_type": cand.get("gene_type", ""),
                    "terminal": cand.get("terminal", ""),
                    "species": f"{cand.get('genus', '')} {cand.get('species', '')}",
                    "selected": "yes" if cand is best else "no",
                    "mp_id": cand["mp_id"],
                }
            )

    return chosen, report_rows


def pair_nc_terminals(
    hits: List[Dict], min_distance: int, max_distance: int
) -> Tuple[List[Dict], List[Dict]]:
    """
    Pair N-terminal and C-terminal hits into genes.

    Pairing criteria:
    - Same seqid (chromosome/scaffold)
    - N/C terms are min_distance - max_distance bp apart
    - Same gene type (e.g., MaSp, Flag, MiSp)
    - Same or closely related species (same family preferred)

    Returns (paired_genes, unpaired_hits)
    """
    # Separate N-terminal and C-terminal hits
    ntd_hits = [h for h in hits if h.get("terminal") == "NTD"]
    ctd_hits = [h for h in hits if h.get("terminal") == "CTD"]
    other_hits = [h for h in hits if h.get("terminal") not in ("NTD", "CTD")]

    paired_genes = []
    used_ntd = set()
    used_ctd = set()

    # Group by seqid and gene_type for efficient pairing
    ntd_by_chr_gene = defaultdict(list)
    ctd_by_chr_gene = defaultdict(list)

    for h in ntd_hits:
        key = (h["seqid"], h.get("gene_type", ""))
        ntd_by_chr_gene[key].append(h)

    for h in ctd_hits:
        key = (h["seqid"], h.get("gene_type", ""))
        ctd_by_chr_gene[key].append(h)

    # Sort hits within each group by position
    for key in ntd_by_chr_gene:
        ntd_by_chr_gene[key].sort(key=lambda x: x["start"])
    for key in ctd_by_chr_gene:
        ctd_by_chr_gene[key].sort(key=lambda x: x["start"])

    # Try to pair NTD with CTD
    for key, ntd_list in ntd_by_chr_gene.items():
        seqid, gene_type = key
        ctd_list = ctd_by_chr_gene.get(key, [])

        for ntd in ntd_list:
            if id(ntd) in used_ntd:
                continue

            best_ctd = None
            best_distance = float("inf")

            for ctd in ctd_list:
                if id(ctd) in used_ctd:
                    continue

                # Calculate distance between hits
                # Distance should be between the end of one and start of another
                if ntd["end"] < ctd["start"]:
                    distance = ctd["start"] - ntd["end"]
                elif ctd["end"] < ntd["start"]:
                    distance = ntd["start"] - ctd["end"]
                else:
                    # Overlapping - skip
                    continue

                # Check distance criteria
                if min_distance <= distance <= max_distance:
                    # Prefer same family/genus
                    same_family = ntd.get("family", "") == ctd.get("family", "")
                    same_genus = ntd.get("genus", "") == ctd.get("genus", "")

                    # Score: prefer closer distance and same taxonomy
                    score = distance
                    if same_genus:
                        score -= 10000
                    elif same_family:
                        score -= 5000

                    if score < best_distance:
                        best_distance = score
                        best_ctd = ctd

            if best_ctd:
                # Create paired gene
                gene_start = min(ntd["start"], best_ctd["start"])
                gene_end = max(ntd["end"], best_ctd["end"])

                paired_genes.append(
                    {
                        "seqid": seqid,
                        "start": gene_start,
                        "end": gene_end,
                        "gene_type": gene_type,
                        "strand": ntd["strand"],
                        "ntd_hit": ntd,
                        "ctd_hit": best_ctd,
                        "ntd_id": ntd["mp_id"],
                        "ctd_id": best_ctd["mp_id"],
                        "ntd_species": f"{ntd.get('genus', '')} {ntd.get('species', '')}",
                        "ctd_species": f"{best_ctd.get('genus', '')} {best_ctd.get('species', '')}",
                        "distance": abs(best_ctd["start"] - ntd["end"])
                        if ntd["end"] < best_ctd["start"]
                        else abs(ntd["start"] - best_ctd["end"]),
                    }
                )

                used_ntd.add(id(ntd))
                used_ctd.add(id(best_ctd))

    # Collect unpaired hits
    unpaired = []
    for h in ntd_hits:
        if id(h) not in used_ntd:
            unpaired.append(h)
    for h in ctd_hits:
        if id(h) not in used_ctd:
            unpaired.append(h)
    unpaired.extend(other_hits)

    return paired_genes, unpaired


def apply_padding_and_clamp(
    start: int, end: int, pad: int, seqid: str, contig_lengths: Dict[str, int]
) -> Tuple[int, int]:
    """Apply padding and clamp to contig boundaries"""
    padded_start = max(1, start - pad)
    padded_end = end + pad
    if seqid in contig_lengths:
        padded_end = min(padded_end, contig_lengths[seqid])
    return padded_start, padded_end


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
        w = csv.DictWriter(
            f, fieldnames=cols, delimiter="\t", extrasaction="ignore"
        )
        w.writeheader()
        for h in sorted(hits, key=lambda x: (x["seqid"], x["start"])):
            w.writerow(h)


def write_paired_bed(
    paired_genes: List[Dict],
    pad: int,
    contig_lengths: Dict[str, int],
    out_path: str,
):
    """Write paired genes as BED file with padding"""
    with open(out_path, "w") as f:
        f.write(
            "#seqid\tstart\tend\tgene_name\tscore\tstrand\tntd_id\tctd_id\tdistance\tntd_species\tctd_species\n"
        )
        for g in sorted(paired_genes, key=lambda x: (x["seqid"], x["start"])):
            padded_start, padded_end = apply_padding_and_clamp(
                g["start"], g["end"], pad, g["seqid"], contig_lengths
            )
            gene_name = f"{g['gene_type']}_paired"
            score = 0  # Could calculate combined score
            f.write(
                f"{g['seqid']}\t{padded_start}\t{padded_end}\t{gene_name}\t{score}\t{g['strand']}\t"
                f"{g['ntd_id']}\t{g['ctd_id']}\t{g['distance']}\t{g['ntd_species']}\t{g['ctd_species']}\n"
            )


def write_unpaired_bed(hits: List[Dict], out_path: str):
    """Write unpaired hits as BED file without padding"""
    with open(out_path, "w") as f:
        f.write("#seqid\tstart\tend\tgene_name\tscore\tstrand\tmp_id\tterminal\tspecies\n")
        for h in sorted(hits, key=lambda x: (x["seqid"], x["start"])):
            gene_name = f"{h.get('gene_type', 'unknown')}_{h.get('terminal', 'unknown')}"
            species = f"{h.get('genus', '')} {h.get('species', '')}"
            f.write(
                f"{h['seqid']}\t{h['start']}\t{h['end']}\t{gene_name}\t{h['score']}\t{h['strand']}\t"
                f"{h['mp_id']}\t{h.get('terminal', '')}\t{species}\n"
            )


def main():
    ap = argparse.ArgumentParser(
        description="Process miniprot GFF output for spidroin gene annotation."
    )
    ap.add_argument(
        "--miniprot-gff", required=True, help="Miniprot mRNA GFF output file."
    )
    ap.add_argument(
        "--out-dedup", required=True, help="Output: deduplicated hits TSV."
    )
    ap.add_argument(
        "--out-paired",
        required=True,
        help="Output: paired genes BED (with padding).",
    )
    ap.add_argument(
        "--out-unpaired",
        required=True,
        help="Output: unpaired hits BED (no padding).",
    )
    ap.add_argument(
        "--bin-size",
        type=int,
        default=500,
        help="MROUND bin size for deduplication (default 500).",
    )
    ap.add_argument(
        "--pad",
        type=int,
        default=1000,
        help="Up/downstream padding for paired genes (default 1000).",
    )
    ap.add_argument(
        "--min-distance",
        type=int,
        default=15000,
        help="Minimum N/C terminal distance for pairing (default 15000).",
    )
    ap.add_argument(
        "--max-distance",
        type=int,
        default=90000,
        help="Maximum N/C terminal distance for pairing (default 90000).",
    )
    ap.add_argument(
        "--ref-fasta", default=None, help="Reference FASTA (to clamp coordinates)."
    )
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

    # Build contig lengths if reference provided
    contig_lengths = (
        build_contig_lengths(args.ref_fasta) if args.ref_fasta else {}
    )

    # Deduplicate by binned position
    print(f"Deduplicating hits (bin size: {args.bin_size} bp)...")
    dedup_hits, _ = select_best_hits(hits, args.bin_size)
    print(f"  {len(dedup_hits)} hits after deduplication")

    # Pair N/C terminals
    print(
        f"Pairing N/C terminals (distance: {args.min_distance}-{args.max_distance} bp)..."
    )
    paired_genes, unpaired_hits = pair_nc_terminals(
        dedup_hits, args.min_distance, args.max_distance
    )
    print(f"  {len(paired_genes)} paired genes")
    print(f"  {len(unpaired_hits)} unpaired hits")

    # Write outputs
    write_dedup_report(dedup_hits, args.out_dedup)
    print(f"Wrote deduplicated hits to: {args.out_dedup}")

    write_paired_bed(paired_genes, args.pad, contig_lengths, args.out_paired)
    print(f"Wrote paired genes (with {args.pad}bp padding) to: {args.out_paired}")

    write_unpaired_bed(unpaired_hits, args.out_unpaired)
    print(f"Wrote unpaired hits (no padding) to: {args.out_unpaired}")


if __name__ == "__main__":
    main()
