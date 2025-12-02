#!/usr/bin/env python3
# filename: miniprot_pair.py
# Purpose: Pair N/C terminal hits from deduplicated miniprot results:
# - Read deduplicated hits TSV (from miniprot_dedup.py)
# - Pair N-terminal and C-terminal hits into genes (15-90kb apart, same gene type)
# - Output: paired genes with padding, unpaired hits without padding
#
# Usage:
#   python miniprot_pair.py \
#       --input dedup_hits.tsv \
#       --out-paired paired_genes.bed \
#       --out-unpaired unpaired_hits.bed \
#       --pad 1000 \
#       --min-distance 1000 \
#       --max-distance 100000 \
#       --ref-fasta /path/to/reference.fa.gz
#
# Notes:
# - Input: deduplicated hits TSV from miniprot_dedup.py
# - Pairing criteria:
#   - Same chromosome/scaffold
#   - Same gene type (MaSp, Flag, MiSp, etc.)
#   - N/C terminals 15,000-90,000 bp apart
#   - Prefer closer distance when multiple candidates
# - Chromosome sizes: auto-detect .chrom.sizes file, generate if not exists

import argparse
import csv
import gzip
import os
from collections import defaultdict
from typing import List, Dict, Tuple

from natsort import natsort_keygen

# Natural sort key for chromosome names (Chr1, Chr2, ..., Chr10, Chr11)
natural_sort_key = natsort_keygen()

try:
    from Bio import SeqIO
except ImportError:
    SeqIO = None


def read_dedup_tsv(path: str) -> List[Dict]:
    """Read deduplicated hits TSV file"""
    rows = []
    with open(path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            # Convert numeric fields
            row["start"] = int(row["start"])
            row["end"] = int(row["end"])
            row["score"] = float(row["score"])
            row["identity"] = float(row["identity"])
            row["positive"] = float(row["positive"])
            rows.append(row)
    return rows


def load_chrom_sizes(chrom_sizes_path: str) -> Dict[str, int]:
    """Load chromosome sizes from .chrom.sizes file (tsv: name\\tlength)"""
    contig_len = {}
    with open(chrom_sizes_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                contig_len[parts[0]] = int(parts[1])
    return contig_len


def build_contig_lengths(ref_fasta: str) -> Dict[str, int]:
    """Build contig length dictionary from reference FASTA (supports .gz)"""
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


def save_chrom_sizes(contig_len: Dict[str, int], out_path: str):
    """Save chromosome sizes to .chrom.sizes file"""
    with open(out_path, "w") as f:
        for name, length in sorted(contig_len.items()):
            f.write(f"{name}\t{length}\n")


def pair_nc_terminals(
    hits: List[Dict], min_distance: int, max_distance: int
) -> Tuple[List[Dict], List[Dict]]:
    """
    Pair N-terminal and C-terminal hits into genes.

    Pairing criteria:
    - Same seqid (chromosome/scaffold)
    - N/C terms are min_distance - max_distance bp apart
    - Relaxed gene_type matching (only by chromosome, not gene_type)
    - Prefer closer distance when multiple candidates exist

    Returns (paired_genes, unpaired_hits)
    """
    # Separate N-terminal and C-terminal hits
    ntd_hits = [h for h in hits if h.get("terminal") == "NTD"]
    ctd_hits = [h for h in hits if h.get("terminal") == "CTD"]
    other_hits = [h for h in hits if h.get("terminal") not in ("NTD", "CTD")]

    paired_genes = []
    used_ntd = set()
    used_ctd = set()

    # Group by seqid only (relaxed matching - ignore gene_type)
    ntd_by_chr = defaultdict(list)
    ctd_by_chr = defaultdict(list)

    for h in ntd_hits:
        ntd_by_chr[h["seqid"]].append(h)

    for h in ctd_hits:
        ctd_by_chr[h["seqid"]].append(h)

    # Sort hits within each group by position
    for seqid in ntd_by_chr:
        ntd_by_chr[seqid].sort(key=lambda x: x["start"])
    for seqid in ctd_by_chr:
        ctd_by_chr[seqid].sort(key=lambda x: x["start"])

    # Try to pair NTD with CTD (by chromosome only)
    for seqid, ntd_list in ntd_by_chr.items():
        ctd_list = ctd_by_chr.get(seqid, [])

        for ntd in ntd_list:
            if id(ntd) in used_ntd:
                continue

            best_ctd = None
            best_distance = float("inf")

            for ctd in ctd_list:
                if id(ctd) in used_ctd:
                    continue

                # Calculate distance as span from first start to last end
                span_start = min(ntd["start"], ctd["start"])
                span_end = max(ntd["end"], ctd["end"])
                distance = span_end - span_start

                # Skip if overlapping (one contains another)
                if ntd["start"] <= ctd["start"] <= ctd["end"] <= ntd["end"]:
                    continue
                if ctd["start"] <= ntd["start"] <= ntd["end"] <= ctd["end"]:
                    continue

                # Check distance criteria, prefer closer distance
                if min_distance <= distance <= max_distance:
                    if distance < best_distance:
                        best_distance = distance
                        best_ctd = ctd

            if best_ctd:
                # Create paired gene
                gene_start = min(ntd["start"], best_ctd["start"])
                gene_end = max(ntd["end"], best_ctd["end"])

                # Determine gene_type by higher positive value
                ntd_type = ntd.get("gene_type", "")
                ctd_type = best_ctd.get("gene_type", "")
                ntd_positive = ntd.get("positive", 0)
                ctd_positive = best_ctd.get("positive", 0)
                gene_type = ntd_type if ntd_positive >= ctd_positive else ctd_type

                paired_genes.append(
                    {
                        "seqid": seqid,
                        "start": gene_start,
                        "end": gene_end,
                        "gene_type": gene_type,
                        "ntd_type": ntd_type,
                        "ctd_type": ctd_type,
                        "strand": ntd["strand"],
                        "ntd_hit": ntd,
                        "ctd_hit": best_ctd,
                        "ntd_id": ntd["mp_id"],
                        "ctd_id": best_ctd["mp_id"],
                        "ntd_species": f"{ntd.get('genus', '')} {ntd.get('species', '')}",
                        "ctd_species": f"{best_ctd.get('genus', '')} {best_ctd.get('species', '')}",
                        "ntd_positive": ntd_positive,
                        "ctd_positive": ctd_positive,
                        "distance": best_distance,
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


def write_paired_bed(
    paired_genes: List[Dict],
    pad: int,
    contig_lengths: Dict[str, int],
    out_path: str,
):
    """Write paired genes as BED file with padding"""
    with open(out_path, "w") as f:
        f.write(
            "#seqid\tstart\tend\tgene_name\tscore\tstrand\tntd_id\tctd_id\tlength\t"
            "ntd_type\tctd_type\tntd_species\tctd_species\n"
        )
        for g in sorted(paired_genes, key=lambda x: (natural_sort_key(x["seqid"]), x["start"])):
            padded_start, padded_end = apply_padding_and_clamp(
                g["start"], g["end"], pad, g["seqid"], contig_lengths
            )
            gene_name = g['gene_type']
            score = 0
            f.write(
                f"{g['seqid']}\t{padded_start}\t{padded_end}\t{gene_name}\t{score}\t{g['strand']}\t"
                f"{g['ntd_id']}\t{g['ctd_id']}\t{g['distance']}\t"
                f"{g['ntd_type']}\t{g['ctd_type']}\t{g['ntd_species']}\t{g['ctd_species']}\n"
            )


def write_unpaired_bed(hits: List[Dict], out_path: str):
    """Write unpaired hits as BED file without padding"""
    with open(out_path, "w") as f:
        f.write("#seqid\tstart\tend\tgene_name\tscore\tstrand\tmp_id\tterminal\tspecies\n")
        for h in sorted(hits, key=lambda x: (natural_sort_key(x["seqid"]), x["start"])):
            gene_name = f"{h.get('gene_type', 'unknown')}_{h.get('terminal', 'unknown')}"
            species = f"{h.get('genus', '')} {h.get('species', '')}"
            f.write(
                f"{h['seqid']}\t{h['start']}\t{h['end']}\t{gene_name}\t{h['score']}\t{h['strand']}\t"
                f"{h['mp_id']}\t{h.get('terminal', '')}\t{species}\n"
            )


def main():
    ap = argparse.ArgumentParser(
        description="Pair N/C terminal hits from deduplicated miniprot results."
    )
    ap.add_argument(
        "--input", required=True, help="Input: deduplicated hits TSV from miniprot_dedup.py."
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
        "--pad",
        type=int,
        default=1000,
        help="Up/downstream padding for paired genes (default 1000).",
    )
    ap.add_argument(
        "--min-distance",
        type=int,
        default=1000,
        help="Minimum N/C terminal distance for pairing (default 15000).",
    )
    ap.add_argument(
        "--max-distance",
        type=int,
        default=100000,
        help="Maximum N/C terminal distance for pairing (default 100000).",
    )
    ap.add_argument(
        "--ref-fasta", default=None, help="Reference FASTA (to clamp coordinates)."
    )
    ap.add_argument(
        "--chrom-sizes",
        default=None,
        help="Chromosome sizes file (tsv: name\\tlength). Faster than --ref-fasta.",
    )
    args = ap.parse_args()

    # Read deduplicated hits
    print(f"Reading deduplicated hits: {args.input}")
    hits = read_dedup_tsv(args.input)
    print(f"  Found {len(hits)} hits")

    # Build contig lengths (check chrom-sizes first, then ref-fasta)
    contig_lengths = {}
    chrom_sizes_path = args.chrom_sizes
    if args.ref_fasta and not chrom_sizes_path:
        # Auto-derive chrom.sizes path from ref-fasta
        chrom_sizes_path = args.ref_fasta.replace(".fa.gz", "").replace(".fa", "") + ".chrom.sizes"

    if chrom_sizes_path and os.path.exists(chrom_sizes_path):
        print(f"Loading chromosome sizes from: {chrom_sizes_path}")
        contig_lengths = load_chrom_sizes(chrom_sizes_path)
    elif args.ref_fasta:
        print(f"Building chromosome sizes from: {args.ref_fasta}")
        contig_lengths = build_contig_lengths(args.ref_fasta)
        # Auto-save chrom.sizes for future use
        if chrom_sizes_path:
            save_chrom_sizes(contig_lengths, chrom_sizes_path)
            print(f"  Saved chromosome sizes to: {chrom_sizes_path}")

    # Pair N/C terminals
    print(
        f"Pairing N/C terminals (distance: {args.min_distance}-{args.max_distance} bp)..."
    )
    paired_genes, unpaired_hits = pair_nc_terminals(
        hits, args.min_distance, args.max_distance
    )
    print(f"  {len(paired_genes)} paired genes")
    print(f"  {len(unpaired_hits)} unpaired hits")

    # Write outputs
    write_paired_bed(paired_genes, args.pad, contig_lengths, args.out_paired)
    print(f"Wrote paired genes (with {args.pad}bp padding) to: {args.out_paired}")

    write_unpaired_bed(unpaired_hits, args.out_unpaired)
    print(f"Wrote unpaired hits (no padding) to: {args.out_unpaired}")


if __name__ == "__main__":
    main()
