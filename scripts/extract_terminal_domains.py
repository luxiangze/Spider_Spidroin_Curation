#!/usr/bin/env python3
"""
Extract non-repetitive CTD/NTD terminal domain sequences from spidroin sequences.

Algorithm:
1. Auto-detect repeat unit length using autocorrelation analysis
2. Find repeat region boundary based on detected unit length
3. Extract terminal domain (non-repetitive region)
4. Fallback to triplet-based detection if no repeat unit found

Output:
- Sequences are grouped by spidroin type and domain type
- Output files: {SpidroinType}_{DomainType}.faa (e.g., MaSp1_CTD.faa, MiSp_NTD.faa)
- TuSp and CySp are merged into TuSp_CySp
- A processing report (processing_report.tsv) is generated

Usage:
    python extract_terminal_domains.py input.faa -o output_dir [options]

Options:
    -o, --output DIR       Output directory (required)
    -s, --skip N           Skip first N sequences (default: 0)
    -l, --min-length N     Minimum extracted length (default: 50)
    --similarity F         Repeat unit similarity threshold (default: 0.5)
    --no-simplify          Keep original headers

Example:
    python extract_terminal_domains.py seeds.faa -o output/ -s 19

Requires:
    pip install biopython

Author: Auto-generated
Date: 2025-11-30
"""

import argparse
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Spidroin type mapping (TuSp and CySp are the same)
SPIDROIN_TYPE_MAP = {
    'MASP': 'MaSp',
    'MASP1': 'MaSp1',
    'MASP2': 'MaSp2',
    'MASP3': 'MaSp3',
    'MASP3B': 'MaSp3B',
    'MISP': 'MiSp',
    'ACSP': 'AcSp',
    'ACSP1': 'AcSp1',
    'PYSP': 'PySp',
    'TUSP': 'TuSp_CySp',
    'CYSP': 'TuSp_CySp',
    'FLAG': 'Flag',
    'PFLAG': 'pFlag',
    'AGSP': 'AgSp',
    'AGSP1': 'AgSp1',
    'AGSP2': 'AgSp2',
    'CRSP': 'CrSp',
    'SPIDROIN': 'Other',
}


def get_spidroin_type(header: str) -> str:
    """Extract spidroin type from header and normalize it."""
    header_upper = header.upper()

    # Remove domain suffix
    for suffix in ['_CTD', '_NTD']:
        if suffix in header_upper:
            header_upper = header_upper.replace(suffix, '')

    # Try to find spidroin type in header
    parts = header_upper.replace('-', '_').split('_')

    for part in reversed(parts):  # Check from end
        if part in SPIDROIN_TYPE_MAP:
            return SPIDROIN_TYPE_MAP[part]

    # Check for partial matches
    for key, value in SPIDROIN_TYPE_MAP.items():
        if key in header_upper:
            return value

    return 'Other'


def get_domain_type(header: str) -> str:
    """Get domain type from header."""
    header_upper = header.upper()
    if header_upper.endswith('_NTD'):
        return 'NTD'
    elif header_upper.endswith('_CTD'):
        return 'CTD'
    elif '_NTD_' in header_upper:
        return 'NTD'
    elif '_CTD_' in header_upper:
        return 'CTD'
    return 'UNKNOWN'


def find_repeat_unit_length(seq: str, min_len: int = 15, max_len: int = 300,
                             min_similarity: float = 0.4) -> tuple[int | None, float]:
    """
    Auto-detect repeat unit length using autocorrelation.

    Args:
        seq: Input amino acid sequence
        min_len: Minimum repeat unit length to search
        max_len: Maximum repeat unit length to search
        min_similarity: Minimum similarity threshold to consider as repeat

    Returns:
        Tuple of (unit_length, similarity_score) or (None, 0) if not found
    """
    seq_len = len(seq)
    if seq_len < min_len * 2:
        return None, 0

    best_lag = None
    best_score = 0

    for lag in range(min_len, min(max_len, seq_len // 2)):
        overlap_len = seq_len - lag
        if overlap_len < lag:
            break

        # Count character matches between seq[0:n-lag] and seq[lag:n]
        matches = sum(1 for i in range(overlap_len) if seq[i] == seq[lag + i])
        score = matches / overlap_len

        if score >= min_similarity and score > best_score:
            best_score = score
            best_lag = lag

    return best_lag, best_score


def find_repeat_boundary_ctd(seq: str, unit_len: int | None,
                              similarity: float = 0.5) -> int:
    """
    Find where repeat region ends (for CTD extraction).

    Args:
        seq: Input amino acid sequence
        unit_len: Detected repeat unit length (None for fallback)
        similarity: Similarity threshold for repeat detection

    Returns:
        Position where repeat region ends (CTD starts after this)
    """
    seq_len = len(seq)

    # Fallback: use triplet-based detection if no unit_len
    if unit_len is None:
        window = 30
        for i in range(seq_len - window, -1, -1):
            triplets = Counter(seq[j:j+3] for j in range(i, min(i + window, seq_len - 2)))
            if any(c >= 3 for c in triplets.values()):
                return i + window
        return 0

    # Find the last position that is part of a repeat
    repeat_end = 0
    i = 0
    while i + unit_len <= seq_len:
        if i + unit_len * 2 <= seq_len:
            unit = seq[i:i+unit_len]
            next_unit = seq[i+unit_len:i+unit_len*2]
            matches = sum(1 for a, b in zip(unit, next_unit) if a == b)

            if matches / unit_len >= similarity:
                repeat_end = i + unit_len * 2
                i += unit_len
                continue
        i += 1

    return repeat_end


def find_repeat_boundary_ntd(seq: str, unit_len: int | None,
                              similarity: float = 0.5) -> int:
    """
    Find where repeat region starts (for NTD extraction).

    Args:
        seq: Input amino acid sequence
        unit_len: Detected repeat unit length (None for fallback)
        similarity: Similarity threshold for repeat detection

    Returns:
        Position where repeat region starts (NTD ends before this)
    """
    seq_len = len(seq)

    # Fallback: use triplet-based detection if no unit_len
    if unit_len is None:
        window = 30
        for i in range(seq_len - window + 1):
            triplets = Counter(seq[j:j+3] for j in range(i, min(i + window, seq_len - 2)))
            if any(c >= 3 for c in triplets.values()):
                return i
        return seq_len

    # Find first position that is part of a repeat
    for i in range(seq_len - unit_len * 2 + 1):
        unit = seq[i:i+unit_len]
        next_unit = seq[i+unit_len:i+unit_len*2]
        matches = sum(1 for a, b in zip(unit, next_unit) if a == b)
        if matches / unit_len >= similarity:
            return i

    return seq_len


def extract_terminal_domain(seq: str, domain_type: str,
                            similarity: float = 0.5) -> tuple[str, int | None]:
    """
    Extract terminal domain from sequence.

    Args:
        seq: Input amino acid sequence
        domain_type: 'CTD' or 'NTD'
        similarity: Similarity threshold for repeat detection

    Returns:
        Tuple of (extracted_sequence, detected_unit_length)
    """
    # Auto-detect repeat unit length
    unit_len, score = find_repeat_unit_length(seq)

    if domain_type == 'CTD':
        boundary = find_repeat_boundary_ctd(seq, unit_len, similarity)
        extracted = seq[boundary:]
    elif domain_type == 'NTD':
        boundary = find_repeat_boundary_ntd(seq, unit_len, similarity)
        extracted = seq[:boundary]
    else:
        extracted = seq

    return extracted, unit_len


def simplify_header(header: str) -> str:
    """
    Simplify header to format: number_type_domain
    Example: 20_5532_Araneidae_Neoscona_sp-OTU0257_MiSp_MiSp_CTD -> 20_MiSp_CTD
    """
    parts = header.split('_')
    if len(parts) < 3:
        return header

    seq_num = parts[0]
    domain_type = parts[-1]
    protein_type = parts[-2] if len(parts) > 2 else 'Spidroin'

    return f"{seq_num}_{protein_type}_{domain_type}"


def process_record(
    record: SeqRecord,
    similarity: float,
    simplify: bool
) -> tuple[SeqRecord, dict]:
    """
    Process a single SeqRecord and extract non-repetitive terminal domain.

    Args:
        record: BioPython SeqRecord object
        similarity: Similarity threshold for repeat detection
        simplify: Whether to simplify headers

    Returns:
        Tuple of (new SeqRecord, processing info dict)
    """
    header = record.description if record.description else record.id
    seq = str(record.seq)
    domain_type = get_domain_type(header)
    spidroin_type = get_spidroin_type(header)

    info = {
        'original_header': header,
        'original_length': len(seq),
        'domain_type': domain_type,
        'spidroin_type': spidroin_type,
        'extracted_length': len(seq),
        'unit_length': None,
        'status': 'unchanged',
    }

    if domain_type in ('CTD', 'NTD'):
        extracted, unit_len = extract_terminal_domain(seq, domain_type, similarity)
        info['status'] = 'processed'
        info['unit_length'] = unit_len
    else:
        info['status'] = 'unknown_domain'
        extracted = seq

    info['extracted_length'] = len(extracted)
    new_id = simplify_header(header) if simplify else header
    info['new_header'] = new_id

    new_record = SeqRecord(
        Seq(extracted),
        id=new_id,
        description=''
    )

    return new_record, info


def write_report(report_path: str, infos: list[dict], args) -> None:
    """Write processing report to file."""
    with open(report_path, 'w') as f:
        f.write("# Terminal Domain Extraction Report\n")
        f.write(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"# Input: {args.input}\n")
        f.write(f"# Parameters: similarity={args.similarity}, min_length={args.min_length}\n")
        f.write("# Algorithm: Auto-detect repeat unit length + boundary detection\n")
        f.write(f"# Skipped: {args.skip} sequences\n")
        f.write("#\n")

        # Summary by spidroin type
        type_stats = defaultdict(lambda: {'count': 0, 'total_orig': 0, 'total_new': 0})
        for info in infos:
            sp_type = info['spidroin_type']
            type_stats[sp_type]['count'] += 1
            type_stats[sp_type]['total_orig'] += info['original_length']
            type_stats[sp_type]['total_new'] += info['extracted_length']

        f.write("# Summary by Spidroin Type:\n")
        f.write("# Type\tCount\tAvg_Orig_Len\tAvg_New_Len\tAvg_Ratio\n")
        for sp_type in sorted(type_stats.keys()):
            stats = type_stats[sp_type]
            avg_orig = stats['total_orig'] / stats['count']
            avg_new = stats['total_new'] / stats['count']
            avg_ratio = avg_new / avg_orig * 100 if avg_orig > 0 else 0
            f.write(f"# {sp_type}\t{stats['count']}\t{avg_orig:.1f}\t{avg_new:.1f}\t{avg_ratio:.1f}%\n")
        f.write("#\n")

        # Detailed records
        f.write("# Detailed Processing Log:\n")
        f.write("Original_Header\tNew_Header\tSpidroin_Type\tDomain_Type\tOrig_Len\tNew_Len\tUnit_Len\tRatio\tStatus\n")
        for info in infos:
            ratio = info['extracted_length'] / info['original_length'] * 100 if info['original_length'] > 0 else 0
            unit_len = info.get('unit_length', '-')
            unit_str = str(unit_len) if unit_len else 'fallback'
            f.write(f"{info['original_header']}\t{info['new_header']}\t{info['spidroin_type']}\t"
                    f"{info['domain_type']}\t{info['original_length']}\t{info['extracted_length']}\t"
                    f"{unit_str}\t{ratio:.1f}%\t{info['status']}\n")


def main():
    parser = argparse.ArgumentParser(
        description='Extract non-repetitive CTD/NTD sequences from spidroin sequences.'
    )
    parser.add_argument(
        'input',
        help='Input FASTA file'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output directory (sequences grouped by spidroin type)'
    )
    parser.add_argument(
        '-s', '--skip',
        type=int,
        default=0,
        help='Number of sequences to skip (already processed)'
    )
    parser.add_argument(
        '-l', '--min-length',
        type=int,
        default=50,
        help='Minimum length threshold for extracted sequences (default: 50)'
    )
    parser.add_argument(
        '--similarity',
        type=float,
        default=0.5,
        help='Repeat unit similarity threshold (default: 0.5)'
    )
    parser.add_argument(
        '--no-simplify',
        action='store_true',
        help='Do not simplify headers'
    )

    args = parser.parse_args()

    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Parse input using BioPython
    records = list(SeqIO.parse(args.input, 'fasta'))
    print(f"Loaded {len(records)} sequences from {args.input}", file=sys.stderr)

    # Process sequences and group by spidroin type and domain type
    results_by_type = defaultdict(list)  # key: (spidroin_type, domain_type)
    all_infos = []

    for i, record in enumerate(records):
        if i < args.skip:
            # Keep skipped sequences as-is
            header = record.description if record.description else record.id
            spidroin_type = get_spidroin_type(header)
            domain_type = get_domain_type(header)
            results_by_type[(spidroin_type, domain_type)].append(record)
            all_infos.append({
                'original_header': header,
                'new_header': header,
                'original_length': len(record.seq),
                'extracted_length': len(record.seq),
                'domain_type': domain_type,
                'spidroin_type': spidroin_type,
                'unit_length': None,
                'status': 'skipped',
            })
        else:
            new_record, info = process_record(
                record,
                similarity=args.similarity,
                simplify=not args.no_simplify
            )
            # Filter by minimum length
            if len(new_record.seq) >= args.min_length:
                key = (info['spidroin_type'], info['domain_type'])
                results_by_type[key].append(new_record)
            else:
                info['status'] = 'filtered'
            all_infos.append(info)

    # Write output files by spidroin type and domain type
    all_ctd = []
    all_ntd = []
    for (spidroin_type, domain_type), records_list in results_by_type.items():
        output_file = output_dir / f"{spidroin_type}_{domain_type}.faa"
        SeqIO.write(records_list, output_file, 'fasta')
        print(f"  {spidroin_type}_{domain_type}: {len(records_list)} sequences -> {output_file}", file=sys.stderr)
        # Collect for Pan files
        if domain_type == 'CTD':
            all_ctd.extend(records_list)
        elif domain_type == 'NTD':
            all_ntd.extend(records_list)

    # Write Pan files (all CTD and all NTD combined)
    if all_ctd:
        pan_ctd_file = output_dir / "PanSpidroin_CTD.faa"
        SeqIO.write(all_ctd, pan_ctd_file, 'fasta')
        print(f"  PanSpidroin_CTD: {len(all_ctd)} sequences -> {pan_ctd_file}", file=sys.stderr)
    if all_ntd:
        pan_ntd_file = output_dir / "PanSpidroin_NTD.faa"
        SeqIO.write(all_ntd, pan_ntd_file, 'fasta')
        print(f"  PanSpidroin_NTD: {len(all_ntd)} sequences -> {pan_ntd_file}", file=sys.stderr)

    # Write processing report
    report_path = output_dir / "processing_report.tsv"
    write_report(report_path, all_infos, args)
    print(f"  Report: {report_path}", file=sys.stderr)

    print(f"Done. Output written to {output_dir}/", file=sys.stderr)


if __name__ == '__main__':
    main()
