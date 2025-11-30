#!/usr/bin/env python3
"""
Extract non-repetitive CTD/NTD terminal domain sequences from spidroin sequences.

Algorithm:
- For CTD: scan from C-terminus (end) toward N-terminus, stop at repeat region boundary
- For NTD: scan from N-terminus (start) toward C-terminus, stop at repeat region boundary

Repeat detection:
- Use sliding window (size = sequence_length * ratio) to detect local repetitive regions
- A window is repetitive if any triplet appears >= threshold times within the window
- Default: ratio=0.03, threshold=3 (optimized for ~100aa terminal domains)

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
    -r, --max-repeat N     Triplet repeat threshold (default: 3)
    -w, --window-ratio F   Window size as ratio of seq length (default: 0.03)
    -l, --min-length N     Minimum extracted length (default: 50)
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


def is_repetitive_window(sequence: str, start: int, window_size: int = 30, max_repeat: int = 2) -> bool:
    """
    Check if the window starting at 'start' is repetitive.
    A window is repetitive if any triplet appears >= max_repeat times within it.
    """
    end = min(start + window_size, len(sequence))
    if end - start < 6:  # Need at least 6 aa to have 2 overlapping triplets
        return False

    triplet_counts = Counter()
    for i in range(start, end - 2):
        triplet = sequence[i:i+3]
        triplet_counts[triplet] += 1
        if triplet_counts[triplet] >= max_repeat:
            return True
    return False


def extract_non_repetitive_ctd(sequence: str, max_repeat: int = 3,
                               window_size: int = None, window_ratio: float = 0.03) -> str:
    """
    Extract non-repetitive region from CTD sequence.
    Scan from C-terminus (end) to N-terminus (start).
    Stop when entering a repetitive region.

    Args:
        sequence: Input amino acid sequence
        max_repeat: Triplet must appear this many times in window to be repetitive
        window_size: Fixed size of sliding window (if None, use window_ratio)
        window_ratio: Window size as ratio of sequence length (default: 0.03)
    """
    # Calculate window size
    if window_size is None:
        window_size = max(6, int(len(sequence) * window_ratio))

    if len(sequence) < window_size:
        return sequence

    # Scan from end to start, find where repetitive region ends
    start_pos = 0
    for i in range(len(sequence) - window_size, -1, -1):
        if is_repetitive_window(sequence, i, window_size, max_repeat):
            start_pos = i + window_size
            break

    return sequence[start_pos:]


def extract_non_repetitive_ntd(sequence: str, max_repeat: int = 3,
                               window_size: int = None, window_ratio: float = 0.03) -> str:
    """
    Extract non-repetitive region from NTD sequence.
    Scan from N-terminus (start) to C-terminus (end).
    Stop when entering a repetitive region.

    Args:
        sequence: Input amino acid sequence
        max_repeat: Triplet must appear this many times in window to be repetitive
        window_size: Fixed size of sliding window (if None, use window_ratio)
        window_ratio: Window size as ratio of sequence length (default: 0.03)
    """
    # Calculate window size
    if window_size is None:
        window_size = max(6, int(len(sequence) * window_ratio))

    if len(sequence) < window_size:
        return sequence

    # Scan from start to end, find where repetitive region starts
    end_pos = len(sequence)
    for i in range(len(sequence) - window_size + 1):
        if is_repetitive_window(sequence, i, window_size, max_repeat):
            end_pos = i
            break

    return sequence[:end_pos]


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
    max_repeat: int,
    window_ratio: float,
    simplify: bool
) -> tuple[SeqRecord, dict]:
    """
    Process a single SeqRecord and extract non-repetitive terminal domain.

    Args:
        record: BioPython SeqRecord object
        max_repeat: Maximum allowed triplet repeat count
        window_ratio: Window size as ratio of sequence length
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
        'status': 'unchanged',
    }

    if domain_type == 'CTD':
        extracted = extract_non_repetitive_ctd(seq, max_repeat, window_ratio=window_ratio)
        info['status'] = 'processed'
    elif domain_type == 'NTD':
        extracted = extract_non_repetitive_ntd(seq, max_repeat, window_ratio=window_ratio)
        info['status'] = 'processed'
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
        f.write(f"# Parameters: max_repeat={args.max_repeat}, window_ratio={args.window_ratio}\n")
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
        f.write("Original_Header\tNew_Header\tSpidroin_Type\tDomain_Type\tOrig_Len\tNew_Len\tRatio\tStatus\n")
        for info in infos:
            ratio = info['extracted_length'] / info['original_length'] * 100 if info['original_length'] > 0 else 0
            f.write(f"{info['original_header']}\t{info['new_header']}\t{info['spidroin_type']}\t"
                    f"{info['domain_type']}\t{info['original_length']}\t{info['extracted_length']}\t"
                    f"{ratio:.1f}%\t{info['status']}\n")


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
        '-r', '--max-repeat',
        type=int,
        default=3,
        help='Maximum allowed triplet repeat count in window (default: 3)'
    )
    parser.add_argument(
        '-w', '--window-ratio',
        type=float,
        default=0.03,
        help='Window size as ratio of sequence length (default: 0.03)'
    )
    parser.add_argument(
        '-l', '--min-length',
        type=int,
        default=50,
        help='Minimum length threshold for extracted sequences (default: 50)'
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
                'status': 'skipped',
            })
        else:
            new_record, info = process_record(
                record,
                max_repeat=args.max_repeat,
                window_ratio=args.window_ratio,
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
    for (spidroin_type, domain_type), records_list in results_by_type.items():
        output_file = output_dir / f"{spidroin_type}_{domain_type}.faa"
        SeqIO.write(records_list, output_file, 'fasta')
        print(f"  {spidroin_type}_{domain_type}: {len(records_list)} sequences -> {output_file}", file=sys.stderr)

    # Write processing report
    report_path = output_dir / "processing_report.tsv"
    write_report(report_path, all_infos, args)
    print(f"  Report: {report_path}", file=sys.stderr)

    print(f"Done. Output written to {output_dir}/", file=sys.stderr)


if __name__ == '__main__':
    main()
