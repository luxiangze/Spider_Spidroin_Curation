#!/usr/bin/env python3
"""
Spidroin Repeat Motif Visualization Tool

This script identifies and visualizes repetitive motifs in spidroin protein sequences.
Based on common spider silk protein repeat patterns.

Usage:
    python spidroin_repeat_visual.py -i input.fasta -o output_dir [options]
"""

import argparse
import re
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

from Bio import SeqIO
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse, Rectangle


# ============================================================================
# Motif Pattern Definitions (ordered by priority - longer/more specific first)
# ============================================================================

@dataclass
class MotifPattern:
    """Definition of a repeat motif pattern"""
    name: str
    regex: str
    color: str
    priority: int  # Lower number = higher priority


# Define motifs based on the reference image, ordered by priority
SPIDROIN_MOTIFS = OrderedDict([
    # High priority: longer/more specific patterns first
    ("GAAA", MotifPattern("GAAA", r"GA{3,}", "#9ACD32", 1)),  # Yellow-green
    ("(GA)n", MotifPattern("(GA)n", r"(?:GA){2,}", "#87CEEB", 2)),  # Light blue
    ("(GS)n", MotifPattern("(GS)n", r"(?:GS){2,}", "#DDA0DD", 3)),  # Light purple
    ("PolyA", MotifPattern("PolyA/S", r"A{4,}", "#FF69B4", 4)),  # Pink - polyA
    ("PolyS", MotifPattern("PolyA/S", r"S{4,}", "#FF69B4", 5)),  # Pink - polyS
    ("PolyT", MotifPattern("polyT", r"T{4,}", "#4169E1", 6)),  # Blue
    ("AnSn", MotifPattern("AnSn", r"(?:A{1,3}S{1,3}){1,}", "#8B0000", 7)),  # Dark red
    ("QnX", MotifPattern("QnX", r"Q{2,}[A-Z]", "#FF0000", 8)),  # Red
    ("SnX", MotifPattern("SnX", r"S{2,}[A-Z]", "#FF4500", 9)),  # Red-orange
    ("TnX", MotifPattern("TnX", r"T{2,}[A-Z]", "#FFA500", 10)),  # Orange
    ("AnX", MotifPattern("AnX", r"A{2,}[A-Z]", "#FFD700", 11)),  # Yellow
    ("GGX", MotifPattern("GGX", r"GG[A-Z]", "#00008B", 12)),  # Dark blue
    ("GA", MotifPattern("GA", r"GA", "#228B22", 13)),  # Green
    ("GX", MotifPattern("GX", r"G[A-Z]", "#800080", 14)),  # Purple
])

# Color for unmatched regions (spacer)
SPACER_COLOR = "#FFFFFF"  # White
TERMINAL_COLOR = "#17202A"  # Dark for N/C terminal


# ============================================================================
# Sequence Parsing
# ============================================================================

def read_sequences(fasta_path: Path) -> List[Tuple[str, str]]:
    """
    Read sequences from a FASTA file.

    Returns:
        List of (sequence_id, sequence) tuples
    """
    sequences = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        sequences.append((record.id, str(record.seq)))
    return sequences


# ============================================================================
# Motif Detection
# ============================================================================

@dataclass
class MotifMatch:
    """A matched motif region in the sequence"""
    start: int
    end: int
    motif_name: str
    sequence: str
    color: str


def find_all_motifs(sequence: str, motifs: OrderedDict = SPIDROIN_MOTIFS) -> List[MotifMatch]:
    """
    Find all motif matches in a sequence with priority-based resolution.

    Higher priority motifs (lower priority number) take precedence when overlapping.
    """
    # Track which positions are already claimed
    claimed = [False] * len(sequence)
    matches = []

    # Process motifs in priority order
    for motif_key, motif in motifs.items():
        pattern = re.compile(motif.regex)
        for match in pattern.finditer(sequence):
            start, end = match.start(), match.end()

            # Check if this region is already claimed
            if any(claimed[start:end]):
                continue

            # Claim this region
            for i in range(start, end):
                claimed[i] = True

            matches.append(MotifMatch(
                start=start,
                end=end,
                motif_name=motif.name,
                sequence=match.group(),
                color=motif.color
            ))

    # Sort matches by position
    matches.sort(key=lambda x: x.start)

    return matches


def segment_sequence(sequence: str, matches: List[MotifMatch]) -> List[MotifMatch]:
    """
    Segment the entire sequence including unmatched regions (spacers).
    """
    segments = []
    current_pos = 0

    for match in matches:
        # Add spacer region if there's a gap
        if match.start > current_pos:
            spacer_seq = sequence[current_pos:match.start]
            segments.append(MotifMatch(
                start=current_pos,
                end=match.start,
                motif_name="spacer",
                sequence=spacer_seq,
                color=SPACER_COLOR
            ))

        segments.append(match)
        current_pos = match.end

    # Add final spacer if needed
    if current_pos < len(sequence):
        spacer_seq = sequence[current_pos:]
        segments.append(MotifMatch(
            start=current_pos,
            end=len(sequence),
            motif_name="spacer",
            sequence=spacer_seq,
            color=SPACER_COLOR
        ))

    return segments


def analyze_sequence(sequence: str) -> Tuple[List[MotifMatch], Dict[str, int]]:
    """
    Analyze a sequence and return segments with motif statistics.
    """
    matches = find_all_motifs(sequence)
    segments = segment_sequence(sequence, matches)

    # Count motif occurrences
    motif_counts = {}
    for seg in segments:
        if seg.motif_name != "spacer":
            motif_counts[seg.motif_name] = motif_counts.get(seg.motif_name, 0) + 1

    return segments, motif_counts


# ============================================================================
# Visualization
# ============================================================================

def plot_single_sequence(
    seq_id: str,
    segments: List[MotifMatch],
    output_path: Path,
    figwidth: int = 20,
    figheight: int = 2
):
    """
    Create a linear sequence visualization with colored ellipses for motifs.
    Style similar to published spidroin repeat diagrams.
    """
    plt.rcParams["pdf.fonttype"] = 42

    # Calculate total sequence length
    total_len = sum(len(seg.sequence) for seg in segments)
    if total_len == 0:
        return

    fig, ax = plt.subplots(figsize=(figwidth, figheight))

    # Background bar (dark brown color like in the reference image)
    bar_height = 0.5
    bar_y = 0.5
    bg_color = "#8B4513"  # Saddle brown
    ax.add_patch(Rectangle(
        (0, bar_y - bar_height / 2), total_len, bar_height,
        facecolor=bg_color, edgecolor="black", linewidth=1.5
    ))

    # Draw motifs as ellipses
    ellipse_height = bar_height * 0.9

    for seg in segments:
        if seg.motif_name == "spacer":
            continue

        # Position at center of motif
        x_center = (seg.start + seg.end) / 2
        width = seg.end - seg.start

        # Minimum width for visibility, make ellipses rounder
        min_width = total_len * 0.004
        ellipse_width = max(width * 1.2, min_width)

        ellipse = Ellipse(
            (x_center, bar_y),
            width=ellipse_width,
            height=ellipse_height,
            facecolor=seg.color,
            edgecolor="black",
            linewidth=0.3
        )
        ax.add_patch(ellipse)

    # Configure axes
    ax.set_xlim(-total_len * 0.15, total_len * 1.02)
    ax.set_ylim(0, 1)
    ax.set_aspect("auto")
    ax.axis("off")

    # Add sequence ID as label on the left
    ax.text(
        -total_len * 0.01, bar_y, seq_id,
        ha="right", va="center", fontsize=11, fontweight="bold"
    )

    # Add scale bar
    scale_len = round(total_len / 10, -2)  # Round to nearest 100
    if scale_len < 100:
        scale_len = 100
    scale_x = total_len * 0.02
    scale_y = 0.12
    ax.plot([scale_x, scale_x + scale_len], [scale_y, scale_y], "k-", linewidth=2)
    ax.text(
        scale_x + scale_len / 2, scale_y - 0.06,
        f"{int(scale_len)} aa", ha="center", va="top", fontsize=9
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_legend(motifs: OrderedDict, output_path: Path, motif_counts: Dict[str, int] = None):
    """
    Create a legend showing all motif types with their colors.
    """
    plt.rcParams["pdf.fonttype"] = 42

    # Get unique motif names and colors
    unique_motifs = []
    seen_names = set()
    for motif in motifs.values():
        if motif.name not in seen_names:
            count = motif_counts.get(motif.name, 0) if motif_counts else 0
            unique_motifs.append((motif.name, motif.color, count))
            seen_names.add(motif.name)

    # Add spacer
    spacer_count = motif_counts.get("spacer", 0) if motif_counts else 0
    unique_motifs.append(("spacer", SPACER_COLOR, spacer_count))

    fig, ax = plt.subplots(figsize=(4, len(unique_motifs) * 0.5))

    for i, (name, color, count) in enumerate(unique_motifs):
        if count > 0:
            label = f"{name} ({count})"
        else:
            label = name
        ax.barh(i, 1, color=color, edgecolor="black", linewidth=0.5)
        ax.text(1.1, i, label, va="center", fontsize=10)

    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_xlim(0, 3)
    ax.set_ylim(-0.5, len(unique_motifs) - 0.5)

    # Remove spines
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.set_title("Motif Legend", fontsize=12)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


def output_motif_table(
    seq_id: str,
    segments: List[MotifMatch],
    output_path: Path
):
    """
    Output a TSV file with motif positions and sequences.
    """
    with open(output_path, "w") as f:
        f.write("sequence_id\tstart\tend\tmotif_name\tlength\tsequence\n")
        for seg in segments:
            f.write(f"{seq_id}\t{seg.start}\t{seg.end}\t{seg.motif_name}\t{len(seg.sequence)}\t{seg.sequence}\n")


def output_summary(
    seq_id: str,
    sequence: str,
    motif_counts: Dict[str, int],
    output_path: Path
):
    """
    Output a summary of motif statistics.
    """
    with open(output_path, "w") as f:
        f.write("# Spidroin Repeat Analysis Summary\n")
        f.write(f"# Sequence: {seq_id}\n")
        f.write(f"# Length: {len(sequence)} residues\n\n")
        f.write("motif_name\tcount\n")
        for name, count in sorted(motif_counts.items(), key=lambda x: -x[1]):
            f.write(f"{name}\t{count}\n")


# ============================================================================
# Main Entry Point
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Spidroin repeat motif visualization tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Analyze a single FASTA file
    python spidroin_repeat_visual.py -i spidroin.fasta -o output/

    # Process multiple sequences with custom figure size
    python spidroin_repeat_visual.py -i sequences.fasta -o output/ --width 50 --height 8
        """
    )

    parser.add_argument(
        "-i", "--input",
        type=Path,
        required=True,
        help="Input FASTA file containing spidroin protein sequence(s)"
    )

    parser.add_argument(
        "-o", "--output",
        type=Path,
        required=True,
        help="Output directory for results"
    )

    parser.add_argument(
        "--width",
        type=int,
        default=35,
        help="Figure width in inches (default: 35)"
    )

    parser.add_argument(
        "--height",
        type=int,
        default=6,
        help="Figure height in inches (default: 6)"
    )

    parser.add_argument(
        "--format",
        choices=["pdf", "png", "svg"],
        default="pdf",
        help="Output figure format (default: pdf)"
    )

    args = parser.parse_args()

    # Create output directory
    args.output.mkdir(parents=True, exist_ok=True)

    # Read sequences
    print(f"Reading sequences from {args.input}...")
    sequences = read_sequences(args.input)
    print(f"Found {len(sequences)} sequence(s)")

    # Process each sequence
    all_motif_counts = {}

    for seq_id, sequence in sequences:
        print(f"\nProcessing: {seq_id}")
        print(f"  Length: {len(sequence)} residues")

        # Analyze sequence
        segments, motif_counts = analyze_sequence(sequence)

        # Merge counts
        for name, count in motif_counts.items():
            all_motif_counts[name] = all_motif_counts.get(name, 0) + count

        # Clean sequence ID for filenames
        safe_id = re.sub(r'[^\w\-]', '_', seq_id)

        # Generate outputs
        figure_path = args.output / f"{safe_id}_repeat.{args.format}"
        table_path = args.output / f"{safe_id}_motifs.tsv"
        summary_path = args.output / f"{safe_id}_summary.txt"

        plot_single_sequence(seq_id, segments, figure_path, args.width, args.height)
        output_motif_table(seq_id, segments, table_path)
        output_summary(seq_id, sequence, motif_counts, summary_path)

        print(f"  Motifs found: {sum(motif_counts.values())}")
        print(f"  Output: {figure_path}")

    # Generate combined legend
    legend_path = args.output / f"motif_legend.{args.format}"
    plot_legend(SPIDROIN_MOTIFS, legend_path, all_motif_counts)
    print(f"\nLegend saved to: {legend_path}")

    print("\nDone!")


if __name__ == "__main__":
    main()
