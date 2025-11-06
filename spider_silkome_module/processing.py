"""
Data Processing Function Module
"""

from collections import defaultdict
import re
from typing import List

from natsort import natsort_keygen, ns

from spider_silkome_module.models import GFFData, Position


def extract_positions_from_gff(
    gff_data: List[GFFData], positive_threshold: float = 0.85
) -> List[Position]:
    """
    Extracting the positional information of spidroin from GFF data

    According to the alignment results of the C-terminus and N-terminus, extract the start and end positions of each chromosome-strand combination, as well as the type of spidroin.

    Core logic explanation:
    - Sense strand (+): Gene runs from 5' to 3', N-terminus at the front (start), C-terminus at the back (end).
    - Anti-sense strand (-): Gene runs from 3' to 5', N-terminus at the back (end), C-terminus at the front (start).

    Position recording rules:
    | Type | Strand | Record Position | Reason                           |
    |------|--------|----------|--------------------------------|
    | CTD  | +      | end      | C-terminal domain at the end of the gene, extends backward        |
    | CTD  | -      | start    | C-terminal domain at the start of the gene (smaller coordinate in gene group), extends forward |
    | NTD  | +      | start    | N-terminal domain at the start of the gene, extends forward        |
    | NTD  | -      | end      | N-terminal domain at the end of the gene (larger coordinate in gene group), extends backward |

    Parameters
    ----------
    gff_data : List[GFFData]
        miniprot output GFF data list
    positive_threshold : float, optional
        Quality threshold, filter out alignment results below this threshold, default is 0.85

    Returns
    -------
    List[Position]
        List of positions sorted by chromosome and strand

    Examples
    --------
    >>> positions = extract_positions_from_gff(spidroin_gff_data, positive_threshold=0.85)
    >>> print(f"Found {len(positions)} chromosome-strand combinations")
    """
    # Store position information for each chromosome+strand
    positions_dict = defaultdict(lambda: {"start": defaultdict(int), "end": defaultdict(int)})

    for aln in gff_data:
        # Filter out low quality alignments
        if aln.attributes.Positive < positive_threshold:
            continue

        chr_id = aln.seqid
        strand = aln.strand
        key = (chr_id, strand)

        # Determine if it is C-terminal or N-terminal
        domain_type = aln.attributes.Target[-1].split(" ")[0]

        if domain_type == "CTD":  # C-terminal domain
            if strand == "+":  # Forward: C-terminal at the end, record end position
                pos_value = aln.end
                positions_dict[key]["end"][pos_value] += 1
            else:  # Reverse: C-terminal at the start (smaller coordinate in gene group), record start position
                pos_value = aln.start
                positions_dict[key]["start"][pos_value] += 1

        elif domain_type == "NTD":  # N-terminal domain
            if strand == "+":  # Forward: N-terminal at the start, record start position
                pos_value = aln.start
                positions_dict[key]["start"][pos_value] += 1
            else:  # Reverse: N-terminal at the end (larger coordinate in gene group), record end position
                pos_value = aln.end
                positions_dict[key]["end"][pos_value] += 1

    # Convert to Position object list
    positions = []
    for (chr_id, strand), pos_data in positions_dict.items():
        positions.append(
            Position(
                chr=chr_id, strand=strand, start=dict(pos_data["start"]), end=dict(pos_data["end"])
            )
        )

    # Sort by chromosome and strand with custom priority (chr > ctg > others)
    nkey = natsort_keygen(alg=ns.IGNORECASE)

    def chr_sort_key(x):
        name = x.chr
        if re.match(r"^chr", name, flags=re.IGNORECASE):
            priority = 0
        elif re.match(r"^ctg", name, flags=re.IGNORECASE):
            priority = 1
        else:
            priority = 2
        return (priority, nkey(name), x.strand)

    positions.sort(key=chr_sort_key)

    return positions
