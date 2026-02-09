#!/usr/bin/env python3
"""
Merge original spidroin GFF with Augustus gene predictions.
"""
import re
from pathlib import Path


def parse_augustus_gff(augustus_gff_path: Path) -> dict:
    """
    Parse Augustus GFF output and extract gene predictions.
    Returns a dict mapping sequence_id to list of features.
    """
    predictions = {}
    if not augustus_gff_path.exists():
        return predictions

    with open(augustus_gff_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            seq_id, source, feature_type, start, end, score, strand, phase, attrs = parts
            if seq_id not in predictions:
                predictions[seq_id] = []
            predictions[seq_id].append({
                'type': feature_type,
                'start': int(start),
                'end': int(end),
                'score': score,
                'strand': strand,
                'phase': phase,
                'attrs': attrs
            })
    return predictions


def local_to_genomic(local_start: int, local_end: int, gene_start: int, gene_end: int, strand: str) -> tuple:
    """
    Convert local coordinates (relative to extracted sequence) to genomic coordinates.

    For + strand: genomic_pos = gene_start + local_pos - 1
    For - strand: genomic_pos = gene_end - local_pos + 1 (reversed)
    """
    if strand == '+':
        genomic_start = gene_start + local_start - 1
        genomic_end = gene_start + local_end - 1
    else:
        # For minus strand, the extracted sequence is reverse complemented
        genomic_start = gene_end - local_end + 1
        genomic_end = gene_end - local_start + 1
    return genomic_start, genomic_end


def merge_gff_files(original_gff_path: Path, augustus_preds: dict, output_gff_path: Path) -> None:
    """
    Merge original spidroin GFF with Augustus predictions.

    Args:
        original_gff_path: Path to original spidroin analysis GFF
        augustus_preds: Dict of Augustus predictions by sequence ID
        output_gff_path: Path for merged output GFF
    """
    original_genes = {}

    if original_gff_path.exists():
        with open(original_gff_path) as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                attrs = parts[8]
                id_match = re.search(r'ID=([^;]+)', attrs)
                parent_match = re.search(r'Parent=([^;]+)', attrs)

                if parts[2] == 'gene' and id_match:
                    gene_id = id_match.group(1)
                    original_genes[gene_id] = {
                        'gene_line': line.strip(),
                        'features': [],
                        'scaffold': parts[0],
                        'strand': parts[6],
                        'gene_start': int(parts[3]),
                        'gene_end': int(parts[4])
                    }
                elif parent_match:
                    parent_id = parent_match.group(1)
                    if parent_id in original_genes:
                        original_genes[parent_id]['features'].append(line.strip())

    with open(output_gff_path, 'w') as f:
        f.write("##gff-version 3\n")

        for gene_id, gene_data in original_genes.items():
            scaffold = gene_data['scaffold']
            strand = gene_data['strand']
            gene_start = gene_data['gene_start']
            gene_end = gene_data['gene_end']

            f.write(gene_data['gene_line'] + '\n')

            # Write original domain features (NTD, CTD)
            for feat in gene_data['features']:
                f.write(feat + '\n')

            # Add Augustus predictions if available
            if gene_id in augustus_preds:
                aug_feats = augustus_preds[gene_id]
                mrna_id = f"{gene_id}_mRNA"

                cds_features = [ft for ft in aug_feats if ft['type'] == 'CDS']
                if cds_features:
                    # Convert all CDS coordinates to genomic
                    genomic_cds = []
                    for cds in cds_features:
                        g_start, g_end = local_to_genomic(
                            cds['start'], cds['end'], gene_start, gene_end, strand
                        )
                        genomic_cds.append({
                            'start': min(g_start, g_end),
                            'end': max(g_start, g_end),
                            'score': cds['score'],
                            'phase': cds['phase']
                        })

                    # Sort CDS by start position
                    genomic_cds.sort(key=lambda x: x['start'])

                    mrna_start = min(c['start'] for c in genomic_cds)
                    mrna_end = max(c['end'] for c in genomic_cds)

                    f.write(f"{scaffold}\tAugustus\tmRNA\t{mrna_start}\t{mrna_end}\t.\t{strand}\t.\tID={mrna_id};Parent={gene_id}\n")

                    for i, cds in enumerate(genomic_cds, 1):
                        cds_id = f"{gene_id}_CDS_{i}"
                        f.write(f"{scaffold}\tAugustus\tCDS\t{cds['start']}\t{cds['end']}\t{cds['score']}\t{strand}\t{cds['phase']}\tID={cds_id};Parent={mrna_id}\n")

    print(f"Wrote merged GFF: {output_gff_path}")
