"""
junction_validation.py — validate predicted spidroin intron boundaries with
RNA-seq splice junctions extracted from BAM files (N CIGAR ops).

Pipeline per species:
  1. Load typing_results/<species>/<species>.gff (genomic spidroin_gene rows)
     and protein_confirmation/<species>/fullprot_miniprot.gff (intron rows on
     spidroin_full_length sequences = relative coordinates).
  2. For each spidroin_gene, convert the relative intron coordinates to
     genomic coordinates (handling strand).
  3. Open the species BAM and extract introns from N CIGAR ops within the
     spidroin_gene window.
  4. Match predicted introns against extracted junctions with ±tolerance bp,
     count read depth supporting each.
  5. Emit per-spidroin support metrics.

Outputs:
  rna_junction/<species>/junctions.bed    (chrom, start, end, depth, strand)
  rna_junction/<species>/intron_validation.tsv
"""

from collections import defaultdict
from pathlib import Path

from loguru import logger
import polars as pl
import pysam
from tqdm import tqdm
import typer

from spider_silkome_module.config import PROCESSED_DATA_DIR, RAW_DATA_DIR

app = typer.Typer()

DEFAULT_TASK_NAME = "protein_confirmation"
TOLERANCE_BP = 3


# ── helpers ──────────────────────────────────────────────────────────────────


def _parse_attrs(attr_str: str) -> dict[str, str]:
    out = {}
    for part in attr_str.strip().split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def parse_typing_genes(typing_gff: Path) -> dict[str, dict]:
    """Return {spidroin_id: {chrom,start,end,strand}} from the typing GFF."""
    genes: dict[str, dict] = {}
    if not typing_gff.exists():
        return genes
    with open(typing_gff) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 9 or parts[2] != "spidroin_gene":
                continue
            attrs = _parse_attrs(parts[8])
            sid = attrs.get("ID")
            if not sid:
                continue
            genes[sid] = {
                "chrom": parts[0],
                "start": int(parts[3]),
                "end": int(parts[4]),
                "strand": parts[6],
            }
    return genes


def parse_relative_introns(fullprot_gff: Path) -> dict[str, list[tuple[int, int]]]:
    """Return {spidroin_id: [(rel_start, rel_end), ...]} from fullprot_miniprot.gff."""
    out: dict[str, list[tuple[int, int]]] = defaultdict(list)
    if not fullprot_gff.exists():
        return out
    with open(fullprot_gff) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 9 or parts[2] != "intron":
                continue
            sid = parts[0]
            out[sid].append((int(parts[3]), int(parts[4])))
    return out


def relative_to_genome(
    rel_start: int, rel_end: int,
    gene: dict,
) -> tuple[int, int]:
    """Convert spidroin_full_length-relative intron coords to genomic coords."""
    if gene["strand"] == "+":
        return gene["start"] + rel_start - 1, gene["start"] + rel_end - 1
    # '-' strand: full_length sequence is RC'd, so intron rel coords flip.
    return gene["end"] - rel_end + 1, gene["end"] - rel_start + 1


# ── BAM junction extraction ──────────────────────────────────────────────────


def extract_junctions(
    bam: pysam.AlignmentFile,
    chrom: str,
    start: int,
    end: int,
) -> dict[tuple[int, int], int]:
    """
    Walk reads overlapping [start, end] (1-based inclusive) and tally junctions
    derived from N CIGAR operations.  Returns {(intron_start, intron_end): depth}
    in 1-based inclusive coords.
    """
    counts: dict[tuple[int, int], int] = defaultdict(int)
    try:
        iterator = bam.fetch(chrom, max(0, start - 1), end)
    except (ValueError, KeyError):
        return counts
    for read in iterator:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.cigartuples is None:
            continue
        ref_pos = read.reference_start  # 0-based
        for op, length in read.cigartuples:
            if op in (0, 7, 8, 2):  # M, =, X, D consume reference
                ref_pos += length
            elif op == 3:  # N
                istart = ref_pos + 1                # 1-based intron start
                iend = ref_pos + length             # 1-based intron end
                if iend >= start and istart <= end:
                    counts[(istart, iend)] += 1
                ref_pos += length
            # other ops (I, S, H, P) consume neither or only query
    return counts


def match_predicted_to_observed(
    predicted: list[tuple[int, int]],
    observed: dict[tuple[int, int], int],
    tolerance: int = TOLERANCE_BP,
) -> list[tuple[int, int, int]]:
    """
    For each predicted intron, find the best-matching observed junction
    within tolerance bp and return [(predicted_start, predicted_end, depth)]
    where depth=0 indicates no supporting junction.
    """
    results: list[tuple[int, int, int]] = []
    for ps, pe in predicted:
        best_depth = 0
        for (os_, oe), d in observed.items():
            if abs(os_ - ps) <= tolerance and abs(oe - pe) <= tolerance:
                if d > best_depth:
                    best_depth = d
        results.append((ps, pe, best_depth))
    return results


# ── BAM discovery ────────────────────────────────────────────────────────────


def find_species_bam(species_id: str, bam_root: Path) -> Path | None:
    """Look up sorted BAM under bam_root/<species_id>/outs/sorted.bam."""
    for layout in (
        bam_root / species_id / "outs" / "sorted.bam",
        bam_root / species_id / "sorted.bam",
    ):
        if layout.exists() and Path(str(layout) + ".bai").exists():
            return layout
        if layout.exists() and Path(str(layout) + ".csi").exists():
            return layout
    return None


# ── per-species orchestration ────────────────────────────────────────────────


def process_species(
    species_name: str,
    confirmation_dir: Path,
    typing_dir: Path,
    bam_root: Path,
    output_species_dir: Path,
    tolerance: int = TOLERANCE_BP,
) -> int:
    """Validate introns for one species; return number of spidroins processed."""
    bam_path = find_species_bam(species_name, bam_root)
    if bam_path is None:
        logger.warning(f"[{species_name}] no indexed BAM, skipping")
        return 0

    typing_gff = typing_dir / species_name / f"{species_name}.gff"
    fullprot_gff = confirmation_dir / species_name / "fullprot_miniprot.gff"
    if not typing_gff.exists() or not fullprot_gff.exists():
        logger.warning(f"[{species_name}] missing GFFs, skipping")
        return 0

    genes = parse_typing_genes(typing_gff)
    rel_introns = parse_relative_introns(fullprot_gff)

    output_species_dir.mkdir(parents=True, exist_ok=True)
    bed_lines: list[str] = []
    rows: list[dict] = []

    bam = pysam.AlignmentFile(str(bam_path), "rb")
    try:
        for sid, gene in genes.items():
            preds_rel = rel_introns.get(sid, [])
            if not preds_rel:
                rows.append({
                    "spidroin_id": sid,
                    "n_introns_predicted": 0,
                    "n_introns_rnaseq_supported": 0,
                    "junction_support_rate": None,
                    "min_junction_depth": None,
                    "median_junction_depth": None,
                })
                continue
            preds_genomic = [relative_to_genome(s, e, gene) for s, e in preds_rel]
            observed = extract_junctions(bam, gene["chrom"], gene["start"], gene["end"])
            matched = match_predicted_to_observed(preds_genomic, observed, tolerance)
            depths = [d for _, _, d in matched]
            n_supported = sum(1 for d in depths if d > 0)
            for (gs, ge), d in zip(preds_genomic, depths):
                bed_lines.append(
                    f"{gene['chrom']}\t{gs - 1}\t{ge}\t{sid}\t{d}\t{gene['strand']}"
                )
            rows.append({
                "spidroin_id": sid,
                "n_introns_predicted": len(preds_rel),
                "n_introns_rnaseq_supported": n_supported,
                "junction_support_rate": round(n_supported / len(preds_rel), 4),
                "min_junction_depth": min(depths) if depths else 0,
                "median_junction_depth": int(sorted(depths)[len(depths) // 2]) if depths else 0,
            })
    finally:
        bam.close()

    (output_species_dir / "junctions.bed").write_text("\n".join(bed_lines) + "\n")
    pl.DataFrame(rows).write_csv(
        output_species_dir / "intron_validation.tsv", separator="\t"
    )
    return len(rows)


# ── CLI ──────────────────────────────────────────────────────────────────────


@app.command()
def main(
    confirmation_dir: Path = PROCESSED_DATA_DIR / DEFAULT_TASK_NAME / "confirmation",
    typing_dir: Path = PROCESSED_DATA_DIR / "typing_results",
    bam_root: Path = RAW_DATA_DIR / "Spider_anno",
    output_dir: Path = PROCESSED_DATA_DIR / DEFAULT_TASK_NAME / "rna_junction",
    species: str | None = None,
    tolerance: int = TOLERANCE_BP,
):
    """Validate predicted intron boundaries with RNA-seq splice junctions."""
    species_dirs = sorted(d for d in confirmation_dir.iterdir() if d.is_dir())
    if species is not None:
        species_dirs = [d for d in species_dirs if d.name == species]
    if not species_dirs:
        logger.error(f"No species under {confirmation_dir} matching filter")
        raise typer.Exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)
    total = 0
    for sp_dir in tqdm(species_dirs, desc="Validating junctions"):
        n = process_species(
            sp_dir.name, confirmation_dir, typing_dir,
            bam_root, output_dir / sp_dir.name, tolerance,
        )
        logger.info(f"[{sp_dir.name}] {n} spidroins junction-validated")
        total += n
    logger.success(f"Junction validation complete for {total} spidroins")


if __name__ == "__main__":
    app()
