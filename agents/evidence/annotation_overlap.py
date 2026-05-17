"""
annotation_overlap.py — check overlap between spidroin candidates and genome annotations.

Two annotation versions are supported:
  - New (spider_anno2): <anno_new_dir>/<NNN>spd/gene.gff  (plain GFF3, source=ReName)
  - Old (01.ref_gff):   <anno_old_dir>/<NNN>.<species>/*.gene.gff3.gz  (gzipped GFF3, source=Braker3)

Both annotation sources are loaded and merged into a single gene index per species.
This ensures that 5′- and 3′-boundary predictions from either annotation version
are available as fallback positions for start/stop codon verification.

Species IDs follow the project convention: "{NNN}.{species}" (e.g., "001.Allagelena_difficilis").
The NNN prefix maps to spider_anno2 subdirectory "<NNN>spd".
"""

import gzip
from pathlib import Path

from loguru import logger


def _parse_gff_genes(path: Path, compressed: bool = False) -> list[dict]:
    """
    Read gene-level features from a GFF3 file.
    Returns list of dicts: {seqid, start, end, gene_id}.
    Input coordinates are 1-based closed; returned as-is.
    """
    open_fn = gzip.open if compressed else open
    genes = []
    try:
        with open_fn(path, "rt") as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split("\t")
                if len(parts) < 9:
                    continue
                if parts[2] != "gene":
                    continue
                seqid = parts[0]
                try:
                    start = int(parts[3])
                    end = int(parts[4])
                except ValueError:
                    continue
                # Extract gene ID from attributes
                gene_id = ""
                for attr in parts[8].split(";"):
                    attr = attr.strip()
                    if attr.startswith("ID="):
                        gene_id = attr[3:].strip()
                        break
                genes.append({"seqid": seqid, "start": start, "end": end, "gene_id": gene_id})
    except Exception as exc:
        logger.warning(f"Failed to parse annotation file {path}: {exc}")
    return genes


def _find_anno_paths(
    species_id: str,
    anno_new_dir: Path | None,
    anno_old_dir: Path | None,
) -> list[tuple[Path, bool]]:
    """
    Locate all available gene annotation files for a species.
    Returns a list of (path, is_gzipped) for each source that exists.

    Both new (spider_anno2) and old (01.ref_gff) annotations are returned when
    available so that gene boundaries from both sources contribute to the index.
    """
    nnn = species_id.split(".")[0]
    results: list[tuple[Path, bool]] = []

    if anno_new_dir is not None:
        new_path = anno_new_dir / f"{nnn}spd" / "gene.gff"
        if new_path.exists():
            results.append((new_path, False))

    if anno_old_dir is not None:
        old_species_dir = anno_old_dir / species_id
        if old_species_dir.is_dir():
            gz_files = list(old_species_dir.glob("*.gene.gff3.gz"))
            if gz_files:
                results.append((gz_files[0], True))

    return results


def build_gene_index(
    species_id: str,
    anno_new_dir: Path | None,
    anno_old_dir: Path | None,
) -> dict[str, list[dict]]:
    """
    Build a per-seqid gene index for a species by merging genes from all
    available annotation sources (both new and old).

    Returns {seqid: [{"start": int, "end": int, "gene_id": str}, ...]}.
    Returns empty dict if no annotation file is found.
    """
    paths = _find_anno_paths(species_id, anno_new_dir, anno_old_dir)
    if not paths:
        logger.debug(f"No annotation found for species {species_id}")
        return {}

    all_genes: list[dict] = []
    for path, compressed in paths:
        genes = _parse_gff_genes(path, compressed=compressed)
        logger.debug(f"Loaded {len(genes)} gene records for {species_id} from {path}")
        all_genes.extend(genes)

    index: dict[str, list[dict]] = {}
    for g in all_genes:
        index.setdefault(g["seqid"], []).append(
            {"start": g["start"], "end": g["end"], "gene_id": g["gene_id"]}
        )
    return index


def query_overlap(
    seqid: str,
    locus_start: int,
    locus_end: int,
    gene_index: dict[str, list[dict]],
    min_overlap_frac: float = 0.0,
) -> tuple[str, float, int | None, int | None]:
    """
    Find the gene with the greatest overlap with the locus region.
    Coordinates are 1-based closed intervals.

    Returns (gene_id, overlap_frac, gene_start, gene_end) of the best-overlapping gene.
    Returns ("", 0.0, None, None) if no overlapping gene is found.
    min_overlap_frac: minimum overlap fraction required to count (default 0 = any overlap).
    """
    genes = gene_index.get(seqid, [])
    locus_len = locus_end - locus_start + 1
    best_gene_id = ""
    best_frac = 0.0
    best_start: int | None = None
    best_end: int | None = None

    for gene in genes:
        overlap_start = max(locus_start, gene["start"])
        overlap_end = min(locus_end, gene["end"])
        overlap_len = max(0, overlap_end - overlap_start + 1)
        if overlap_len == 0:
            continue
        # Overlap fraction relative to the locus length
        frac = overlap_len / locus_len
        if frac > best_frac:
            best_frac = frac
            best_gene_id = gene["gene_id"]
            best_start = gene["start"]
            best_end = gene["end"]

    if best_frac < min_overlap_frac:
        return "", 0.0, None, None
    return best_gene_id, round(best_frac, 4), best_start, best_end


def annotate_loci(
    loci: list[dict],
    anno_new_dir: Path | None,
    anno_old_dir: Path | None,
) -> list[dict]:
    """
    Add annotation overlap columns to each locus dict.
    Modifies loci in-place and returns the list.
    Added keys: anno_gene_id (str), anno_overlap_frac (float),
                anno_start (int | None), anno_end (int | None).
    """
    # Build gene indices per species (cache to avoid re-parsing)
    gene_indices: dict[str, dict] = {}

    for locus in loci:
        species_id = locus["species"]
        if species_id not in gene_indices:
            gene_indices[species_id] = build_gene_index(species_id, anno_new_dir, anno_old_dir)

        gene_id, overlap_frac, anno_start, anno_end = query_overlap(
            seqid=locus["seqid"],
            locus_start=locus["start"],
            locus_end=locus["end"],
            gene_index=gene_indices[species_id],
        )
        locus["anno_gene_id"] = gene_id
        locus["anno_overlap_frac"] = overlap_frac
        locus["anno_start"] = anno_start
        locus["anno_end"] = anno_end

    return loci


def upgrade_loci_from_annotation(
    loci: list[dict],
    max_ntd_offset_bp: int = 5000,
    min_anno_extension_bp: int = 10000,
) -> list[dict]:
    """
    Upgrade N-terminal-only loci to Full_length when a gene annotation
    covers the locus from its 5' end and extends substantially to the 3' end.

    Conditions (all must be met):
      1. Locus is N-terminal only: has_ntd=True, has_ctd=False.
      2. A gene annotation is found: anno_gene_id is set (populated by annotate_loci).
      3. The NTD is at the 5' end of the annotation:
           '+' strand → ntd_start is within max_ntd_offset_bp of anno_start
           '-' strand → ntd_end   is within max_ntd_offset_bp of anno_end
      4. The annotation extends at least min_anno_extension_bp beyond the NTD:
           '+' strand → anno_end  - ntd_end   ≥ min_anno_extension_bp
           '-' strand → ntd_start - anno_start ≥ min_anno_extension_bp

    When conditions are met:
      - Sets ctd_start and ctd_end to a 3-nt proxy at the annotation's 3' boundary
        so the downstream codon-check step can verify the stop codon.  The proxy
        is positioned to cover the stop codon whether the gene caller includes it
        inside the gene model or just outside it:
            '+' strand: ctd_start = ctd_end = anno_end  - 3  (1-based)
            '-' strand: ctd_start = ctd_end = anno_start + 3  (1-based)
      - Extends locus end ('+') or start ('-') to the annotation boundary.
      - Sets completeness = "Full_length", needs_review = True.
      - Sets ctd_profile = "annotation_proxy" to distinguish from nHMMER evidence.
      - Appends note "completeness_from_anno".

    Requires annotate_loci() to have run first (anno_start / anno_end must be set).
    """
    for locus in loci:
        if not (locus.get("has_ntd") and not locus.get("has_ctd")):
            continue
        if not locus.get("anno_gene_id"):
            continue

        anno_start = locus.get("anno_start")
        anno_end = locus.get("anno_end")
        if anno_start is None or anno_end is None:
            continue

        ntd_start = locus.get("ntd_start")
        ntd_end = locus.get("ntd_end")
        if ntd_start is None or ntd_end is None:
            continue

        strand = locus.get("strand", "+")

        if strand == "+":
            ntd_near_start = abs(ntd_start - anno_start) <= max_ntd_offset_bp
            anno_extends = (anno_end - ntd_end) >= min_anno_extension_bp
            if not (ntd_near_start and anno_extends):
                continue
            # Proxy CTD at annotation 3' end (covers stop codon for both
            # gene-model conventions: stop codon included vs excluded).
            proxy = anno_end - 3
            locus["ctd_start"] = proxy
            locus["ctd_end"] = proxy
            locus["end"] = anno_end
        else:
            ntd_near_end = abs(ntd_end - anno_end) <= max_ntd_offset_bp
            anno_extends = (ntd_start - anno_start) >= min_anno_extension_bp
            if not (ntd_near_end and anno_extends):
                continue
            proxy = anno_start + 3
            locus["ctd_start"] = proxy
            locus["ctd_end"] = proxy
            locus["start"] = anno_start

        locus["ctd_profile"] = "annotation_proxy"
        locus["ctd_type_prefix"] = None
        locus["ctd_e_value"] = None
        locus["completeness"] = "Full_length"
        locus["needs_review"] = True
        existing = locus.get("notes") or ""
        sep = "; " if existing else ""
        locus["notes"] = existing + sep + "completeness_from_anno"

    return loci
