"""
rna_support.py — query RNA-seq BigWig signal over spidroin candidate loci.

Data sources (MinIO via mc):
  BGI Illumina (102 species):
    rustfs/spider-genome/BGI_RNA_align/Spider_anno/{NNN}.{species}/outs/bam.bw
  ONT mRNA (10 species):
    rustfs/spider-genome/ONT_RNA_align/bigwig/{species}.bw

BigWig files are downloaded on demand to a local cache directory (default: cache/bw/).
Downloaded files are skipped if already present (cache-first).

Species ID format in this project: "{NNN}.{species}" (e.g., "064.Araneus_ventricosus").
"""

from pathlib import Path
import subprocess

from loguru import logger
import pyBigWig

# MinIO remote paths
_BGI_SPIDER_ANNO_PREFIX = "rustfs/spider-genome/BGI_RNA_align/Spider_anno"
_ONT_BIGWIG_PREFIX = "rustfs/spider-genome/ONT_RNA_align/bigwig"


def _run_mc(args: list[str]) -> str:
    """Run a mc subcommand, return stdout. Raises on non-zero exit."""
    result = subprocess.run(
        ["mc"] + args,
        capture_output=True,
        text=True,
        timeout=120,
    )
    if result.returncode != 0:
        raise RuntimeError(f"mc {' '.join(args)} failed: {result.stderr.strip()}")
    return result.stdout


def build_bw_index(bw_cache_dir: Path) -> dict[str, dict[str, Path | None]]:
    """
    Discover available BigWig files for all species.
    Queries the MinIO remote to build a {species_id: {BGI: remote_path, ONT: remote_path}} index.

    Returns a dict mapping local species_id → {
        "BGI": remote mc path (str) or None,
        "ONT": remote mc path (str) or None,
    }.
    """
    index: dict[str, dict[str, str | None]] = {}

    # --- BGI: Spider_anno/{NNN}.{species}/outs/bam.bw (102 species) ---
    try:
        listing = _run_mc(["ls", f"{_BGI_SPIDER_ANNO_PREFIX}/"])
        for line in listing.splitlines():
            # mc ls output format: "[date] [time] [size] [name]/"
            parts = line.strip().split()
            if not parts:
                continue
            name = parts[-1].rstrip("/")  # e.g., "001.Allagelena_difficilis"
            if not name or "." not in name:
                continue
            remote_bw = f"{_BGI_SPIDER_ANNO_PREFIX}/{name}/outs/bam.bw"
            index.setdefault(name, {"BGI": None, "ONT": None})["BGI"] = remote_bw
        logger.info(f"BGI Spider_anno index: {len(index)} species found")
    except Exception as exc:
        logger.warning(f"Could not list BGI Spider_anno on MinIO: {exc}")

    # --- ONT: bigwig/{species}.bw (10 species, pure species name without NNN prefix) ---
    try:
        listing = _run_mc(["ls", f"{_ONT_BIGWIG_PREFIX}/"])
        for line in listing.splitlines():
            parts = line.strip().split()
            if not parts:
                continue
            filename = parts[-1]  # e.g., "Araneus_ventricosus.bw"
            if not filename.endswith(".bw"):
                continue
            species_name = filename[:-3]  # strip ".bw" → "Araneus_ventricosus"
            remote_bw = f"{_ONT_BIGWIG_PREFIX}/{filename}"
            # Match to NNN-prefixed species ID in existing index
            matched = False
            for species_id in index:
                if species_id.split(".", 1)[-1] == species_name:
                    index[species_id]["ONT"] = remote_bw
                    matched = True
                    break
            if not matched:
                # ONT species not in BGI index; add with NNN=unknown
                index[species_name] = {"BGI": None, "ONT": remote_bw}
        logger.info(f"ONT bigwig index: entries with ONT data = {sum(1 for v in index.values() if v['ONT'])}")
    except Exception as exc:
        logger.warning(f"Could not list ONT bigwig on MinIO: {exc}")

    return index


def download_bw(
    species_id: str,
    remote_path: str,
    local_path: Path,
) -> bool:
    """
    Download a BigWig file from MinIO if not already cached.
    Returns True on success (or if already cached), False on failure.
    """
    if local_path.exists():
        logger.debug(f"Cache hit: {local_path}")
        return True

    local_path.parent.mkdir(parents=True, exist_ok=True)
    logger.info(f"Downloading {remote_path} → {local_path}")
    try:
        _run_mc(["cp", remote_path, str(local_path)])
        return True
    except Exception as exc:
        logger.warning(f"Failed to download {remote_path}: {exc}")
        return False


def query_bw_signal(
    bw_path: Path,
    chrom: str,
    start: int,
    end: int,
) -> float | None:
    """
    Query mean BigWig signal over [start, end) (0-based half-open coordinates).
    Input locus coordinates are 1-based closed; converted internally.
    Returns None on failure (missing chromosome, file error, etc.).
    """
    # Convert 1-based closed → 0-based half-open
    bw_start = start - 1
    bw_end = end

    try:
        bw = pyBigWig.open(str(bw_path))
        chroms = bw.chroms()
        if chrom not in chroms:
            bw.close()
            return None
        chrom_len = chroms[chrom]
        bw_end = min(bw_end, chrom_len)
        if bw_start >= bw_end:
            bw.close()
            return None
        values = bw.stats(chrom, bw_start, bw_end, type="mean")
        bw.close()
        return values[0] if values and values[0] is not None else 0.0
    except Exception as exc:
        logger.debug(f"BigWig query failed ({bw_path}, {chrom}:{start}-{end}): {exc}")
        return None


def annotate_loci(
    loci: list[dict],
    bw_cache_dir: Path,
    bw_index: dict[str, dict] | None = None,
) -> list[dict]:
    """
    Add RNA-seq support columns to each locus dict.
    Downloads BigWig files on demand; skips species without available data.

    Modifies loci in-place and returns the list.
    Added keys:
      bgi_rna_support (bool | None) — True if mean BGI signal > 0 over locus
      ont_rna_support (bool | None) — True if mean ONT signal > 0 over locus
      rna_mean_signal (float | None) — BGI mean signal (primary)
    """
    if bw_index is None:
        bw_index = build_bw_index(bw_cache_dir)

    # Cache downloaded paths per species
    local_paths: dict[str, dict[str, Path | None]] = {}

    for locus in loci:
        species_id = locus["species"]

        if species_id not in local_paths:
            info = bw_index.get(species_id, {})
            bgi_local = ont_local = None

            if info.get("BGI"):
                bgi_dest = bw_cache_dir / "BGI" / f"{species_id}.bw"
                if download_bw(species_id, info["BGI"], bgi_dest):
                    bgi_local = bgi_dest

            if info.get("ONT"):
                ont_dest = bw_cache_dir / "ONT" / f"{species_id}.bw"
                if download_bw(species_id, info["ONT"], ont_dest):
                    ont_local = ont_dest

            local_paths[species_id] = {"BGI": bgi_local, "ONT": ont_local}

        paths = local_paths[species_id]
        seqid = locus["seqid"]
        start = locus["start"]
        end = locus["end"]

        # BGI signal
        bgi_signal: float | None = None
        if paths["BGI"] is not None:
            bgi_signal = query_bw_signal(paths["BGI"], seqid, start, end)
        locus["bgi_rna_support"] = (bgi_signal > 0) if bgi_signal is not None else None
        locus["rna_mean_signal"] = bgi_signal

        # ONT signal
        ont_signal: float | None = None
        if paths["ONT"] is not None:
            ont_signal = query_bw_signal(paths["ONT"], seqid, start, end)
        locus["ont_rna_support"] = (ont_signal > 0) if ont_signal is not None else None

    return loci
