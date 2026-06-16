from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
from typing import Any

from loguru import logger

from spider_silkome_module.config import PROCESSED_DATA_DIR, RAW_DATA_DIR

FASTA_PATTERNS = ("*.fa", "*.fna", "*.fasta", "*.fa.gz", "*.fna.gz", "*.fasta.gz")
NON_GENOME_KEYWORDS = ("pep", "prot", "mito", "cds", "transcript")


@dataclass(frozen=True)
class Interval:
    chrom: str
    start0: int
    end0: int
    label: str

    @property
    def start1(self) -> int:
        return self.start0 + 1


@dataclass(frozen=True)
class Window:
    chrom: str
    start0: int
    end0: int
    labels: str = "."
    n_raw_intervals: int = 0

    @property
    def start1(self) -> int:
        return self.start0 + 1

    @property
    def window_id(self) -> str:
        return f"{self.chrom}_{self.start1}-{self.end0}:."

    @property
    def region(self) -> str:
        return f"{self.chrom}:{self.start1}-{self.end0}"


def parse_attrs(attr_str: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in attr_str.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value
    return attrs


def parse_float(value: Any) -> float | None:
    if value in (None, "", "."):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def parse_int(value: Any) -> int | None:
    if value in (None, "", "."):
        return None
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def read_fai(fai_path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    with fai_path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                lengths[fields[0]] = int(fields[1])
    return lengths


def _is_genome_candidate(path: Path) -> bool:
    name = path.name.lower()
    return (
        any(path.match(pattern) for pattern in FASTA_PATTERNS)
        and not any(keyword in name for keyword in NON_GENOME_KEYWORDS)
    )


def find_genome_fasta(species: str, raw_ref_dir: Path = RAW_DATA_DIR / "01.ref_gff") -> Path:
    species_dir = raw_ref_dir / species
    if not species_dir.exists():
        raise FileNotFoundError(f"Missing species raw directory: {species_dir}")
    candidates = [path for path in species_dir.iterdir() if path.is_file() and _is_genome_candidate(path)]
    if not candidates:
        raise FileNotFoundError(f"No genome FASTA found in {species_dir}")
    preferred = [path for path in candidates if "genome" in path.name.lower()]
    chosen = sorted(preferred or candidates)[0]
    logger.debug(f"Selected genome FASTA for {species}: {chosen}")
    return chosen


def discover_species(
    raw_ref_dir: Path = RAW_DATA_DIR / "01.ref_gff",
    miniprot_evidence_dir: Path = PROCESSED_DATA_DIR / "miniprot_output",
) -> list[str]:
    raw_species = {path.name for path in raw_ref_dir.iterdir() if path.is_dir()}
    miniprot_species = {path.stem for path in miniprot_evidence_dir.glob("*.gff")}
    return sorted(raw_species & miniprot_species)


def load_windows_bed(path: Path) -> list[Window]:
    windows: list[Window] = []
    if not path.exists() or path.stat().st_size == 0:
        return windows
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            labels = fields[3] if len(fields) > 3 else "."
            n_raw = parse_int(fields[4]) if len(fields) > 4 else 0
            windows.append(Window(fields[0], int(fields[1]), int(fields[2]), labels, n_raw or 0))
    return windows


def overlaps(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return a_end >= b_start and a_start <= b_end

