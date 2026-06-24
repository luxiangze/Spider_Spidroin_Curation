#!/usr/bin/env python3
"""ONT transcript spidroin screening via nucleotide nhmmer.

Search ONT isoform transcripts with nucleotide HMM profiles (NTD/CTD),
filter transcripts containing both terminal domains, and determine the
CDS sequence bounded by start/stop codons.
"""

from __future__ import annotations

from collections import Counter, defaultdict
import csv
from dataclasses import dataclass
from pathlib import Path
import re
import shlex
import shutil
import subprocess
from typing import Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from loguru import logger
import polars as pl
from tqdm import tqdm
import typer

from spider_silkome_module.config import (
    INTERIM_DATA_DIR,
    PROCESSED_DATA_DIR,
    PROJ_ROOT,
    REFERENCES_DIR,
)
from spider_silkome_module.parse_nhmmer_results import parse_tbl_file
from spider_silkome_module.run_nhmmer import press_hmm_models

app = typer.Typer(help="Screen ONT isoform transcripts for spidroins via nhmmer.")

DEFAULT_ISOFORM_DIR = PROJ_ROOT / "workflow" / "ONT_RNA_align" / "results" / "isoforms"
DEFAULT_HMM_DIR = REFERENCES_DIR / "2025_Schoneberg_data" / "hmmer_nucl_profile_trimmed"
DEFAULT_INTERIM_DIR = INTERIM_DATA_DIR / "ont_nhmmer_spidroin"
DEFAULT_PROCESSED_DIR = PROCESSED_DATA_DIR / "ont_nhmmer_spidroin"

DOMAIN_RE = re.compile(r"^(.+)_(NTD|CTD)(?:\.hmm)?$")
STOP_CODONS = {"TAA", "TAG", "TGA"}
START_CODON = "ATG"


@dataclass(frozen=True)
class SpeciesRecord:
    species: str
    transcript_fasta: Path


@dataclass
class DomainHit:
    """A single nhmmer domain hit on a transcript."""

    species: str
    target_name: str
    target_len: int
    query_name: str
    hmm_from: int
    hmm_to: int
    ali_from: int
    ali_to: int
    env_from: int
    env_to: int
    strand: str
    e_value: float
    score: float
    bias: float
    spidroin_type: str
    domain_type: str

    @property
    def hmm_coverage(self) -> float:
        # HMM profile length is unknown per-hit; use ali span vs env span as proxy.
        # filter_hits already enforces coverage via its own hmm_length estimate.
        return 0.0


def parse_profile_name(profile: str) -> tuple[str, str]:
    """Extract (spidroin_type, domain_type) from a profile name like MaSp1_NTD."""
    match = DOMAIN_RE.match(Path(profile).stem)
    if not match:
        return "", ""
    return match.group(1), match.group(2)


def species_from_fasta(path: Path) -> str:
    name = path.name
    for suffix in (".isoforms.fa", ".isoforms.fasta", ".transcripts.fa", ".fa", ".fasta"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def discover_species_records(
    isoform_dir: Path,
    species_filter: Optional[set[str]] = None,
) -> tuple[list[SpeciesRecord], list[dict[str, str]]]:
    """Find ONT transcript FASTA files and build species records."""
    records: list[SpeciesRecord] = []
    rows: list[dict[str, str]] = []
    patterns = ("*.isoforms.fa", "*.isoforms.fasta", "*.transcripts.fa")
    seen: set[Path] = set()
    for pattern in patterns:
        for fasta in sorted(isoform_dir.glob(pattern)):
            if fasta in seen:
                continue
            seen.add(fasta)
            species = species_from_fasta(fasta)
            if species_filter and species not in species_filter:
                continue
            rows.append(
                {
                    "species": species,
                    "transcript_fasta": str(fasta),
                    "status": "matched",
                    "reason": "",
                }
            )
            records.append(SpeciesRecord(species=species, transcript_fasta=fasta))
    return records, rows


def write_tsv(rows: list[dict], output_path: Path, fieldnames: Optional[list[str]] = None) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "" if row.get(key) is None else row.get(key) for key in fieldnames})


def run_command(cmd: list[str], outputs: list[Path], force: bool = False) -> None:
    if outputs and all(p.exists() and p.stat().st_size > 0 for p in outputs) and not force:
        logger.info(f"Skipping existing output: {outputs[0]}")
        return
    for output in outputs:
        output.parent.mkdir(parents=True, exist_ok=True)
    logger.info(shlex.join([str(part) for part in cmd]))
    subprocess.run([str(part) for part in cmd], check=True)


def concatenate_hmms(hmm_dir: Path, output_hmm: Path, force: bool = False) -> list[Path]:
    """Concatenate all *_TD.hmm files into a single multi-model HMM file."""
    hmm_files = sorted(hmm_dir.glob("*TD.hmm"))
    if not hmm_files:
        raise FileNotFoundError(f"No *_TD.hmm files found in {hmm_dir}")
    if output_hmm.exists() and output_hmm.stat().st_size > 0 and not force:
        return hmm_files
    output_hmm.parent.mkdir(parents=True, exist_ok=True)
    with output_hmm.open("wb") as out:
        for hmm_file in hmm_files:
            out.write(hmm_file.read_bytes())
            out.write(b"\n")
    return hmm_files


def run_nhmmer_for_species(
    combined_hmm: Path,
    transcript_fasta: Path,
    tblout: Path,
    textout: Path,
    threads: int,
    evalue: float,
    force: bool = False,
) -> None:
    """Run nhmmer on a transcript FASTA against the combined nucleotide HMM."""
    if tblout.exists() and tblout.stat().st_size > 0 and not force:
        logger.info(f"Skipping existing nhmmer output: {tblout}")
        return
    tblout.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "nhmmer",
        "--cpu",
        str(threads),
        "--noali",
        "-E",
        str(evalue),
        "--tblout",
        tblout,
        combined_hmm,
        transcript_fasta,
    ]
    with textout.open("w") as out:
        subprocess.run([str(c) for c in cmd], stdout=out, check=True)


def filter_hits_per_model(
    df: pl.DataFrame,
    e_value_threshold: float,
    coverage_threshold: float,
) -> pl.DataFrame:
    """Filter nhmmer hits with per-model HMM length.

    Unlike parse_nhmmer_results.filter_hits which uses a single global HMM
    length, this groups hits by query_name so each HMM model's coverage is
    computed against its own profile length. This is required when a single
    tblout contains hits from multiple concatenated HMM models (NTD + CTD)
    of different lengths.
    """
    if df.is_empty():
        return df
    parts: list[pl.DataFrame] = []
    for (query_name,), group in df.group_by("query_name", maintain_order=True):
        hmm_length = group["hmm_to"].max()
        if hmm_length == 0:
            continue
        filtered = group.filter(
            (pl.col("e_value") < e_value_threshold)
            & (((pl.col("hmm_to") - pl.col("hmm_from") + 1) / hmm_length) >= coverage_threshold)
            & (pl.col("score") > pl.col("bias"))
        )
        if not filtered.is_empty():
            parts.append(filtered)
    return pl.concat(parts) if parts else df.head(0)


def parse_hits_for_species(
    tblout: Path,
    species: str,
    e_value_threshold: float,
    coverage_threshold: float,
) -> list[DomainHit]:
    """Parse and filter nhmmer tblout into DomainHit objects."""
    if not tblout.exists() or tblout.stat().st_size == 0:
        return []
    df = parse_tbl_file(tblout, species)
    if df.is_empty():
        return []
    filtered = filter_hits_per_model(df, e_value_threshold, coverage_threshold)
    hits: list[DomainHit] = []
    for row in filtered.iter_rows(named=True):
        spidroin_type, domain_type = parse_profile_name(row["query_name"])
        if not spidroin_type or not domain_type:
            continue
        hits.append(
            DomainHit(
                species=species,
                target_name=row["target_name"],
                target_len=row["sq_len"],
                query_name=row["query_name"],
                hmm_from=row["hmm_from"],
                hmm_to=row["hmm_to"],
                ali_from=row["ali_from"],
                ali_to=row["ali_to"],
                env_from=row["env_from"],
                env_to=row["env_to"],
                strand=row["strand"],
                e_value=row["e_value"],
                score=row["score"],
                bias=row["bias"],
                spidroin_type=spidroin_type,
                domain_type=domain_type,
            )
        )
    return hits


def find_start_codon(seq: str, frame: int, upstream_pos: int) -> Optional[int]:
    """Scan upstream from upstream_pos in the given frame for the first ATG.

    Returns the 0-based start index of the ATG, or None if not found.
    frame is the reading frame (0, 1, or 2) relative to the transcript.
    """
    # Align upstream_pos to frame, then step backwards in codons.
    pos = upstream_pos - ((upstream_pos - frame) % 3)
    while pos >= frame:
        codon = seq[pos:pos + 3]
        if len(codon) < 3:
            break
        if codon == START_CODON:
            return pos
        pos -= 3
    return None


def find_stop_codon(seq: str, frame: int, downstream_pos: int) -> Optional[int]:
    """Scan downstream from downstream_pos in the given frame for the first stop codon.

    Returns the 0-based start index of the stop codon, or None if not found.
    """
    pos = downstream_pos - ((downstream_pos - frame) % 3)
    seq_len = len(seq)
    while pos + 3 <= seq_len:
        codon = seq[pos:pos + 3]
        if codon in STOP_CODONS:
            return pos
        pos += 3
    return None


def determine_cds(
    transcript_seq: str,
    ntd_env_from: int,
    ctd_env_to: int,
) -> tuple[int, int, list[str]]:
    """Determine CDS boundaries from NTD/CTD envelope coordinates.

    env_from/env_to are 1-based inclusive nhmmer coordinates.
    Returns (cds_start_1based, cds_end_1based, flags).
    """
    flags: list[str] = []
    seq = str(transcript_seq)
    seq_len = len(seq)

    # Convert to 0-based for internal processing.
    ntd_start_0 = max(0, ntd_env_from - 1)
    ctd_end_0 = min(seq_len, ctd_env_to)

    frame = ntd_start_0 % 3

    # Find start codon upstream of NTD.
    start_idx = find_start_codon(seq, frame, ntd_start_0)
    if start_idx is None:
        flags.append("missing_start_codon")
        start_idx = ntd_start_0

    # Find stop codon downstream of CTD.
    stop_idx = find_stop_codon(seq, frame, ctd_end_0)
    if stop_idx is None:
        flags.append("missing_stop_codon")
        end_idx = ctd_end_0
    else:
        # CDS ends just before the stop codon (stop codon excluded from CDS).
        end_idx = stop_idx

    # Ensure length is a multiple of 3.
    cds_len = end_idx - start_idx
    if cds_len % 3 != 0:
        end_idx = start_idx + (cds_len // 3) * 3

    return start_idx + 1, end_idx, flags


def classify_and_extract(
    hits: list[DomainHit],
    transcript_seqs: dict[str, str],
    e_value_threshold: float,
) -> tuple[list[dict], list[dict]]:
    """Group hits by transcript, classify, and extract CDS for full-length candidates.

    Returns (result_rows, audit_rows).
    """
    grouped: dict[tuple[str, str], list[DomainHit]] = defaultdict(list)
    for hit in hits:
        grouped[(hit.species, hit.target_name)].append(hit)

    rows: list[dict] = []
    audit_rows: list[dict] = []

    for (species, target_name), target_hits in sorted(grouped.items()):
        ntd_hits = [h for h in target_hits if h.domain_type == "NTD"]
        ctd_hits = [h for h in target_hits if h.domain_type == "CTD"]

        # Strand check: transcripts are spliced mRNA, all hits should be on '+'.
        strand_flags: list[str] = []
        neg_strands = [h for h in target_hits if h.strand == "-"]
        if neg_strands:
            strand_flags.append("strand_mismatch")
            logger.warning(
                f"{species}|{target_name}: {len(neg_strands)} hits on '-' strand "
                f"(expected '+' for transcripts)"
            )

        # Select best NTD (lowest e_value) and best CTD.
        best_ntd = sorted(ntd_hits, key=lambda h: h.e_value)[0] if ntd_hits else None
        best_ctd = sorted(ctd_hits, key=lambda h: h.e_value)[0] if ctd_hits else None

        # Classification.
        if best_ntd and best_ctd:
            classification = "full_length"
        elif best_ntd:
            classification = "ntd_only"
        elif best_ctd:
            classification = "ctd_only"
        else:
            continue

        # Spidroin type confirmation.
        type_flags: list[str] = []
        if best_ntd and best_ctd:
            if best_ntd.spidroin_type == best_ctd.spidroin_type:
                spidroin_type = best_ntd.spidroin_type
            else:
                type_flags.append("type_mismatch")
                spidroin_type = f"{best_ntd.spidroin_type}/{best_ctd.spidroin_type}"
        elif best_ntd:
            spidroin_type = best_ntd.spidroin_type
        else:
            spidroin_type = best_ctd.spidroin_type if best_ctd else "Unknown"

        # CDS determination for full_length candidates.
        cds_start = ""
        cds_end = ""
        cds_len = ""
        cds_flags: list[str] = []
        if classification == "full_length":
            seq = transcript_seqs.get(target_name, "")
            if seq:
                cds_start_1, cds_end_1, cds_flags = determine_cds(
                    seq,
                    ntd_env_from=best_ntd.env_from,
                    ctd_env_to=best_ctd.env_to,
                )
                cds_start = cds_start_1
                cds_end = cds_end_1
                cds_len = cds_end_1 - cds_start_1 + 1
            else:
                cds_flags.append("missing_transcript_sequence")

        all_flags = strand_flags + type_flags + cds_flags
        notes = ",".join(all_flags) if all_flags else ""

        row = {
            "species": species,
            "target_id": target_name,
            "transcript_len": target_hits[0].target_len,
            "classification": classification,
            "spidroin_type": spidroin_type,
            "ntd_profile": best_ntd.query_name if best_ntd else "",
            "ntd_spidroin_type": best_ntd.spidroin_type if best_ntd else "",
            "ntd_evalue": best_ntd.e_value if best_ntd else None,
            "ntd_score": best_ntd.score if best_ntd else None,
            "ntd_env_from": best_ntd.env_from if best_ntd else None,
            "ntd_env_to": best_ntd.env_to if best_ntd else None,
            "ntd_strand": best_ntd.strand if best_ntd else "",
            "ctd_profile": best_ctd.query_name if best_ctd else "",
            "ctd_spidroin_type": best_ctd.spidroin_type if best_ctd else "",
            "ctd_evalue": best_ctd.e_value if best_ctd else None,
            "ctd_score": best_ctd.score if best_ctd else None,
            "ctd_env_from": best_ctd.env_from if best_ctd else None,
            "ctd_env_to": best_ctd.env_to if best_ctd else None,
            "ctd_strand": best_ctd.strand if best_ctd else "",
            "cds_start": cds_start,
            "cds_end": cds_end,
            "cds_len": cds_len,
            "num_ntd_hits": len(ntd_hits),
            "num_ctd_hits": len(ctd_hits),
            "notes": notes,
        }
        rows.append(row)

        if classification == "full_length" and best_ntd and best_ctd:
            audit_rows.append(
                {
                    "species": species,
                    "target_id": target_name,
                    "spidroin_type": spidroin_type,
                    "ntd_profile": best_ntd.query_name,
                    "ctd_profile": best_ctd.query_name,
                    "ntd_evalue": best_ntd.e_value,
                    "ctd_evalue": best_ctd.e_value,
                    "cds_start": cds_start,
                    "cds_end": cds_end,
                    "cds_len": cds_len,
                    "notes": notes,
                }
            )

    return rows, audit_rows


RESULT_COLUMNS = [
    "species",
    "target_id",
    "transcript_len",
    "classification",
    "spidroin_type",
    "ntd_profile",
    "ntd_spidroin_type",
    "ntd_evalue",
    "ntd_score",
    "ntd_env_from",
    "ntd_env_to",
    "ntd_strand",
    "ctd_profile",
    "ctd_spidroin_type",
    "ctd_evalue",
    "ctd_score",
    "ctd_env_from",
    "ctd_env_to",
    "ctd_strand",
    "cds_start",
    "cds_end",
    "cds_len",
    "num_ntd_hits",
    "num_ctd_hits",
    "notes",
]


def load_transcripts(fasta_path: Path) -> dict[str, str]:
    """Load transcript sequences as {id: sequence_string}."""
    return {rec.id: str(rec.seq) for rec in SeqIO.parse(fasta_path, "fasta")}


def write_cds_fastas(
    full_rows: list[dict],
    transcript_seqs: dict[str, str],
    processed_dir: Path,
) -> int:
    """Write CDS nucleotide and translated protein FASTAs for full-length candidates."""
    cds_records: list[SeqRecord] = []
    protein_records: list[SeqRecord] = []
    by_type: dict[str, list[SeqRecord]] = defaultdict(list)

    for row in full_rows:
        target_id = row["target_id"]
        seq = transcript_seqs.get(target_id)
        if not seq or not row["cds_start"] or not row["cds_end"]:
            continue
        start = int(row["cds_start"]) - 1  # convert to 0-based
        end = int(row["cds_end"])
        cds_seq = seq[start:end]
        if len(cds_seq) == 0 or len(cds_seq) % 3 != 0:
            continue

        cds_rec = SeqRecord(
            Seq(cds_seq),
            id=target_id,
            description=(
                f"spidroin_type={row['spidroin_type']} "
                f"cds={row['cds_start']}-{row['cds_end']} "
                f"len={len(cds_seq)} notes={row.get('notes', '')}"
            ),
        )
        cds_records.append(cds_rec)
        by_type[row["spidroin_type"]].append(cds_rec)

        # Translate to protein.
        protein_seq = Seq(cds_seq).translate()
        protein_str = str(protein_seq)
        if protein_str.endswith("*"):
            protein_str = protein_str[:-1]
        protein_str = protein_str.replace("*", "X")
        protein_rec = SeqRecord(
            Seq(protein_str),
            id=target_id,
            description=f"spidroin_type={row['spidroin_type']} cds_len={len(cds_seq)}",
        )
        protein_records.append(protein_rec)

    processed_dir.mkdir(parents=True, exist_ok=True)
    SeqIO.write(cds_records, processed_dir / "full_length_spidroins_cds.fasta", "fasta")
    SeqIO.write(protein_records, processed_dir / "full_length_spidroins_protein.fasta", "fasta")

    by_type_dir = processed_dir / "by_type"
    by_type_dir.mkdir(parents=True, exist_ok=True)
    for spidroin_type, type_records in sorted(by_type.items()):
        safe_type = re.sub(r"[^A-Za-z0-9_.-]+", "_", spidroin_type)
        SeqIO.write(type_records, by_type_dir / f"{safe_type}.cds.fasta", "fasta")

    return len(cds_records)


def write_pivot(full_rows: list[dict], species_names: list[str], output_path: Path) -> None:
    counts: dict[tuple[str, str], int] = Counter(
        (row["species"], row["spidroin_type"]) for row in full_rows
    )
    types = sorted({row["spidroin_type"] for row in full_rows})
    fieldnames = ["species", *types]
    rows = [
        {"species": species, **{t: counts[(species, t)] for t in types}}
        for species in sorted(species_names)
    ]
    write_tsv(rows, output_path, fieldnames=fieldnames)


def validate_processed_outputs(processed_dir: Path) -> list[dict[str, str]]:
    full_tsv = processed_dir / "full_length_spidroins.tsv"
    cds_fasta = processed_dir / "full_length_spidroins_cds.fasta"
    checks: list[dict[str, str]] = []

    full_rows = []
    if full_tsv.exists():
        with full_tsv.open() as handle:
            full_rows = list(csv.DictReader(handle, delimiter="\t"))

    cds_ids = [rec.id for rec in SeqIO.parse(cds_fasta, "fasta")] if cds_fasta.exists() else []
    cds_lens = [len(rec.seq) for rec in SeqIO.parse(cds_fasta, "fasta")] if cds_fasta.exists() else []
    non_multiple_of_3 = sum(1 for length in cds_lens if length % 3 != 0)
    same_type_failures = [
        row["target_id"]
        for row in full_rows
        if row.get("ntd_spidroin_type") and row.get("ctd_spidroin_type")
        and row.get("ntd_spidroin_type") != row.get("ctd_spidroin_type")
    ]
    count_matches = len(full_rows) == len(cds_ids)

    checks.extend(
        [
            {
                "check": "cds_fasta_count_matches_tsv",
                "status": str(count_matches),
                "value": f"{len(cds_ids)}/{len(full_rows)}",
            },
            {
                "check": "cds_length_multiple_of_3",
                "status": str(non_multiple_of_3 == 0),
                "value": str(non_multiple_of_3),
            },
            {
                "check": "full_length_same_type_ntd_ctd",
                "status": str(not same_type_failures),
                "value": str(len(same_type_failures)),
            },
        ]
    )
    write_tsv(checks, processed_dir / "validation_summary.tsv")
    return checks


def run_classification(
    records: list[SpeciesRecord],
    hmmer_dir: Path,
    processed_dir: Path,
    e_value: float,
    coverage_threshold: float,
) -> None:
    """Parse nhmmer hits, classify, and write outputs."""
    all_hits: list[DomainHit] = []
    transcript_seqs: dict[str, str] = {}

    for record in records:
        tblout = hmmer_dir / f"{record.species}.tblout"
        hits = parse_hits_for_species(tblout, record.species, e_value, coverage_threshold)
        all_hits.extend(hits)
        transcript_seqs.update(load_transcripts(record.transcript_fasta))

    rows, audit_rows = classify_and_extract(all_hits, transcript_seqs, e_value)
    full_rows = [row for row in rows if row["classification"] == "full_length"]
    partial_rows = [row for row in rows if row["classification"] != "full_length"]

    processed_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(rows, processed_dir / "target_classification.tsv", RESULT_COLUMNS)
    write_tsv(full_rows, processed_dir / "full_length_spidroins.tsv", RESULT_COLUMNS)
    write_tsv(partial_rows, processed_dir / "partial_spidroins.tsv", RESULT_COLUMNS)
    write_tsv(
        audit_rows,
        processed_dir / "full_length_pair_audit.tsv",
        [
            "species", "target_id", "spidroin_type", "ntd_profile", "ctd_profile",
            "ntd_evalue", "ctd_evalue", "cds_start", "cds_end", "cds_len", "notes",
        ],
    )
    write_pivot(full_rows, [r.species for r in records], processed_dir / "spidroin_pivot_table.tsv")
    fasta_count = write_cds_fastas(full_rows, transcript_seqs, processed_dir)
    validate_processed_outputs(processed_dir)
    logger.success(
        f"Classified {len(rows)} transcripts: {len(full_rows)} full-length, "
        f"{len(partial_rows)} partial; wrote {fasta_count} CDS FASTA records."
    )


def preflight_checks(isoform_dir: Path, hmm_dir: Path) -> list[dict[str, str]]:
    checks: list[dict[str, str]] = []
    for tool in ("nhmmer", "hmmpress"):
        path = shutil.which(tool)
        checks.append({"check": f"executable:{tool}", "status": str(path is not None), "value": path or ""})
    isoform_count = 0
    for pattern in ("*.isoforms.fa", "*.isoforms.fasta", "*.transcripts.fa"):
        isoform_count += len(list(isoform_dir.glob(pattern)))
    hmm_count = len(list(hmm_dir.glob("*TD.hmm")))
    checks.extend(
        [
            {"check": "isoform_dir_exists", "status": str(isoform_dir.exists()), "value": str(isoform_dir)},
            {"check": "hmm_dir_exists", "status": str(hmm_dir.exists()), "value": str(hmm_dir)},
            {"check": "isoform_fasta_count", "status": str(isoform_count > 0), "value": str(isoform_count)},
            {"check": "hmm_count", "status": str(hmm_count > 0), "value": str(hmm_count)},
        ]
    )
    return checks


@app.command()
def preflight(
    isoform_dir: Path = typer.Option(DEFAULT_ISOFORM_DIR, help="ONT transcript FASTA directory."),
    hmm_dir: Path = typer.Option(DEFAULT_HMM_DIR, help="Nucleotide HMM profile directory."),
    output: Optional[Path] = typer.Option(None, help="Optional TSV output path."),
) -> None:
    """Check required tools and input directories."""
    checks = preflight_checks(isoform_dir, hmm_dir)
    if output:
        write_tsv(checks, output, ["check", "status", "value"])
    for row in checks:
        logger.info(f"{row['check']}: {row['status']} {row['value']}")
    if not all(row["status"] == "True" for row in checks):
        raise typer.Exit(code=1)


@app.command()
def discover(
    isoform_dir: Path = typer.Option(DEFAULT_ISOFORM_DIR, help="ONT transcript FASTA directory."),
    output: Path = typer.Option(DEFAULT_INTERIM_DIR / "species_manifest.tsv", help="Manifest TSV path."),
    species: Optional[list[str]] = typer.Option(None, "--species", help="Restrict to one species; repeatable."),
) -> None:
    """Discover ONT transcript FASTA inputs."""
    _, rows = discover_species_records(isoform_dir, set(species or []) or None)
    write_tsv(rows, output, ["species", "transcript_fasta", "status", "reason"])
    matched = sum(1 for row in rows if row["status"] == "matched")
    logger.success(f"Wrote {output} with {matched}/{len(rows)} matched species")


@app.command()
def run(
    isoform_dir: Path = typer.Option(DEFAULT_ISOFORM_DIR, help="ONT transcript FASTA directory."),
    hmm_dir: Path = typer.Option(DEFAULT_HMM_DIR, help="Nucleotide HMM profile directory."),
    interim_dir: Path = typer.Option(DEFAULT_INTERIM_DIR, help="Intermediate output directory."),
    processed_dir: Path = typer.Option(DEFAULT_PROCESSED_DIR, help="Processed output directory."),
    species: Optional[list[str]] = typer.Option(None, "--species", help="Restrict to one species; repeatable."),
    threads: int = typer.Option(70, "--threads", min=1, help="Threads for nhmmer."),
    evalue: float = typer.Option(1e-10, "--evalue", help="nhmmer E-value threshold."),
    coverage_threshold: float = typer.Option(0.90, "--coverage", help="Minimum HMM profile coverage."),
    force: bool = typer.Option(False, "--force", help="Regenerate existing outputs."),
) -> None:
    """Run the complete ONT nhmmer spidroin screening workflow."""
    checks = preflight_checks(isoform_dir, hmm_dir)
    failed = [row for row in checks if row["status"] != "True"]
    if failed:
        for row in failed:
            logger.error(f"Preflight failed: {row}")
        raise typer.Exit(code=1)

    records, manifest_rows = discover_species_records(isoform_dir, set(species or []) or None)
    if not records:
        raise typer.BadParameter("No matched species found.")

    interim_dir.mkdir(parents=True, exist_ok=True)
    processed_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(manifest_rows, interim_dir / "species_manifest.tsv", ["species", "transcript_fasta", "status", "reason"])

    hmmer_dir = interim_dir / "hmmer"
    combined_hmm = hmmer_dir / "spidroin_nucl_terminals.all.hmm"

    logger.info("Pressing HMM models...")
    press_hmm_models(hmm_dir, force=force)
    concatenate_hmms(hmm_dir, combined_hmm, force=force)

    for record in tqdm(records, desc="Species"):
        tblout = hmmer_dir / f"{record.species}.tblout"
        textout = hmmer_dir / f"{record.species}.nhmmer.txt"
        run_nhmmer_for_species(
            combined_hmm,
            record.transcript_fasta,
            tblout,
            textout,
            threads=threads,
            evalue=evalue,
            force=force,
        )

    run_classification(records, hmmer_dir, processed_dir, e_value=evalue, coverage_threshold=coverage_threshold)


if __name__ == "__main__":
    app()
