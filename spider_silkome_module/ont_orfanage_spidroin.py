#!/usr/bin/env python3
"""ONT isoform ORFanage protein annotation and spidroin screening."""

from __future__ import annotations

from collections import Counter, defaultdict
import csv
from dataclasses import dataclass
from pathlib import Path
import re
import shlex
import shutil
import subprocess
from typing import Iterable, Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from loguru import logger
from tqdm import tqdm
import typer

from spider_silkome_module.config import (
    INTERIM_DATA_DIR,
    PROCESSED_DATA_DIR,
    PROJ_ROOT,
    RAW_DATA_DIR,
)

app = typer.Typer(help="Annotate ONT isoform ORFs with ORFanage and screen spidroins.")

DEFAULT_ISOFORM_DIR = PROJ_ROOT / "workflow" / "ONT_RNA_align" / "results" / "isoforms"
DEFAULT_GENOME_DIR = RAW_DATA_DIR / "spider_genome"
DEFAULT_HMM_DIR = INTERIM_DATA_DIR / "spider_silkome_20251222" / "hmmbuild_output"
DEFAULT_INTERIM_DIR = INTERIM_DATA_DIR / "ont_orfanage_spidroin"
DEFAULT_PROCESSED_DIR = PROCESSED_DATA_DIR / "ont_orfanage_spidroin"

DOMAIN_RE = re.compile(r"^(.+)_(NTD|CTD)(?:\.hmm)?$")


@dataclass(frozen=True)
class SpeciesRecord:
    species: str
    isoform_gtf: Path
    genome_fasta: Path
    template_gff: Path


@dataclass
class DomainHit:
    species: str
    target_id: str
    tlen: int
    query_id: str
    qlen: int
    full_evalue: float
    full_score: float
    full_bias: float
    dom_index: int
    dom_count: int
    c_evalue: float
    i_evalue: float
    dom_score: float
    dom_bias: float
    hmm_from: int
    hmm_to: int
    ali_from: int
    ali_to: int
    env_from: int
    env_to: int
    acc: float
    spidroin_type: str
    domain_type: str

    @property
    def hmm_coverage(self) -> float:
        return (abs(self.hmm_to - self.hmm_from) + 1) / self.qlen if self.qlen else 0.0

    @property
    def env_coverage(self) -> float:
        return (abs(self.env_to - self.env_from) + 1) / self.qlen if self.qlen else 0.0

    def is_terminal(self, terminal_window: int) -> bool:
        if self.domain_type == "NTD":
            return self.env_from <= terminal_window
        if self.domain_type == "CTD":
            return self.env_to >= max(1, self.tlen - terminal_window + 1)
        return False


def parse_attributes(raw: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in raw.strip().rstrip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key.strip()] = value.strip().strip('"')
        elif " " in item:
            key, value = item.split(" ", 1)
            attrs[key.strip()] = value.strip().strip('"')
    return attrs


def parse_profile_name(profile: str) -> tuple[str, str]:
    profile = Path(profile).stem
    match = DOMAIN_RE.match(profile)
    if not match:
        return "", ""
    return match.group(1), match.group(2)


def species_from_isoform_gtf(path: Path) -> str:
    name = path.name
    suffix = ".isoforms.gtf"
    if not name.endswith(suffix):
        raise ValueError(f"Unexpected isoform GTF name: {path}")
    return name[: -len(suffix)]


def first_existing(paths: Iterable[Path]) -> Optional[Path]:
    for path in paths:
        if path.exists():
            return path
    return None


def find_genome_fasta(genome_dir: Path, species: str) -> Optional[Path]:
    candidates = [
        *genome_dir.glob(f"*/{species}.fa"),
        *genome_dir.glob(f"*/{species}.fasta"),
        *genome_dir.glob(f"*/{species}.fna"),
    ]
    candidates = [
        p
        for p in candidates
        if not any(token in p.name.lower() for token in ("protein", "pep", "cds", "transcript"))
    ]
    return sorted(candidates)[0] if candidates else None


def find_template_gff(genome_dir: Path, species: str) -> Optional[Path]:
    candidates = [
        *genome_dir.glob(f"*/{species}.gff"),
        *genome_dir.glob(f"*/{species}.gff3"),
        *genome_dir.glob(f"*/{species}.gtf"),
    ]
    if not candidates:
        return None
    exact = [p for p in candidates if p.name == f"{species}.gff"]
    if exact:
        return sorted(exact)[0]
    non_fixed = [p for p in candidates if "fixed" not in p.name.lower()]
    return sorted(non_fixed or candidates)[0]


def discover_species_records(
    isoform_dir: Path,
    genome_dir: Path,
    species_filter: Optional[set[str]] = None,
) -> tuple[list[SpeciesRecord], list[dict[str, str]]]:
    rows: list[dict[str, str]] = []
    records: list[SpeciesRecord] = []

    for isoform_gtf in sorted(isoform_dir.glob("*.isoforms.gtf")):
        species = species_from_isoform_gtf(isoform_gtf)
        if species_filter and species not in species_filter:
            continue

        genome_fasta = find_genome_fasta(genome_dir, species)
        template_gff = find_template_gff(genome_dir, species)
        status = "matched"
        reason = ""
        if genome_fasta is None:
            status = "skipped"
            reason = "missing genome FASTA"
        elif template_gff is None:
            status = "skipped"
            reason = "missing template GFF/GTF"

        row = {
            "species": species,
            "isoform_gtf": str(isoform_gtf),
            "genome_fasta": str(genome_fasta or ""),
            "template_gff": str(template_gff or ""),
            "status": status,
            "reason": reason,
        }
        rows.append(row)

        if status == "matched" and genome_fasta and template_gff:
            records.append(
                SpeciesRecord(
                    species=species,
                    isoform_gtf=isoform_gtf,
                    genome_fasta=genome_fasta,
                    template_gff=template_gff,
                )
            )

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
    if outputs and all(path.exists() and path.stat().st_size > 0 for path in outputs) and not force:
        logger.info(f"Skipping existing output: {outputs[0]}")
        return
    for output in outputs:
        output.parent.mkdir(parents=True, exist_ok=True)
    logger.info(shlex.join([str(part) for part in cmd]))
    subprocess.run([str(part) for part in cmd], check=True)


def count_feature(gtf_path: Path, feature_type: str) -> int:
    if not gtf_path.exists():
        return 0
    count = 0
    with gtf_path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[2] == feature_type:
                count += 1
    return count


def record_id_from_attrs(attrs: dict[str, str]) -> Optional[str]:
    return attrs.get("ID") or attrs.get("transcript_id")


def parent_ids_from_attrs(attrs: dict[str, str]) -> list[str]:
    parent = attrs.get("Parent") or attrs.get("transcript_id")
    if not parent:
        return []
    return [item.strip() for item in parent.split(",") if item.strip()]


def clean_template_gff(input_gff: Path, output_gff: Path, force: bool = False) -> tuple[int, int]:
    """Write a template GFF with invalid transcript coordinates repaired.

    Some Braker3 GFF files contain mRNA rows where end < start while child
    CDS/exon/UTR records are valid. ORFanage rejects those records before it
    can use the valid CDS evidence. For transcript-like rows, repair the span
    from child feature coordinates. For other rare malformed rows, swap start
    and end so the file remains parseable.
    """
    if output_gff.exists() and output_gff.stat().st_size > 0 and not force:
        return 0, 0

    output_gff.parent.mkdir(parents=True, exist_ok=True)
    raw_lines = input_gff.read_text().splitlines()
    child_spans: dict[str, list[int]] = {}

    for line in raw_lines:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 9:
            continue
        try:
            start = int(parts[3])
            end = int(parts[4])
        except ValueError:
            continue
        if end < start:
            continue
        attrs = parse_attributes(parts[8])
        for parent_id in parent_ids_from_attrs(attrs):
            span = child_spans.setdefault(parent_id, [start, end])
            span[0] = min(span[0], start)
            span[1] = max(span[1], end)

    repaired_transcripts = 0
    swapped_features = 0
    cleaned_lines: list[str] = []
    for line in raw_lines:
        if line.startswith("#") or not line.strip():
            cleaned_lines.append(line)
            continue
        parts = line.split("\t")
        if len(parts) < 9:
            cleaned_lines.append(line)
            continue
        try:
            start = int(parts[3])
            end = int(parts[4])
        except ValueError:
            cleaned_lines.append(line)
            continue
        if end < start:
            attrs = parse_attributes(parts[8])
            feature_id = record_id_from_attrs(attrs)
            if parts[2] in {"mRNA", "transcript"} and feature_id in child_spans:
                start, end = child_spans[feature_id]
                repaired_transcripts += 1
            else:
                start, end = end, start
                swapped_features += 1
            parts[3] = str(start)
            parts[4] = str(end)
            line = "\t".join(parts)
        cleaned_lines.append(line)

    output_gff.write_text("\n".join(cleaned_lines) + "\n")
    return repaired_transcripts, swapped_features


def parse_gtf_metadata(gtf_path: Path) -> dict[str, dict[str, str]]:
    metadata: dict[str, dict[str, str]] = {}
    with gtf_path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            attrs = parse_attributes(parts[8])
            transcript_id = attrs.get("transcript_id") or attrs.get("ID") or attrs.get("Parent")
            if not transcript_id:
                continue
            entry = metadata.setdefault(transcript_id, {})
            if parts[2] in {"transcript", "mRNA"}:
                entry.update(
                    {
                        "scaffold": parts[0],
                        "source": parts[1],
                        "strand": parts[6],
                        "start": parts[3],
                        "end": parts[4],
                    }
                )
            if attrs.get("gene_id"):
                entry["gene_id"] = attrs["gene_id"]
            elif attrs.get("Parent") and "gene_id" not in entry:
                entry["gene_id"] = attrs["Parent"]
            for key in ("orfanage_status", "orfanage_template"):
                if key in attrs:
                    entry[key] = attrs[key]
    return metadata


def sanitize_protein_sequence(seq: Seq) -> Seq:
    text = str(seq).strip()
    if text.endswith("*"):
        text = text[:-1]
    text = text.replace("*", "X")
    return Seq(text)


def normalize_protein_fasta(raw_fasta: Path, output_fasta: Path, species: str, gtf_path: Path) -> int:
    metadata = parse_gtf_metadata(gtf_path)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    records: list[SeqRecord] = []
    for record in SeqIO.parse(raw_fasta, "fasta"):
        original_id = record.id.split()[0]
        normalized_id = f"{species}|{original_id}"
        meta = metadata.get(original_id, {})
        description_parts = [
            f"original_id={original_id}",
            f"gene_id={meta.get('gene_id', '')}",
            f"scaffold={meta.get('scaffold', '')}",
            f"strand={meta.get('strand', '')}",
            f"start={meta.get('start', '')}",
            f"end={meta.get('end', '')}",
        ]
        records.append(
            SeqRecord(
                sanitize_protein_sequence(record.seq),
                id=normalized_id,
                description=" ".join(description_parts),
            )
        )

    SeqIO.write(records, output_fasta, "fasta")
    return len(records)


def concatenate_hmms(hmm_dir: Path, output_hmm: Path, force: bool = False) -> list[Path]:
    hmm_files = sorted(hmm_dir.glob("*.hmm"))
    if not hmm_files:
        raise FileNotFoundError(f"No .hmm files found in {hmm_dir}")
    if output_hmm.exists() and output_hmm.stat().st_size > 0 and not force:
        return hmm_files
    output_hmm.parent.mkdir(parents=True, exist_ok=True)
    with output_hmm.open("wb") as out:
        for hmm_file in hmm_files:
            out.write(hmm_file.read_bytes())
            out.write(b"\n")
    return hmm_files


def run_orfanage_for_species(
    record: SpeciesRecord,
    output_gtf: Path,
    stats_tsv: Path,
    threads: int,
    force: bool = False,
) -> None:
    cmd = [
        "orfanage",
        "--reference",
        record.genome_fasta,
        "--query",
        record.isoform_gtf,
        "--output",
        output_gtf,
        "--stats",
        stats_tsv,
        "--threads",
        str(threads),
        "--use_id",
        "--cleanq",
        "--cleant",
        "--rescue",
        record.template_gff,
    ]
    run_command(cmd, [output_gtf, stats_tsv], force=force)


def translate_for_species(
    genome_fasta: Path,
    orfanage_gtf: Path,
    raw_faa: Path,
    normalized_faa: Path,
    species: str,
    force: bool = False,
) -> int:
    run_command(
        ["gffread", orfanage_gtf, "-g", genome_fasta, "-y", raw_faa],
        [raw_faa],
        force=force,
    )
    if normalized_faa.exists() and normalized_faa.stat().st_size > 0 and not force:
        return sum(1 for _ in SeqIO.parse(normalized_faa, "fasta"))
    return normalize_protein_fasta(raw_faa, normalized_faa, species, orfanage_gtf)


def run_hmmer_for_species(
    combined_hmm: Path,
    protein_faa: Path,
    tblout: Path,
    domtblout: Path,
    textout: Path,
    threads: int,
    evalue: float,
    force: bool = False,
) -> None:
    if not any(SeqIO.parse(protein_faa, "fasta")):
        logger.warning(f"No protein sequences in {protein_faa}; skipping HMMER")
        domtblout.write_text("")
        tblout.write_text("")
        textout.write_text("")
        return
    cmd = [
        "hmmsearch",
        "--cpu",
        str(threads),
        "--noali",
        "-E",
        str(evalue),
        "--domE",
        str(evalue),
        "--incE",
        str(evalue),
        "--incdomE",
        str(evalue),
        "-o",
        textout,
        "--tblout",
        tblout,
        "--domtblout",
        domtblout,
        combined_hmm,
        protein_faa,
    ]
    run_command(cmd, [domtblout], force=force)


def parse_domtblout(
    domtblout: Path,
    species: str,
    c_evalue_threshold: float,
    min_hmm_coverage: float,
    require_score_gt_bias: bool = True,
) -> list[DomainHit]:
    hits: list[DomainHit] = []
    if not domtblout.exists() or domtblout.stat().st_size == 0:
        return hits
    with domtblout.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) < 22:
                continue
            spidroin_type, domain_type = parse_profile_name(parts[3])
            if not spidroin_type or not domain_type:
                logger.debug(f"Skipping unrecognized HMM profile: {parts[3]}")
                continue
            try:
                hit = DomainHit(
                    species=species,
                    target_id=parts[0],
                    tlen=int(parts[2]),
                    query_id=parts[3],
                    qlen=int(parts[5]),
                    full_evalue=float(parts[6]),
                    full_score=float(parts[7]),
                    full_bias=float(parts[8]),
                    dom_index=int(parts[9]),
                    dom_count=int(parts[10]),
                    c_evalue=float(parts[11]),
                    i_evalue=float(parts[12]),
                    dom_score=float(parts[13]),
                    dom_bias=float(parts[14]),
                    hmm_from=int(parts[15]),
                    hmm_to=int(parts[16]),
                    ali_from=int(parts[17]),
                    ali_to=int(parts[18]),
                    env_from=int(parts[19]),
                    env_to=int(parts[20]),
                    acc=float(parts[21]),
                    spidroin_type=spidroin_type,
                    domain_type=domain_type,
                )
            except ValueError:
                continue
            if hit.c_evalue > c_evalue_threshold:
                continue
            if require_score_gt_bias and hit.dom_score <= hit.dom_bias:
                continue
            if hit.hmm_coverage < min_hmm_coverage:
                continue
            hits.append(hit)
    return hits


def hit_sort_key(hit: DomainHit) -> tuple[float, float, float]:
    return (-hit.dom_score, hit.c_evalue, -hit.hmm_coverage)


def best_hit(hits: list[DomainHit], terminal_window: int) -> Optional[DomainHit]:
    if not hits:
        return None
    terminal_hits = [hit for hit in hits if hit.is_terminal(terminal_window)]
    return sorted(terminal_hits or hits, key=hit_sort_key)[0]


def split_target_id(target_id: str) -> tuple[str, str]:
    if "|" in target_id:
        species, transcript_id = target_id.split("|", 1)
        return species, transcript_id
    return "", target_id


def hit_columns(prefix: str, hit: Optional[DomainHit], terminal_window: int) -> dict:
    if hit is None:
        return {
            f"{prefix}_profile": "",
            f"{prefix}_spidroin_type": "",
            f"{prefix}_end_ok": False,
            f"{prefix}_c_evalue": None,
            f"{prefix}_dom_score": None,
            f"{prefix}_env_from": None,
            f"{prefix}_env_to": None,
            f"{prefix}_hmm_from": None,
            f"{prefix}_hmm_to": None,
            f"{prefix}_hmm_coverage": None,
            f"{prefix}_env_coverage": None,
        }
    return {
        f"{prefix}_profile": hit.query_id,
        f"{prefix}_spidroin_type": hit.spidroin_type,
        f"{prefix}_end_ok": hit.is_terminal(terminal_window),
        f"{prefix}_c_evalue": hit.c_evalue,
        f"{prefix}_dom_score": hit.dom_score,
        f"{prefix}_env_from": hit.env_from,
        f"{prefix}_env_to": hit.env_to,
        f"{prefix}_hmm_from": hit.hmm_from,
        f"{prefix}_hmm_to": hit.hmm_to,
        f"{prefix}_hmm_coverage": hit.hmm_coverage,
        f"{prefix}_env_coverage": hit.env_coverage,
    }


def load_metadata_for_species(records: list[SpeciesRecord], orfanage_dir: Path) -> dict[str, dict]:
    all_metadata: dict[str, dict] = {}
    for record in records:
        gtf_path = orfanage_dir / f"{record.species}.orfanage.gtf"
        if not gtf_path.exists():
            continue
        metadata = parse_gtf_metadata(gtf_path)
        for transcript_id, meta in metadata.items():
            all_metadata[f"{record.species}|{transcript_id}"] = meta
    return all_metadata


def classify_hits(
    hits: list[DomainHit],
    metadata: dict[str, dict],
    terminal_window: int,
) -> tuple[list[dict], list[dict]]:
    grouped: dict[tuple[str, str], list[DomainHit]] = defaultdict(list)
    for hit in hits:
        grouped[(hit.species, hit.target_id)].append(hit)

    rows: list[dict] = []
    audit_rows: list[dict] = []

    for (species, target_id), target_hits in sorted(grouped.items()):
        ntd_hits = [hit for hit in target_hits if hit.domain_type == "NTD"]
        ctd_hits = [hit for hit in target_hits if hit.domain_type == "CTD"]
        terminal_ntd = [hit for hit in ntd_hits if hit.is_terminal(terminal_window)]
        terminal_ctd = [hit for hit in ctd_hits if hit.is_terminal(terminal_window)]

        candidate_pairs = [
            (ntd, ctd)
            for ntd in terminal_ntd
            for ctd in terminal_ctd
            if ntd.spidroin_type == ctd.spidroin_type
        ]

        selected_ntd: Optional[DomainHit]
        selected_ctd: Optional[DomainHit]
        summed_score: Optional[float] = None
        if candidate_pairs:
            selected_ntd, selected_ctd = sorted(
                candidate_pairs,
                key=lambda pair: (
                    -(pair[0].dom_score + pair[1].dom_score),
                    max(pair[0].c_evalue, pair[1].c_evalue),
                ),
            )[0]
            classification = "full_length"
            spidroin_type = selected_ntd.spidroin_type
            summed_score = selected_ntd.dom_score + selected_ctd.dom_score
            for ntd, ctd in candidate_pairs:
                audit_rows.append(
                    {
                        "species": species,
                        "target_id": target_id,
                        "spidroin_type": ntd.spidroin_type,
                        "ntd_profile": ntd.query_id,
                        "ctd_profile": ctd.query_id,
                        "summed_dom_score": ntd.dom_score + ctd.dom_score,
                        "selected": ntd is selected_ntd and ctd is selected_ctd,
                    }
                )
        else:
            selected_ntd = best_hit(ntd_hits, terminal_window)
            selected_ctd = best_hit(ctd_hits, terminal_window)
            if terminal_ntd and not terminal_ctd:
                classification = "ntd_only"
                spidroin_type = selected_ntd.spidroin_type if selected_ntd else "Unknown"
            elif terminal_ctd and not terminal_ntd:
                classification = "ctd_only"
                spidroin_type = selected_ctd.spidroin_type if selected_ctd else "Unknown"
            else:
                classification = "both_hits_not_full"
                ntd_type = selected_ntd.spidroin_type if selected_ntd else ""
                ctd_type = selected_ctd.spidroin_type if selected_ctd else ""
                spidroin_type = ntd_type if ntd_type == ctd_type else "/".join(t for t in (ntd_type, ctd_type) if t)
                if not spidroin_type:
                    spidroin_type = "Unknown"

        _, transcript_id = split_target_id(target_id)
        meta = metadata.get(target_id, {})
        row = {
            "species": species,
            "target_id": target_id,
            "tlen": target_hits[0].tlen,
            "classification": classification,
            "spidroin_type": spidroin_type,
            "num_ntd_hits": len(ntd_hits),
            "num_ctd_hits": len(ctd_hits),
            "original_transcript_id": transcript_id,
            "gene_id": meta.get("gene_id", ""),
            "scaffold": meta.get("scaffold", ""),
            "strand": meta.get("strand", ""),
            "start": meta.get("start", ""),
            "end": meta.get("end", ""),
            "orfanage_status": meta.get("orfanage_status", ""),
            "orfanage_template": meta.get("orfanage_template", ""),
            "summed_dom_score": summed_score,
        }
        row.update(hit_columns("ntd", selected_ntd, terminal_window))
        row.update(hit_columns("ctd", selected_ctd, terminal_window))
        rows.append(row)

    return rows, audit_rows


RESULT_COLUMNS = [
    "species",
    "target_id",
    "tlen",
    "classification",
    "spidroin_type",
    "ntd_profile",
    "ctd_profile",
    "ntd_spidroin_type",
    "ctd_spidroin_type",
    "ntd_end_ok",
    "ctd_end_ok",
    "ntd_c_evalue",
    "ntd_dom_score",
    "ntd_env_from",
    "ntd_env_to",
    "ntd_hmm_from",
    "ntd_hmm_to",
    "ntd_hmm_coverage",
    "ntd_env_coverage",
    "ctd_c_evalue",
    "ctd_dom_score",
    "ctd_env_from",
    "ctd_env_to",
    "ctd_hmm_from",
    "ctd_hmm_to",
    "ctd_hmm_coverage",
    "ctd_env_coverage",
    "num_ntd_hits",
    "num_ctd_hits",
    "original_transcript_id",
    "gene_id",
    "scaffold",
    "strand",
    "start",
    "end",
    "orfanage_status",
    "orfanage_template",
    "summed_dom_score",
]


def write_pivot(full_rows: list[dict], species_names: list[str], output_path: Path) -> None:
    counts: dict[tuple[str, str], int] = Counter(
        (row["species"], row["spidroin_type"]) for row in full_rows
    )
    types = sorted({row["spidroin_type"] for row in full_rows})
    fieldnames = ["species", *types]
    rows = []
    for species in sorted(species_names):
        rows.append({"species": species, **{spid_type: counts[(species, spid_type)] for spid_type in types}})
    write_tsv(rows, output_path, fieldnames=fieldnames)


def load_proteins(protein_dir: Path) -> dict[str, SeqRecord]:
    proteins: dict[str, SeqRecord] = {}
    for fasta in sorted(protein_dir.glob("*.normalized.faa")):
        for record in SeqIO.parse(fasta, "fasta"):
            proteins[record.id] = record
    return proteins


def write_full_length_fastas(full_rows: list[dict], protein_dir: Path, processed_dir: Path) -> int:
    proteins = load_proteins(protein_dir)
    records: list[SeqRecord] = []
    by_type: dict[str, list[SeqRecord]] = defaultdict(list)
    missing: list[str] = []

    for row in full_rows:
        target_id = row["target_id"]
        source = proteins.get(target_id)
        if source is None:
            missing.append(target_id)
            continue
        record = SeqRecord(
            source.seq,
            id=target_id,
            description=(
                f"spidroin_type={row['spidroin_type']} "
                f"gene_id={row.get('gene_id', '')} "
                f"scaffold={row.get('scaffold', '')}:{row.get('start', '')}-{row.get('end', '')}"
            ),
        )
        records.append(record)
        by_type[row["spidroin_type"]].append(record)

    processed_dir.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, processed_dir / "full_length_spidroins.fasta", "fasta")

    by_type_dir = processed_dir / "by_type"
    by_type_dir.mkdir(parents=True, exist_ok=True)
    for spidroin_type, type_records in sorted(by_type.items()):
        safe_type = re.sub(r"[^A-Za-z0-9_.-]+", "_", spidroin_type)
        SeqIO.write(type_records, by_type_dir / f"{safe_type}.fasta", "fasta")

    if missing:
        write_tsv([{"target_id": target_id} for target_id in missing], processed_dir / "missing_fasta_records.tsv")
        logger.warning(f"Missing {len(missing)} full-length FASTA records")
    return len(records)


def validate_processed_outputs(processed_dir: Path) -> list[dict[str, str]]:
    full_tsv = processed_dir / "full_length_spidroins.tsv"
    full_fasta = processed_dir / "full_length_spidroins.fasta"
    checks: list[dict[str, str]] = []

    full_rows = []
    if full_tsv.exists():
        with full_tsv.open() as handle:
            full_rows = list(csv.DictReader(handle, delimiter="\t"))

    fasta_ids = [record.id for record in SeqIO.parse(full_fasta, "fasta")] if full_fasta.exists() else []
    duplicate_ids = len(fasta_ids) - len(set(fasta_ids))
    same_type_failures = [
        row["target_id"]
        for row in full_rows
        if row.get("ntd_spidroin_type") != row.get("ctd_spidroin_type")
    ]
    count_matches = len(full_rows) == len(fasta_ids)

    checks.extend(
        [
            {"check": "full_length_fasta_count_matches_tsv", "status": str(count_matches), "value": f"{len(fasta_ids)}/{len(full_rows)}"},
            {"check": "full_length_fasta_duplicate_ids", "status": str(duplicate_ids == 0), "value": str(duplicate_ids)},
            {"check": "full_length_same_type_ntd_ctd", "status": str(not same_type_failures), "value": str(len(same_type_failures))},
        ]
    )
    write_tsv(checks, processed_dir / "validation_summary.tsv")
    return checks


def run_classification(
    records: list[SpeciesRecord],
    orfanage_dir: Path,
    protein_dir: Path,
    hmmer_dir: Path,
    processed_dir: Path,
    c_evalue: float,
    min_hmm_coverage: float,
    terminal_window: int,
) -> None:
    metadata = load_metadata_for_species(records, orfanage_dir)
    all_hits: list[DomainHit] = []
    for record in records:
        all_hits.extend(
            parse_domtblout(
                hmmer_dir / f"{record.species}.domtblout",
                species=record.species,
                c_evalue_threshold=c_evalue,
                min_hmm_coverage=min_hmm_coverage,
            )
        )

    rows, audit_rows = classify_hits(all_hits, metadata, terminal_window)
    full_rows = [row for row in rows if row["classification"] == "full_length"]
    partial_rows = [row for row in rows if row["classification"] != "full_length"]

    processed_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(rows, processed_dir / "target_classification.tsv", RESULT_COLUMNS)
    write_tsv(full_rows, processed_dir / "full_length_spidroins.tsv", RESULT_COLUMNS)
    write_tsv(partial_rows, processed_dir / "partial_spidroins.tsv", RESULT_COLUMNS)
    write_tsv(
        audit_rows,
        processed_dir / "full_length_pair_audit.tsv",
        ["species", "target_id", "spidroin_type", "ntd_profile", "ctd_profile", "summed_dom_score", "selected"],
    )
    write_pivot(full_rows, [record.species for record in records], processed_dir / "spidroin_pivot_table.tsv")
    fasta_count = write_full_length_fastas(full_rows, protein_dir, processed_dir)
    validate_processed_outputs(processed_dir)
    logger.success(
        f"Classified {len(rows)} targets: {len(full_rows)} full-length, "
        f"{len(partial_rows)} partial/non-full; wrote {fasta_count} FASTA records."
    )


def write_qc_summary(records: list[SpeciesRecord], interim_dir: Path, processed_dir: Path) -> None:
    rows = []
    for record in records:
        orfanage_gtf = interim_dir / "orfanage" / f"{record.species}.orfanage.gtf"
        raw_faa = interim_dir / "proteins" / f"{record.species}.orfanage.proteins.faa"
        normalized_faa = interim_dir / "proteins" / f"{record.species}.orfanage.proteins.normalized.faa"
        rows.append(
            {
                "species": record.species,
                "orfanage_cds_features": count_feature(orfanage_gtf, "CDS"),
                "raw_proteins": sum(1 for _ in SeqIO.parse(raw_faa, "fasta")) if raw_faa.exists() else 0,
                "normalized_proteins": sum(1 for _ in SeqIO.parse(normalized_faa, "fasta")) if normalized_faa.exists() else 0,
            }
        )
    write_tsv(rows, processed_dir / "orfanage_protein_qc.tsv")


def preflight_checks(isoform_dir: Path, genome_dir: Path, hmm_dir: Path) -> list[dict[str, str]]:
    checks: list[dict[str, str]] = []
    for tool in ("orfanage", "gffread", "hmmsearch"):
        path = shutil.which(tool)
        checks.append({"check": f"executable:{tool}", "status": str(path is not None), "value": path or ""})
    checks.extend(
        [
            {"check": "isoform_dir_exists", "status": str(isoform_dir.exists()), "value": str(isoform_dir)},
            {"check": "genome_dir_exists", "status": str(genome_dir.exists()), "value": str(genome_dir)},
            {"check": "hmm_dir_exists", "status": str(hmm_dir.exists()), "value": str(hmm_dir)},
            {"check": "isoform_gtf_count", "status": str(len(list(isoform_dir.glob('*.isoforms.gtf'))) > 0), "value": str(len(list(isoform_dir.glob('*.isoforms.gtf'))))},
            {"check": "hmm_count", "status": str(len(list(hmm_dir.glob('*.hmm'))) > 0), "value": str(len(list(hmm_dir.glob('*.hmm'))))},
        ]
    )
    return checks


@app.command()
def preflight(
    isoform_dir: Path = typer.Option(DEFAULT_ISOFORM_DIR, help="ONT isoform GTF directory."),
    genome_dir: Path = typer.Option(DEFAULT_GENOME_DIR, help="Genome directory."),
    hmm_dir: Path = typer.Option(DEFAULT_HMM_DIR, help="HMM profile directory."),
    output: Optional[Path] = typer.Option(None, help="Optional TSV output path."),
) -> None:
    """Check required tools and input directories."""
    checks = preflight_checks(isoform_dir, genome_dir, hmm_dir)
    if output:
        write_tsv(checks, output, ["check", "status", "value"])
    for row in checks:
        logger.info(f"{row['check']}: {row['status']} {row['value']}")
    if not all(row["status"] == "True" for row in checks):
        raise typer.Exit(code=1)


@app.command()
def discover(
    isoform_dir: Path = typer.Option(DEFAULT_ISOFORM_DIR, help="ONT isoform GTF directory."),
    genome_dir: Path = typer.Option(DEFAULT_GENOME_DIR, help="Genome directory."),
    output: Path = typer.Option(DEFAULT_INTERIM_DIR / "species_manifest.tsv", help="Manifest TSV path."),
    species: Optional[list[str]] = typer.Option(None, "--species", help="Restrict to one species; repeatable."),
) -> None:
    """Discover matched isoform/genome/template annotation inputs."""
    _, rows = discover_species_records(isoform_dir, genome_dir, set(species or []) or None)
    write_tsv(rows, output, ["species", "isoform_gtf", "genome_fasta", "template_gff", "status", "reason"])
    matched = sum(1 for row in rows if row["status"] == "matched")
    logger.success(f"Wrote {output} with {matched}/{len(rows)} matched species")


@app.command()
def run(
    isoform_dir: Path = typer.Option(DEFAULT_ISOFORM_DIR, help="ONT isoform GTF directory."),
    genome_dir: Path = typer.Option(DEFAULT_GENOME_DIR, help="Genome directory."),
    hmm_dir: Path = typer.Option(DEFAULT_HMM_DIR, help="Protein HMM profile directory."),
    interim_dir: Path = typer.Option(DEFAULT_INTERIM_DIR, help="Intermediate output directory."),
    processed_dir: Path = typer.Option(DEFAULT_PROCESSED_DIR, help="Processed output directory."),
    species: Optional[list[str]] = typer.Option(None, "--species", help="Restrict to one species; repeatable."),
    threads: int = typer.Option(8, "--threads", min=1, help="Threads for ORFanage/HMMER."),
    evalue: float = typer.Option(1e-4, "--evalue", help="HMMER reporting and classification E-value."),
    min_hmm_coverage: float = typer.Option(0.90, "--min-hmm-coverage", help="Minimum HMM profile coverage."),
    terminal_window: int = typer.Option(300, "--terminal-window", help="Amino-acid distance allowed from protein termini."),
    force: bool = typer.Option(False, "--force", help="Regenerate existing outputs."),
) -> None:
    """Run the complete ONT ORFanage spidroin screening workflow."""
    checks = preflight_checks(isoform_dir, genome_dir, hmm_dir)
    failed = [row for row in checks if row["status"] != "True"]
    if failed:
        for row in failed:
            logger.error(f"Preflight failed: {row}")
        raise typer.Exit(code=1)

    records, manifest_rows = discover_species_records(isoform_dir, genome_dir, set(species or []) or None)
    if not records:
        raise typer.BadParameter("No matched species found.")

    interim_dir.mkdir(parents=True, exist_ok=True)
    processed_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(
        manifest_rows,
        interim_dir / "species_manifest.tsv",
        ["species", "isoform_gtf", "genome_fasta", "template_gff", "status", "reason"],
    )

    orfanage_dir = interim_dir / "orfanage"
    template_gff_dir = interim_dir / "template_gff_cleaned"
    protein_dir = interim_dir / "proteins"
    hmmer_dir = interim_dir / "hmmer"
    combined_hmm = hmmer_dir / "spidroin_terminals.all.hmm"
    concatenate_hmms(hmm_dir, combined_hmm, force=force)

    for record in tqdm(records, desc="Species"):
        cleaned_template_gff = template_gff_dir / f"{record.species}.cleaned.gff"
        repaired, swapped = clean_template_gff(record.template_gff, cleaned_template_gff, force=force)
        if repaired or swapped:
            logger.info(
                f"{record.species}: repaired {repaired} transcript rows and "
                f"swapped {swapped} other malformed feature rows in template GFF"
            )
        cleaned_record = SpeciesRecord(
            species=record.species,
            isoform_gtf=record.isoform_gtf,
            genome_fasta=record.genome_fasta,
            template_gff=cleaned_template_gff,
        )

        orfanage_gtf = orfanage_dir / f"{record.species}.orfanage.gtf"
        orfanage_stats = orfanage_dir / f"{record.species}.stats.tsv"
        raw_faa = protein_dir / f"{record.species}.orfanage.proteins.faa"
        normalized_faa = protein_dir / f"{record.species}.orfanage.proteins.normalized.faa"
        tblout = hmmer_dir / f"{record.species}.tblout"
        domtblout = hmmer_dir / f"{record.species}.domtblout"
        textout = hmmer_dir / f"{record.species}.hmmsearch.txt"

        run_orfanage_for_species(cleaned_record, orfanage_gtf, orfanage_stats, threads, force=force)
        protein_count = translate_for_species(
            record.genome_fasta,
            orfanage_gtf,
            raw_faa,
            normalized_faa,
            record.species,
            force=force,
        )
        logger.info(f"{record.species}: {protein_count} translated proteins")
        run_hmmer_for_species(
            combined_hmm,
            normalized_faa,
            tblout,
            domtblout,
            textout,
            threads=threads,
            evalue=evalue,
            force=force,
        )

    write_qc_summary(records, interim_dir, processed_dir)
    run_classification(
        records,
        orfanage_dir,
        protein_dir,
        hmmer_dir,
        processed_dir,
        c_evalue=evalue,
        min_hmm_coverage=min_hmm_coverage,
        terminal_window=terminal_window,
    )


if __name__ == "__main__":
    app()
