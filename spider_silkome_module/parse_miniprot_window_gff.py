from __future__ import annotations

import csv
from dataclasses import asdict, dataclass, field
import json
from pathlib import Path
from typing import Any

from Bio import SeqIO
from Bio.Seq import Seq
import typer

app = typer.Typer()


@dataclass
class MiniprotCdsBlock:
    parent_id: str
    window_id: str
    local_start: int
    local_end: int
    genomic_seqid: str
    genomic_start: int
    genomic_end: int
    score: float | None
    strand: str
    phase: str
    rank: int | None
    identity: float | None
    target_protein: str | None
    target_start: int | None
    target_end: int | None
    nucleotide_sequence: str
    protein_sequence: str


@dataclass
class MiniprotWindowGene:
    mrna_id: str
    window_id: str
    genomic_seqid: str
    window_start: int
    window_end: int
    local_start: int
    local_end: int
    genomic_start: int
    genomic_end: int
    score: float | None
    strand: str
    rank: int | None
    identity: float | None
    positive: float | None
    target_protein: str | None
    target_start: int | None
    target_end: int | None
    query_length: int | None
    query_coverage: float | None
    paf_query_coverage: float | None
    paf: dict[str, Any] | None
    exon_count: int
    intron_count: int
    mrna_sequence: str
    spliced_cds_sequence: str
    protein_sequence: str
    cds_blocks: list[MiniprotCdsBlock] = field(default_factory=list)

    @property
    def protein_sequence_length(self) -> int:
        return len(self.protein_sequence)


def parse_attrs(attr_str: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in attr_str.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value
    return attrs


def parse_float(value: str | None) -> float | None:
    if value in (None, "", "."):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def parse_int(value: str | None) -> int | None:
    if value in (None, "", "."):
        return None
    try:
        return int(value)
    except ValueError:
        return None


def parse_target(value: str | None) -> tuple[str | None, int | None, int | None]:
    if not value:
        return None, None, None
    parts = value.split()
    target = parts[0] if parts else None
    start = parse_int(parts[1]) if len(parts) > 1 else None
    end = parse_int(parts[2]) if len(parts) > 2 else None
    return target, start, end


def parse_window_id(window_id: str) -> tuple[str, int, int]:
    try:
        seqid, coord_text = window_id.rsplit("_", 1)
        coord_text = coord_text.rstrip(".").rstrip(":")
        start_text, end_text = coord_text.split("-", 1)
        return seqid, int(start_text), int(end_text)
    except ValueError as exc:
        raise ValueError(f"Cannot parse window id: {window_id}") from exc


def load_fasta_sequences(fasta: Path) -> dict[str, str]:
    return {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(fasta, "fasta")}


def load_fasta_lengths(fasta: Path) -> dict[str, int]:
    return {rec.id: len(rec.seq) for rec in SeqIO.parse(fasta, "fasta")}


def parse_paf_comment(line: str) -> dict[str, Any] | None:
    fields = line.rstrip("\n").split("\t")
    if not fields or fields[0] != "##PAF" or len(fields) < 13:
        return None

    tags: dict[str, str] = {}
    for paf_tag in fields[13:]:
        key = paf_tag.split(":", 1)[0]
        tags[key] = paf_tag

    return {
        "query_name": fields[1],
        "query_length": parse_int(fields[2]),
        "query_start": parse_int(fields[3]),
        "query_end": parse_int(fields[4]),
        "strand": fields[5],
        "target_name": fields[6],
        "target_length": parse_int(fields[7]),
        "target_start": parse_int(fields[8]),
        "target_end": parse_int(fields[9]),
        "residue_matches": parse_int(fields[10]),
        "alignment_block_length": parse_int(fields[11]),
        "mapping_quality": parse_int(fields[12]),
        "tags": tags,
    }


def translate_trimmed(nt_sequence: str) -> str:
    usable_len = len(nt_sequence) - (len(nt_sequence) % 3)
    if usable_len < 3:
        return ""
    return str(Seq(nt_sequence[:usable_len]).translate(to_stop=False))


def extract_interval(seq: str, start: int, end: int, strand: str) -> str:
    subseq = seq[start - 1:end]
    if strand == "-":
        return str(Seq(subseq).reverse_complement())
    return subseq


def transcript_order(blocks: list[dict[str, Any]], strand: str) -> list[dict[str, Any]]:
    reverse = strand == "-"
    return sorted(blocks, key=lambda item: item["local_start"], reverse=reverse)


def build_models(
    miniprot_gff: Path,
    window_fasta: Path,
    query_fasta: Path,
) -> list[MiniprotWindowGene]:
    window_sequences = load_fasta_sequences(window_fasta)
    query_lengths = load_fasta_lengths(query_fasta)

    mrna_records: dict[str, dict[str, Any]] = {}
    cds_records: dict[str, list[dict[str, Any]]] = {}
    pending_paf: dict[str, Any] | None = None

    with miniprot_gff.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("##PAF\t"):
                pending_paf = parse_paf_comment(line)
                continue
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            seqid, source, feature, start, end, score, strand, phase, attr = fields[:9]
            attrs = parse_attrs(attr)
            target, target_start, target_end = parse_target(attrs.get("Target"))
            row = {
                "window_id": seqid,
                "source": source,
                "feature": feature,
                "local_start": int(start),
                "local_end": int(end),
                "score": parse_float(score),
                "strand": strand,
                "phase": phase,
                "attrs": attrs,
                "rank": parse_int(attrs.get("Rank")),
                "identity": parse_float(attrs.get("Identity")),
                "positive": parse_float(attrs.get("Positive")),
                "target_protein": target,
                "target_start": target_start,
                "target_end": target_end,
            }

            if feature == "mRNA":
                mrna_id = attrs.get("ID")
                if not mrna_id:
                    continue
                row["paf"] = pending_paf
                mrna_records[mrna_id] = row
                pending_paf = None
            elif feature == "CDS":
                parent = attrs.get("Parent")
                if parent:
                    cds_records.setdefault(parent, []).append(row)

    models: list[MiniprotWindowGene] = []
    for mrna_id, mrna in mrna_records.items():
        window_id = mrna["window_id"]
        if window_id not in window_sequences:
            raise KeyError(f"{window_id} was not found in {window_fasta}")
        genomic_seqid, window_start, window_end = parse_window_id(window_id)
        window_seq = window_sequences[window_id]

        genomic_start = window_start + mrna["local_start"] - 1
        genomic_end = window_start + mrna["local_end"] - 1
        target_protein = mrna["target_protein"]
        query_length = query_lengths.get(target_protein) if target_protein else None
        if query_length is None and mrna.get("paf"):
            query_length = mrna["paf"].get("query_length")

        query_coverage = None
        if query_length and mrna["target_start"] is not None and mrna["target_end"] is not None:
            query_coverage = (abs(mrna["target_end"] - mrna["target_start"]) + 1) / query_length

        paf_query_coverage = None
        paf = mrna.get("paf")
        if paf and paf.get("query_length") and paf.get("query_start") is not None and paf.get("query_end") is not None:
            paf_query_coverage = (paf["query_end"] - paf["query_start"]) / paf["query_length"]

        raw_cds_blocks = transcript_order(cds_records.get(mrna_id, []), mrna["strand"])
        cds_blocks: list[MiniprotCdsBlock] = []
        spliced_parts: list[str] = []
        for cds in raw_cds_blocks:
            nt_seq = extract_interval(window_seq, cds["local_start"], cds["local_end"], mrna["strand"])
            spliced_parts.append(nt_seq)
            block = MiniprotCdsBlock(
                parent_id=mrna_id,
                window_id=window_id,
                local_start=cds["local_start"],
                local_end=cds["local_end"],
                genomic_seqid=genomic_seqid,
                genomic_start=window_start + cds["local_start"] - 1,
                genomic_end=window_start + cds["local_end"] - 1,
                score=cds["score"],
                strand=cds["strand"],
                phase=cds["phase"],
                rank=cds["rank"],
                identity=cds["identity"],
                target_protein=cds["target_protein"],
                target_start=cds["target_start"],
                target_end=cds["target_end"],
                nucleotide_sequence=nt_seq,
                protein_sequence=translate_trimmed(nt_seq),
            )
            cds_blocks.append(block)

        spliced_cds = "".join(spliced_parts)
        model = MiniprotWindowGene(
            mrna_id=mrna_id,
            window_id=window_id,
            genomic_seqid=genomic_seqid,
            window_start=window_start,
            window_end=window_end,
            local_start=mrna["local_start"],
            local_end=mrna["local_end"],
            genomic_start=genomic_start,
            genomic_end=genomic_end,
            score=mrna["score"],
            strand=mrna["strand"],
            rank=mrna["rank"],
            identity=mrna["identity"],
            positive=mrna["positive"],
            target_protein=target_protein,
            target_start=mrna["target_start"],
            target_end=mrna["target_end"],
            query_length=query_length,
            query_coverage=query_coverage,
            paf_query_coverage=paf_query_coverage,
            paf=paf,
            exon_count=len(cds_blocks),
            intron_count=max(0, len(cds_blocks) - 1),
            mrna_sequence=extract_interval(window_seq, mrna["local_start"], mrna["local_end"], mrna["strand"]),
            spliced_cds_sequence=spliced_cds,
            protein_sequence=translate_trimmed(spliced_cds),
            cds_blocks=cds_blocks,
        )
        models.append(model)

    return models


def summary_row(model: MiniprotWindowGene) -> dict[str, Any]:
    return {
        "mrna_id": model.mrna_id,
        "window_id": model.window_id,
        "genomic_seqid": model.genomic_seqid,
        "window_start": model.window_start,
        "window_end": model.window_end,
        "genomic_start": model.genomic_start,
        "genomic_end": model.genomic_end,
        "strand": model.strand,
        "score": model.score,
        "rank": model.rank,
        "identity": model.identity,
        "positive": model.positive,
        "target_protein": model.target_protein,
        "target_start": model.target_start,
        "target_end": model.target_end,
        "query_length": model.query_length,
        "query_coverage": model.query_coverage,
        "paf_query_coverage": model.paf_query_coverage,
        "exon_count": model.exon_count,
        "intron_count": model.intron_count,
        "mrna_sequence_length": len(model.mrna_sequence),
        "spliced_cds_sequence_length": len(model.spliced_cds_sequence),
        "protein_sequence_length": model.protein_sequence_length,
    }


def write_jsonl(models: list[MiniprotWindowGene], output_jsonl: Path) -> None:
    output_jsonl.parent.mkdir(parents=True, exist_ok=True)
    with output_jsonl.open("w") as out:
        for model in models:
            out.write(json.dumps(asdict(model), ensure_ascii=False) + "\n")


def write_summary(models: list[MiniprotWindowGene], output_summary: Path) -> None:
    output_summary.parent.mkdir(parents=True, exist_ok=True)
    rows = [summary_row(model) for model in models]
    if not rows:
        output_summary.write_text("")
        return

    with output_summary.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


@app.command()
def main(
    miniprot_gff: Path,
    window_fasta: Path,
    query_fasta: Path,
    output_jsonl: Path,
    output_summary: Path,
) -> None:
    models = build_models(miniprot_gff, window_fasta, query_fasta)
    models.sort(key=lambda m: (m.window_id, m.genomic_start, m.genomic_end, m.mrna_id))
    write_jsonl(models, output_jsonl)
    write_summary(models, output_summary)
    typer.echo(f"Parsed {len(models)} miniprot mRNA models")
    typer.echo(f"JSONL: {output_jsonl}")
    typer.echo(f"Summary: {output_summary}")


if __name__ == "__main__":
    app()
