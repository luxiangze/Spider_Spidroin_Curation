from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from loguru import logger
from tqdm import tqdm
import typer

from spider_silkome_module.miniprot_window_screen.common import parse_attrs
from spider_silkome_module.parse_miniprot_window_gff import parse_window_id

app = typer.Typer(help="Export selected MPID models to genome-coordinate GFF3 and protein FASTA.")


EXPORT_FEATURES = {"mRNA", "CDS", "start_codon", "stop_codon"}
STOP_FIELDNAMES = ["terminal_stop_only", "internal_stop_count"]


def load_selected_rows(selected_tsv: Path) -> list[dict[str, str]]:
    if not selected_tsv.exists() or selected_tsv.stat().st_size == 0:
        return []
    with selected_tsv.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_models(models_jsonl: Path) -> dict[str, dict[str, Any]]:
    models: dict[str, dict[str, Any]] = {}
    with models_jsonl.open() as handle:
        for line in handle:
            if line.strip():
                model = json.loads(line)
                models[model["mrna_id"]] = model
    return models


def parse_gff_attr_text(attrs: dict[str, str]) -> str:
    return ";".join(f"{key}={value}" for key, value in attrs.items())


def format_float(value: Any) -> str:
    if value in (None, ""):
        return "NA"
    try:
        return f"{float(value):.4g}"
    except (TypeError, ValueError):
        return str(value)


def fasta_wrap(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def stop_codon_stats(seq: str) -> tuple[bool, int]:
    terminal_stop_only = seq.endswith("*") and seq.count("*") == 1
    seq_without_terminal = seq[:-1] if seq.endswith("*") else seq
    internal_stop_count = seq_without_terminal.count("*")
    return terminal_stop_only, internal_stop_count


def protein_sequence_for_export(seq: str) -> str | None:
    _terminal_stop_only, internal_stop_count = stop_codon_stats(seq)
    if internal_stop_count > 0:
        return None
    return seq[:-1] if seq.endswith("*") else seq


def selected_fieldnames(rows: list[dict[str, str]]) -> list[str]:
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames and key not in STOP_FIELDNAMES:
                fieldnames.append(key)
    fieldnames.extend(field for field in STOP_FIELDNAMES if field not in fieldnames)
    return fieldnames


def write_selected_rows(selected_tsv: Path, rows: list[dict[str, str]]) -> None:
    selected_tsv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = selected_fieldnames(rows)
    with selected_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def annotate_selected_stop_fields(
    selected_rows: list[dict[str, str]],
    models_by_id: dict[str, dict[str, Any]],
) -> list[dict[str, str]]:
    annotated_rows: list[dict[str, str]] = []
    for row in selected_rows:
        annotated = dict(row)
        model = models_by_id.get(row.get("mpid", ""))
        seq = (model or {}).get("protein_sequence") or ""
        terminal_stop_only, internal_stop_count = stop_codon_stats(seq)
        annotated["terminal_stop_only"] = str(terminal_stop_only)
        annotated["internal_stop_count"] = str(internal_stop_count)
        annotated_rows.append(annotated)
    return annotated_rows


def selected_metadata(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    return {row["mpid"]: row for row in rows if row.get("mpid")}


def write_selected_fasta(
    selected_rows: list[dict[str, str]],
    models_by_id: dict[str, dict[str, Any]],
    output_fasta: Path,
) -> int:
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with output_fasta.open("w") as handle:
        for row in selected_rows:
            mpid = row["mpid"]
            model = models_by_id.get(mpid)
            if not model:
                continue
            seq = model.get("protein_sequence") or ""
            if not seq:
                continue
            export_seq = protein_sequence_for_export(seq)
            if export_seq is None:
                continue
            header = (
                f">{row['species']}|{row['window_id']}|{mpid} "
                f"positive={format_float(row.get('positive'))} "
                f"query_coverage={format_float(row.get('query_coverage'))} "
                f"aa_len={row.get('protein_sequence_length', '')} "
                f"terminal_stop_only={row.get('terminal_stop_only', '')} "
                f"internal_stop_count={row.get('internal_stop_count', '')} "
                f"exon_count={row.get('exon_count', '')} "
                f"domain_pair_id={row.get('domain_pair_id', '')} "
                f"typing_spidroin_id={row.get('typing_spidroin_id', '')} "
                f"typing_spidroin_type={row.get('typing_spidroin_type', '')} "
                f"selection_mode={row.get('selection_mode', '')}"
            )
            handle.write(f"{header}\n{fasta_wrap(export_seq)}\n")
            n += 1
    return n


def should_export_feature(feature: str, attrs: dict[str, str], selected_ids: set[str]) -> bool:
    if feature == "mRNA":
        return attrs.get("ID") in selected_ids
    parent = attrs.get("Parent")
    return bool(parent and parent in selected_ids)


def write_selected_gff(
    selected_rows: list[dict[str, str]],
    miniprot_gff: Path,
    output_gff: Path,
) -> int:
    output_gff.parent.mkdir(parents=True, exist_ok=True)
    metadata = selected_metadata(selected_rows)
    selected_ids = set(metadata)
    n = 0
    with miniprot_gff.open() as inp, output_gff.open("w") as out:
        out.write("##gff-version 3\n")
        for line in tqdm(inp, desc=f"Export {output_gff.name}"):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            seqid, source, feature, start_s, end_s, score, strand, phase, attr_s = fields[:9]
            if feature not in EXPORT_FEATURES:
                continue
            attrs = parse_attrs(attr_s)
            if not should_export_feature(feature, attrs, selected_ids):
                continue
            mpid = attrs.get("ID") if feature == "mRNA" else attrs.get("Parent")
            if not mpid:
                continue
            genomic_seqid, window_start, _window_end = parse_window_id(seqid)
            start = window_start + int(start_s) - 1
            end = window_start + int(end_s) - 1
            meta = metadata[mpid]
            if feature == "mRNA":
                attrs["ID"] = mpid
                attrs["positive"] = str(meta.get("positive", ""))
                attrs["identity"] = str(meta.get("identity", ""))
                attrs["query_coverage"] = str(meta.get("query_coverage", ""))
                attrs["protein_sequence_length"] = str(meta.get("protein_sequence_length", ""))
                attrs["domain_pair_id"] = str(meta.get("domain_pair_id", ""))
                attrs["selection_status"] = str(meta.get("selection_status", ""))
                attrs["selection_mode"] = str(meta.get("selection_mode", ""))
                attrs["typing_spidroin_id"] = str(meta.get("typing_spidroin_id", ""))
                attrs["typing_spidroin_type"] = str(meta.get("typing_spidroin_type", ""))
            else:
                attrs["Parent"] = mpid
            out.write(
                "\t".join([
                    genomic_seqid,
                    source,
                    feature,
                    str(start),
                    str(end),
                    score,
                    strand,
                    phase,
                    parse_gff_attr_text(attrs),
                ])
                + "\n"
            )
            n += 1
    return n


@app.command()
def main(
    selected_tsv: Path,
    models_jsonl: Path,
    miniprot_gff: Path,
    output_gff: Path,
    output_fasta: Path,
) -> None:
    logger.info(f"Exporting selected MPIDs from {selected_tsv}")
    selected_rows = load_selected_rows(selected_tsv)
    models_by_id = load_models(models_jsonl)
    selected_rows = annotate_selected_stop_fields(selected_rows, models_by_id)
    write_selected_rows(selected_tsv, selected_rows)
    n_gff = write_selected_gff(selected_rows, miniprot_gff, output_gff)
    n_fasta = write_selected_fasta(selected_rows, models_by_id, output_fasta)
    logger.success(f"Exported {n_gff} GFF rows and {n_fasta} protein sequences")


if __name__ == "__main__":
    app()

