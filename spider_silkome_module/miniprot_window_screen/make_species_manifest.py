from __future__ import annotations

import csv
from pathlib import Path
from typing import Annotated

from loguru import logger
from tqdm import tqdm
import typer

from spider_silkome_module.config import PROCESSED_DATA_DIR, RAW_DATA_DIR
from spider_silkome_module.miniprot_window_screen.batch_common import write_tsv
from spider_silkome_module.miniprot_window_screen.common import discover_species, find_genome_fasta

app = typer.Typer(help="Create a species manifest for miniprot window screening.")


def parse_bool(value: str | None) -> bool:
    return str(value or "").lower() == "true"


def typing_stats(typing_tsv: Path, typing_gff: Path) -> tuple[str, int, int, str]:
    if typing_tsv.exists() and typing_tsv.stat().st_size > 0:
        with typing_tsv.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        n_full = sum(parse_bool(row.get("Full_length")) for row in rows)
        return "present", len(rows), n_full, "typing_guided" if rows else "fallback"

    if typing_gff.exists() and typing_gff.stat().st_size > 0:
        n_loci = 0
        n_full = 0
        with typing_gff.open() as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9 or fields[2] != "spidroin_gene":
                    continue
                n_loci += 1
                if "hint_type=Full_length" in fields[8]:
                    n_full += 1
        return "present", n_loci, n_full, "typing_guided" if n_loci else "fallback"

    return "missing", 0, 0, "fallback"


def species_row(
    task_name: str,
    species: str,
    interim_root: Path,
    processed_root: Path,
    raw_ref_dir: Path,
    miniprot_evidence_dir: Path,
    nhmmer_parsed_dir: Path,
    typing_results_dir: Path,
) -> dict[str, str | int]:
    interim = interim_root / species
    processed = processed_root / species
    pygt = interim / "pygenometracks"
    typing_dir = typing_results_dir / species
    nhmmer_gff = nhmmer_parsed_dir / species / f"{species}.gff"
    typing_tsv = typing_dir / f"{species}.tsv"
    typing_gff = typing_dir / f"{species}.gff"
    miniprot_evidence_gff = miniprot_evidence_dir / f"{species}.gff"
    typing_status, typing_locus_count, typing_full_length_count, selection_mode = typing_stats(typing_tsv, typing_gff)

    return {
        "task_name": task_name,
        "species": species,
        "genome_fasta": str(find_genome_fasta(species, raw_ref_dir)),
        "nhmmer_gff": str(nhmmer_gff),
        "nhmmer_status": "present" if nhmmer_gff.exists() else "missing",
        "typing_tsv": str(typing_tsv),
        "typing_gff": str(typing_gff),
        "typing_status": typing_status,
        "typing_locus_count": typing_locus_count,
        "typing_full_length_count": typing_full_length_count,
        "selection_mode": selection_mode,
        "miniprot_evidence_gff": str(miniprot_evidence_gff),
        "interim_dir": str(interim),
        "processed_dir": str(processed),
        "raw_bed": str(interim / "candidate_windows.raw.bed"),
        "merged_bed": str(interim / "candidate_windows.merged.bed"),
        "window_fasta": str(interim / "candidate_windows.fa"),
        "miniprot_mixed": str(interim / "miniprot_vs_windows.gff"),
        "miniprot_gff": str(interim / "miniprot_vs_windows.only.gff3"),
        "miniprot_paf": str(interim / "miniprot_vs_windows.paf"),
        "miniprot_unclassified": str(interim / "miniprot_vs_windows.unclassified.txt"),
        "models_jsonl": str(interim / "miniprot_window_mrna_models.jsonl"),
        "models_summary": str(interim / "miniprot_window_mrna_summary.tsv"),
        "pygt_commands": str(pygt / "pygenometracks_commands.tsv"),
        "selected_tsv": str(processed / "selected_spidroin_mpid.tsv"),
        "selected_gff": str(processed / "selected_spidroin_models.gff3"),
        "selected_faa": str(processed / "selected_spidroin_proteins.faa"),
    }


@app.command()
def main(
    species_manifest: Path = typer.Option(..., "--species-manifest"),
    interim_root: Path = typer.Option(..., "--interim-root"),
    processed_root: Path = typer.Option(..., "--processed-root"),
    task_name: str = typer.Option(..., "--task-name"),
    raw_ref_dir: Path = RAW_DATA_DIR / "01.ref_gff",
    miniprot_evidence_dir: Path = PROCESSED_DATA_DIR / "miniprot_output",
    nhmmer_parsed_dir: Path = PROCESSED_DATA_DIR / "nhmmer_search_parsed",
    typing_results_dir: Path = PROCESSED_DATA_DIR / "typing_results",
    species_filter: Annotated[list[str] | None, typer.Option("--species-filter")] = None,
    max_species: int | None = typer.Option(None, "--max-species"),
    force: bool = False,
) -> None:
    del force
    logger.info(f"Creating species manifest: {species_manifest}")
    species_all = discover_species(raw_ref_dir, miniprot_evidence_dir)
    if species_filter:
        wanted = set(species_filter)
        species_list = [species for species in species_all if species in wanted]
    else:
        species_list = species_all
    if max_species is not None:
        species_list = species_list[:max_species]

    rows = [
        species_row(
            task_name,
            species,
            interim_root,
            processed_root,
            raw_ref_dir,
            miniprot_evidence_dir,
            nhmmer_parsed_dir,
            typing_results_dir,
        )
        for species in tqdm(species_list, desc="Manifest species")
    ]
    write_tsv(species_manifest, rows)
    logger.success(f"Manifest species={len(rows)} available={len(species_all)} selected={len(species_list)}")


if __name__ == "__main__":
    app()
