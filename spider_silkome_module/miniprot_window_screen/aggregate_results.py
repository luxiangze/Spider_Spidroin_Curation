from __future__ import annotations

from pathlib import Path
from typing import Any

from loguru import logger
from tqdm import tqdm
import typer

from spider_silkome_module.miniprot_window_screen.batch_common import (
    count_fasta_records,
    count_tsv_records,
    read_tsv,
    write_tsv,
)

app = typer.Typer(help="Aggregate selected MPID TSV/FASTA outputs across species.")

SELECTED_FIELD_ORDER = [
    "species",
    "window_id",
    "mpid",
    "selection_status",
    "selection_mode",
    "typing_spidroin_id",
    "typing_spidroin_type",
    "typing_hint_type",
    "typing_full_length",
    "typing_confidence",
    "typing_needs_review",
    "typing_overlap_frac",
    "domain_pair_id",
    "domain_pair_start",
    "domain_pair_end",
    "domain_support_count",
    "domain_support_hits",
    "miniprot_evidence_overlap_count",
    "genomic_seqid",
    "genomic_start",
    "genomic_end",
    "strand",
    "positive",
    "identity",
    "query_coverage",
    "score",
    "protein_sequence_length",
    "terminal_stop_only",
    "internal_stop_count",
    "exon_count",
    "intron_count",
    "target_protein",
    "selection_score",
]

MANIFEST_FIELD_ORDER = [
    "task_name",
    "species",
    "status",
    "reason",
    "nhmmer_status",
    "typing_status",
    "typing_locus_count",
    "typing_full_length_count",
    "selection_mode",
    "n_selected",
    "n_selected_proteins",
    "n_plots",
    "selected_tsv",
    "selected_gff",
    "selected_faa",
]


def read_selected_rows(path: Path) -> list[dict[str, str]]:
    return read_tsv(path)


@app.command()
def main(
    species_manifest: Path = typer.Option(..., "--species-manifest"),
    interim_root: Path = typer.Option(..., "--interim-root"),
    processed_root: Path = typer.Option(..., "--processed-root"),
    force: bool = False,
) -> None:
    del interim_root, force
    manifest_rows = read_tsv(species_manifest)
    all_selected_rows: list[dict[str, Any]] = []
    final_manifest: list[dict[str, Any]] = []

    all_selected_tsv = processed_root / "all_species_selected_spidroin_mpid.tsv"
    all_selected_faa = processed_root / "all_species_selected_spidroin_proteins.faa"
    species_manifest_out = processed_root / "species_run_manifest.tsv"
    processed_root.mkdir(parents=True, exist_ok=True)

    with all_selected_faa.open("w") as faa_out:
        for row in tqdm(manifest_rows, desc="Aggregate species"):
            selected_tsv = Path(row["selected_tsv"])
            selected_faa = Path(row["selected_faa"])
            selected_gff = Path(row["selected_gff"])
            selected_rows = read_selected_rows(selected_tsv) if selected_tsv.exists() else []
            all_selected_rows.extend(selected_rows)
            if selected_faa.exists() and selected_faa.stat().st_size > 0:
                text = selected_faa.read_text()
                faa_out.write(text)
                if not text.endswith("\n"):
                    faa_out.write("\n")
            if selected_rows:
                selection_modes = sorted({item.get("selection_mode", "") for item in selected_rows if item.get("selection_mode")})
                status = "ok"
                reason = ""
                selection_mode = ",".join(selection_modes) or row.get("selection_mode", "")
            elif selected_tsv.exists():
                status = "ok"
                reason = "no selected MPIDs"
                selection_mode = row.get("selection_mode", "")
            else:
                status = "missing_selected"
                reason = "missing selected TSV"
                selection_mode = row.get("selection_mode", "")
            final_manifest.append({
                "task_name": row.get("task_name", ""),
                "species": row["species"],
                "status": status,
                "reason": reason,
                "nhmmer_status": row.get("nhmmer_status", ""),
                "typing_status": row.get("typing_status", ""),
                "typing_locus_count": row.get("typing_locus_count", ""),
                "typing_full_length_count": row.get("typing_full_length_count", ""),
                "selection_mode": selection_mode,
                "n_selected": count_tsv_records(selected_tsv),
                "n_selected_proteins": count_fasta_records(selected_faa),
                "n_plots": len(list((Path(row["interim_dir"]) / "pygenometracks" / "plots").glob("*.png"))),
                "selected_tsv": str(selected_tsv),
                "selected_gff": str(selected_gff),
                "selected_faa": str(selected_faa),
            })

    write_tsv(all_selected_tsv, all_selected_rows, SELECTED_FIELD_ORDER)
    write_tsv(species_manifest_out, final_manifest, MANIFEST_FIELD_ORDER)
    logger.success(
        f"Aggregated selected_rows={len(all_selected_rows)} proteins={count_fasta_records(all_selected_faa)} "
        f"manifest={species_manifest_out}"
    )


if __name__ == "__main__":
    app()
