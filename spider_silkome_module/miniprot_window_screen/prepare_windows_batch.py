from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

from loguru import logger
from tqdm import tqdm
import typer

from spider_silkome_module.miniprot_window_screen.batch_common import (
    existing_path,
    fail_if_any_failed,
    read_tsv,
    summary_path,
    write_tsv,
)
from spider_silkome_module.miniprot_window_screen.prepare_windows import (
    bed_record_count,
    build_species_window_beds,
    extract_species_window_fasta,
    faidx_regions_path,
)

app = typer.Typer(help="Batch prepare candidate windows for miniprot window screening.")


def build_beds_for_row(row: dict[str, str], flank_bp: int, force: bool) -> dict[str, Any]:
    species = row["species"]
    raw_bed = Path(row["raw_bed"])
    merged_bed = Path(row["merged_bed"])
    try:
        if not force and raw_bed.exists() and merged_bed.exists():
            return {
                "species": species,
                "bed_status": "skipped",
                "bed_reason": "BED outputs exist",
                "n_raw_intervals": bed_record_count(raw_bed),
                "n_merged_windows": bed_record_count(merged_bed),
                "raw_bed": str(raw_bed),
                "merged_bed": str(merged_bed),
            }
        result = build_species_window_beds(
            species,
            Path(row["interim_dir"]),
            genome_fasta=Path(row["genome_fasta"]),
            nhmmer_gff=existing_path(row.get("nhmmer_gff")),
            miniprot_evidence_gff=Path(row["miniprot_evidence_gff"]),
            flank_bp=flank_bp,
            force=force,
        )
        return {
            "species": species,
            "bed_status": "ok",
            "bed_reason": "",
            "n_raw_intervals": result.get("n_raw_intervals", ""),
            "n_merged_windows": result.get("n_merged_windows", ""),
            "raw_bed": result.get("candidate_windows_raw_bed", str(raw_bed)),
            "merged_bed": result.get("candidate_windows_merged_bed", str(merged_bed)),
        }
    except Exception as exc:
        return {
            "species": species,
            "bed_status": "failed",
            "bed_reason": repr(exc),
            "n_raw_intervals": "",
            "n_merged_windows": "",
            "raw_bed": str(raw_bed),
            "merged_bed": str(merged_bed),
        }


def overall_status(bed_status: str, extract_status: str) -> str:
    if bed_status == "failed" or extract_status == "failed":
        return "failed"
    if bed_status == "skipped" and extract_status in {"skipped", "no_candidate_windows"}:
        return "skipped"
    return "ok"


@app.command()
def main(
    species_manifest: Path = typer.Option(..., "--species-manifest"),
    interim_root: Path = typer.Option(..., "--interim-root"),
    processed_root: Path = typer.Option(..., "--processed-root"),
    flank_bp: int = typer.Option(50_000, "--flank-bp", min=0),
    window_processes: int = typer.Option(70, "--window-processes", min=1),
    faidx_threads: int = typer.Option(70, "--faidx-threads", min=1),
    seqkit_threads: int | None = typer.Option(None, "--seqkit-threads", min=1),
    force: bool = False,
) -> None:
    if seqkit_threads is not None:
        logger.warning("--seqkit-threads is deprecated; using it as --faidx-threads")
        faidx_threads = seqkit_threads

    del processed_root
    rows = read_tsv(species_manifest)
    logger.info(
        f"Preparing windows for {len(rows)} species; "
        f"window_processes={window_processes}; faidx_threads={faidx_threads}"
    )

    bed_results: dict[str, dict[str, Any]] = {}
    if rows:
        max_workers = min(window_processes, len(rows))
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(build_beds_for_row, row, flank_bp, force): row["species"]
                for row in rows
            }
            for future in tqdm(as_completed(futures), total=len(futures), desc="Build window BEDs"):
                species = futures[future]
                bed_results[species] = future.result()

    summary_rows: list[dict[str, Any]] = []
    for row in tqdm(rows, desc="Extract window FASTA"):
        species = row["species"]
        bed_result = bed_results.get(species, {
            "species": species,
            "bed_status": "failed",
            "bed_reason": "missing BED build result",
            "n_raw_intervals": "",
            "n_merged_windows": "",
            "raw_bed": row.get("raw_bed", ""),
            "merged_bed": row.get("merged_bed", ""),
        })
        window_fasta = Path(row["window_fasta"])
        extract_status = "skipped"
        extract_reason = ""
        n_merged_windows = int(bed_result.get("n_merged_windows") or 0)

        if bed_result["bed_status"] == "failed":
            extract_status = "skipped"
            extract_reason = "BED build failed"
        else:
            try:
                if n_merged_windows == 0:
                    if force or not window_fasta.exists():
                        window_fasta.parent.mkdir(parents=True, exist_ok=True)
                        window_fasta.write_text("")
                    extract_status = "no_candidate_windows"
                    extract_reason = "no candidate windows"
                elif not force and window_fasta.exists():
                    extract_status = "skipped"
                    extract_reason = "FASTA output exists"
                else:
                    extract_species_window_fasta(
                        Path(row["genome_fasta"]),
                        Path(row["merged_bed"]),
                        window_fasta,
                        faidx_threads=faidx_threads,
                        force=force,
                    )
                    extract_status = "ok"
                    extract_reason = ""
            except Exception as exc:
                logger.exception(f"Window FASTA extraction failed for {species}")
                extract_status = "failed"
                extract_reason = repr(exc)

        reason = "; ".join(
            item for item in [str(bed_result.get("bed_reason") or ""), extract_reason] if item
        )
        summary_rows.append({
            "species": species,
            "status": overall_status(str(bed_result["bed_status"]), extract_status),
            "reason": reason,
            "bed_status": bed_result["bed_status"],
            "extract_status": extract_status,
            "raw_bed": bed_result.get("raw_bed", row.get("raw_bed", "")),
            "merged_bed": bed_result.get("merged_bed", row.get("merged_bed", "")),
            "window_fasta": str(window_fasta),
            "n_raw_intervals": bed_result.get("n_raw_intervals", ""),
            "n_merged_windows": bed_result.get("n_merged_windows", ""),
            "window_processes": window_processes,
            "extract_tool": "samtools faidx",
            "faidx_threads": faidx_threads,
            "faidx_regions": str(faidx_regions_path(window_fasta)),
        })

    out = summary_path(interim_root, "prepare_windows")
    write_tsv(out, summary_rows)
    logger.success(f"Window preparation summary: {out}")
    fail_if_any_failed(summary_rows)


if __name__ == "__main__":
    app()
