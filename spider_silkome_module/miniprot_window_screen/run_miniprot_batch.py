from __future__ import annotations

from pathlib import Path
import subprocess
from typing import Any

from loguru import logger
from tqdm import tqdm
import typer

from spider_silkome_module.miniprot_window_screen.batch_common import (
    count_fasta_records,
    count_noncomment_records,
    fail_if_any_failed,
    read_tsv,
    summary_path,
    write_tsv,
)

app = typer.Typer(help="Batch run miniprot --gff on candidate window FASTA files.")


@app.command()
def main(
    species_manifest: Path = typer.Option(..., "--species-manifest"),
    interim_root: Path = typer.Option(..., "--interim-root"),
    processed_root: Path = typer.Option(..., "--processed-root"),
    query_proteins: Path = typer.Option(..., "--query-proteins"),
    threads: int = typer.Option(16, "--threads", min=1),
    force: bool = False,
) -> None:
    del processed_root
    rows = read_tsv(species_manifest)
    summary_rows: list[dict[str, Any]] = []
    logger.info(f"Running miniprot for {len(rows)} species")

    for row in tqdm(rows, desc="Run miniprot"):
        species = row["species"]
        window_fasta = Path(row["window_fasta"])
        output = Path(row["miniprot_mixed"])
        try:
            n_windows = count_fasta_records(window_fasta)
            if n_windows == 0:
                status = "skipped"
                reason = "no candidate windows"
            elif not force and output.exists():
                status = "skipped"
                reason = "output exists"
            else:
                output.parent.mkdir(parents=True, exist_ok=True)
                cmd = [
                    "miniprot",
                    "--gff",
                    "-t",
                    str(threads),
                    "-G",
                    "200000",
                    "--outc",
                    "0.05",
                    "--outs",
                    "0.10",
                    "--outn",
                    "1000",
                    str(window_fasta),
                    str(query_proteins),
                ]
                with output.open("w") as handle:
                    subprocess.run(cmd, check=True, stdout=handle)
                status = "ok"
                reason = ""
            summary_rows.append({
                "species": species,
                "status": status,
                "reason": reason,
                "input": str(window_fasta),
                "output": str(output),
                "n_windows": n_windows,
                "n_output_records": count_noncomment_records(output),
            })
        except Exception as exc:
            logger.exception(f"Miniprot failed for {species}")
            summary_rows.append({
                "species": species,
                "status": "failed",
                "reason": repr(exc),
                "input": str(window_fasta),
                "output": str(output),
                "n_windows": "",
                "n_output_records": "",
            })

    out = summary_path(interim_root, "run_miniprot")
    write_tsv(out, summary_rows)
    logger.success(f"Miniprot summary: {out}")
    fail_if_any_failed(summary_rows)


if __name__ == "__main__":
    app()
