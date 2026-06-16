from __future__ import annotations

from pathlib import Path
from typing import Any

from loguru import logger
from tqdm import tqdm
import typer

from spider_silkome_module.miniprot_window_screen.batch_common import (
    fail_if_any_failed,
    read_tsv,
    summary_path,
    write_tsv,
)
from spider_silkome_module.miniprot_window_screen.split_miniprot_output import (
    split_miniprot_mixed_output,
)

app = typer.Typer(help="Batch split miniprot mixed GFF/PAF output.")


@app.command()
def main(
    species_manifest: Path = typer.Option(..., "--species-manifest"),
    interim_root: Path = typer.Option(..., "--interim-root"),
    processed_root: Path = typer.Option(..., "--processed-root"),
    force: bool = False,
) -> None:
    del processed_root
    rows = read_tsv(species_manifest)
    summary_rows: list[dict[str, Any]] = []
    logger.info(f"Splitting miniprot output for {len(rows)} species")

    for row in tqdm(rows, desc="Split miniprot"):
        species = row["species"]
        mixed = Path(row["miniprot_mixed"])
        gff = Path(row["miniprot_gff"])
        paf = Path(row["miniprot_paf"])
        unclassified = Path(row["miniprot_unclassified"])
        try:
            if not mixed.exists():
                status = "skipped"
                reason = "missing miniprot output"
                counts = {"gff": 0, "paf": 0, "unclassified": 0}
            elif not force and gff.exists() and paf.exists() and unclassified.exists():
                status = "skipped"
                reason = "outputs exist"
                counts = {"gff": "", "paf": "", "unclassified": ""}
            else:
                counts = split_miniprot_mixed_output(mixed, gff, paf, unclassified)
                status = "ok"
                reason = ""
            summary_rows.append({
                "species": species,
                "status": status,
                "reason": reason,
                "input": str(mixed),
                "gff": str(gff),
                "paf": str(paf),
                "unclassified": str(unclassified),
                "n_gff": counts["gff"],
                "n_paf": counts["paf"],
                "n_unclassified": counts["unclassified"],
            })
        except Exception as exc:
            logger.exception(f"Split failed for {species}")
            summary_rows.append({
                "species": species,
                "status": "failed",
                "reason": repr(exc),
                "input": str(mixed),
                "gff": str(gff),
                "paf": str(paf),
                "unclassified": str(unclassified),
                "n_gff": "",
                "n_paf": "",
                "n_unclassified": "",
            })

    out = summary_path(interim_root, "split_miniprot")
    write_tsv(out, summary_rows)
    logger.success(f"Split summary: {out}")
    fail_if_any_failed(summary_rows)


if __name__ == "__main__":
    app()
