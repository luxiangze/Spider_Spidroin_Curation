from __future__ import annotations

import csv
from pathlib import Path
import subprocess
from typing import Any


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = []
        for row in rows:
            for key in row:
                if key not in fieldnames:
                    fieldnames.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def path_or_none(value: str | None) -> Path | None:
    if not value:
        return None
    return Path(value)


def existing_path(value: str | None) -> Path | None:
    path = path_or_none(value)
    return path if path and path.exists() else None


def exists_and_nonempty(path: Path | None) -> bool:
    return bool(path and path.exists() and path.stat().st_size > 0)


def outputs_ready(paths: list[Path], force: bool = False) -> bool:
    return (not force) and all(path.exists() for path in paths)


def count_fasta_records(path: Path | None) -> int:
    if not exists_and_nonempty(path):
        return 0
    assert path is not None
    with path.open() as handle:
        return sum(1 for line in handle if line.startswith(">"))


def count_tsv_records(path: Path | None) -> int:
    if not exists_and_nonempty(path):
        return 0
    assert path is not None
    with path.open() as handle:
        return max(0, sum(1 for _ in handle) - 1)


def count_noncomment_records(path: Path | None) -> int:
    if not exists_and_nonempty(path):
        return 0
    assert path is not None
    with path.open() as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith("#"))


def manifests_dir(interim_root: Path) -> Path:
    path = interim_root / "manifests"
    path.mkdir(parents=True, exist_ok=True)
    return path


def summary_path(interim_root: Path, step_name: str) -> Path:
    return manifests_dir(interim_root) / f"{step_name}_summary.tsv"


def run_shell(command: str) -> None:
    subprocess.run(command, shell=True, check=True)


def fail_if_any_failed(rows: list[dict[str, Any]]) -> None:
    if any(row.get("status") == "failed" for row in rows):
        raise SystemExit(1)
