from __future__ import annotations

import csv
from pathlib import Path
import shlex

from loguru import logger
from tqdm import tqdm
import typer

from spider_silkome_module.config import PROCESSED_DATA_DIR
from spider_silkome_module.miniprot_window_screen.common import (
    Interval,
    find_genome_fasta,
    parse_attrs,
    read_fai,
)
from spider_silkome_module.utils.run_cmd import run_cmd

app = typer.Typer(help="Prepare merged candidate genome windows for miniprot spidroin screening.")


def shell_path(path: Path) -> str:
    return shlex.quote(str(path))


def parse_evidence_intervals(
    gff_path: Path,
    source_label: str,
    features: set[str],
    contig_lengths: dict[str, int],
    flank_bp: int,
) -> list[Interval]:
    intervals: list[Interval] = []
    if not gff_path.exists():
        logger.warning(f"Missing {source_label} GFF: {gff_path}")
        return intervals

    with gff_path.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _src, feature, start_s, end_s, score, strand, _phase, attrs_s = fields[:9]
            if feature not in features or chrom not in contig_lengths:
                continue
            start = int(start_s)
            end = int(end_s)
            start0 = max(0, start - 1 - flank_bp)
            end0 = min(contig_lengths[chrom], end + flank_bp)
            attrs = parse_attrs(attrs_s)
            name = attrs.get("Name") or attrs.get("ID") or attrs.get("Target", ".").split()[0] or feature
            label = f"{source_label}:{feature}:{name}:{strand}:{score}"
            if end0 > start0:
                intervals.append(Interval(chrom, start0, end0, label))
    return intervals


def merge_intervals(intervals: list[Interval]) -> list[tuple[str, int, int, str, int]]:
    if not intervals:
        return []

    merged: list[tuple[str, int, int, set[str], int]] = []
    for interval in sorted(intervals, key=lambda item: (item.chrom, item.start0, item.end0, item.label)):
        if not merged or interval.chrom != merged[-1][0] or interval.start0 > merged[-1][2]:
            merged.append((interval.chrom, interval.start0, interval.end0, {interval.label}, 1))
            continue
        chrom, start0, end0, labels, count = merged[-1]
        labels.add(interval.label)
        merged[-1] = (chrom, start0, max(end0, interval.end0), labels, count + 1)

    return [
        (chrom, start0, end0, ",".join(sorted(labels)), count)
        for chrom, start0, end0, labels, count in merged
    ]


def write_raw_bed(path: Path, intervals: list[Interval]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for interval in intervals:
            writer.writerow([interval.chrom, interval.start0, interval.end0, interval.label])


def write_merged_bed(path: Path, rows: list[tuple[str, int, int, str, int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(rows)


def bed_record_count(path: Path) -> int:
    if not path.exists() or path.stat().st_size == 0:
        return 0
    with path.open() as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith("#"))


def faidx_regions_path(window_fasta: Path) -> Path:
    return window_fasta.with_name("candidate_windows.faidx.regions.txt")


def faidx_header_map_path(window_fasta: Path) -> Path:
    return window_fasta.with_name("candidate_windows.faidx.header_map.tsv")


def faidx_tmp_fasta_path(window_fasta: Path) -> Path:
    return window_fasta.with_name("candidate_windows.faidx.tmp.fa")


def read_merged_bed_records(merged_bed: Path) -> list[dict[str, str | int]]:
    records: list[dict[str, str | int]] = []
    if not merged_bed.exists() or merged_bed.stat().st_size == 0:
        return records

    with merged_bed.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            chrom = fields[0]
            start0 = int(fields[1])
            end0 = int(fields[2])
            labels = fields[3] if len(fields) > 3 and fields[3] else "."
            start1 = start0 + 1
            records.append({
                "chrom": chrom,
                "start0": start0,
                "start1": start1,
                "end0": end0,
                "region": f"{chrom}:{start1}-{end0}",
                "window_id": f"{chrom}_{start1}-{end0}:.",
                "labels": labels,
            })
    return records


def write_faidx_region_files(
    merged_bed: Path,
    regions_path: Path,
    header_map_path: Path,
) -> list[dict[str, str | int]]:
    records = read_merged_bed_records(merged_bed)
    regions_path.parent.mkdir(parents=True, exist_ok=True)
    with regions_path.open("w") as handle:
        for record in records:
            handle.write(f"{record['region']}\n")

    with header_map_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["region", "window_id", "labels"])
        for record in records:
            writer.writerow([record["region"], record["window_id"], record["labels"]])

    return records


def rewrite_faidx_headers(
    faidx_fasta: Path,
    window_fasta: Path,
    records: list[dict[str, str | int]],
) -> None:
    record_iter = iter(records)
    n_headers = 0
    window_fasta.parent.mkdir(parents=True, exist_ok=True)

    with faidx_fasta.open() as src, window_fasta.open("w") as dst:
        for line in src:
            if line.startswith(">"):
                try:
                    record = next(record_iter)
                except StopIteration as exc:
                    raise ValueError(
                        f"samtools faidx returned more FASTA records than expected: {faidx_fasta}"
                    ) from exc
                n_headers += 1
                window_id = str(record["window_id"])
                labels = str(record["labels"])
                suffix = f" {labels}" if labels and labels != "." else ""
                dst.write(f">{window_id}{suffix}\n")
            else:
                dst.write(line)

    if n_headers != len(records):
        raise ValueError(
            f"samtools faidx returned {n_headers} FASTA records, expected {len(records)}"
        )


def build_species_window_beds(
    species: str,
    output_dir: Path,
    genome_fasta: Path | None = None,
    nhmmer_gff: Path | None = None,
    miniprot_evidence_gff: Path | None = None,
    flank_bp: int = 50_000,
    force: bool = False,
) -> dict[str, str | int]:
    del force
    genome_fasta = genome_fasta or find_genome_fasta(species)
    fai_path = Path(f"{genome_fasta}.fai")
    if not fai_path.exists():
        raise FileNotFoundError(f"Missing FASTA index: {fai_path}")

    nhmmer_gff = nhmmer_gff or PROCESSED_DATA_DIR / "nhmmer_search_parsed" / species / f"{species}.gff"
    miniprot_evidence_gff = miniprot_evidence_gff or PROCESSED_DATA_DIR / "miniprot_output" / f"{species}.gff"
    if not miniprot_evidence_gff.exists():
        raise FileNotFoundError(f"Missing miniprot evidence GFF: {miniprot_evidence_gff}")

    output_dir.mkdir(parents=True, exist_ok=True)
    raw_bed = output_dir / "candidate_windows.raw.bed"
    merged_bed = output_dir / "candidate_windows.merged.bed"
    window_fasta = output_dir / "candidate_windows.fa"

    contig_lengths = read_fai(fai_path)
    nhmmer_intervals = parse_evidence_intervals(
        nhmmer_gff,
        "nhmmer",
        {"NTD", "CTD"},
        contig_lengths,
        flank_bp,
    )
    miniprot_intervals = parse_evidence_intervals(
        miniprot_evidence_gff,
        "miniprot",
        {"mRNA"},
        contig_lengths,
        flank_bp,
    )
    intervals = nhmmer_intervals + miniprot_intervals
    merged = merge_intervals(intervals)
    write_raw_bed(raw_bed, intervals)
    write_merged_bed(merged_bed, merged)

    return {
        "species": species,
        "genome_fasta": str(genome_fasta),
        "nhmmer_gff": str(nhmmer_gff) if nhmmer_gff.exists() else "",
        "miniprot_evidence_gff": str(miniprot_evidence_gff),
        "nhmmer_status": "present" if nhmmer_gff.exists() else "missing",
        "n_nhmmer_intervals": len(nhmmer_intervals),
        "n_miniprot_intervals": len(miniprot_intervals),
        "n_raw_intervals": len(intervals),
        "n_merged_windows": len(merged),
        "candidate_windows_raw_bed": str(raw_bed),
        "candidate_windows_merged_bed": str(merged_bed),
        "candidate_windows_fasta": str(window_fasta),
    }


def extract_species_window_fasta(
    genome_fasta: Path,
    merged_bed: Path,
    window_fasta: Path,
    faidx_threads: int = 4,
    force: bool = False,
    seqkit_threads: int | None = None,
) -> dict[str, str | int]:
    if seqkit_threads is not None:
        logger.warning("--seqkit-threads is deprecated; using it as --faidx-threads")
        faidx_threads = seqkit_threads

    n_merged_windows = bed_record_count(merged_bed)
    regions_path = faidx_regions_path(window_fasta)
    if n_merged_windows > 0:
        if window_fasta.exists() and not force:
            return {
                "genome_fasta": str(genome_fasta),
                "candidate_windows_merged_bed": str(merged_bed),
                "candidate_windows_fasta": str(window_fasta),
                "n_merged_windows": n_merged_windows,
                "extract_tool": "samtools faidx",
                "faidx_threads": faidx_threads,
                "faidx_regions": str(regions_path),
            }

        header_map = faidx_header_map_path(window_fasta)
        tmp_fasta = faidx_tmp_fasta_path(window_fasta)
        records = write_faidx_region_files(merged_bed, regions_path, header_map)
        cmd = (
            f"samtools faidx -@ {faidx_threads} -r {shell_path(regions_path)} "
            f"-o {shell_path(tmp_fasta)} {shell_path(genome_fasta)}"
        )
        run_cmd(cmd, [tmp_fasta], force=True)
        rewrite_faidx_headers(tmp_fasta, window_fasta, records)
        tmp_fasta.unlink(missing_ok=True)
    elif force or not window_fasta.exists():
        window_fasta.parent.mkdir(parents=True, exist_ok=True)
        window_fasta.write_text("")
        write_faidx_region_files(merged_bed, regions_path, faidx_header_map_path(window_fasta))

    return {
        "genome_fasta": str(genome_fasta),
        "candidate_windows_merged_bed": str(merged_bed),
        "candidate_windows_fasta": str(window_fasta),
        "n_merged_windows": n_merged_windows,
        "extract_tool": "samtools faidx",
        "faidx_threads": faidx_threads,
        "faidx_regions": str(regions_path),
    }

def prepare_species_windows(
    species: str,
    output_dir: Path,
    genome_fasta: Path | None = None,
    nhmmer_gff: Path | None = None,
    miniprot_evidence_gff: Path | None = None,
    flank_bp: int = 50_000,
    faidx_threads: int = 4,
    force: bool = False,
    seqkit_threads: int | None = None,
) -> dict[str, str | int]:
    if seqkit_threads is not None:
        logger.warning("--seqkit-threads is deprecated; using it as --faidx-threads")
        faidx_threads = seqkit_threads

    bed_row = build_species_window_beds(
        species,
        output_dir,
        genome_fasta=genome_fasta,
        nhmmer_gff=nhmmer_gff,
        miniprot_evidence_gff=miniprot_evidence_gff,
        flank_bp=flank_bp,
        force=force,
    )
    extract_row = extract_species_window_fasta(
        Path(str(bed_row["genome_fasta"])),
        Path(str(bed_row["candidate_windows_merged_bed"])),
        Path(str(bed_row["candidate_windows_fasta"])),
        faidx_threads=faidx_threads,
        force=force,
    )
    return {
        "species": species,
        "genome_fasta": bed_row["genome_fasta"],
        "nhmmer_gff": bed_row["nhmmer_gff"],
        "miniprot_evidence_gff": bed_row["miniprot_evidence_gff"],
        "nhmmer_status": bed_row["nhmmer_status"],
        "n_nhmmer_intervals": bed_row["n_nhmmer_intervals"],
        "n_miniprot_intervals": bed_row["n_miniprot_intervals"],
        "n_raw_intervals": bed_row["n_raw_intervals"],
        "n_merged_windows": bed_row["n_merged_windows"],
        "candidate_windows_fasta": bed_row["candidate_windows_fasta"],
        "extract_tool": extract_row["extract_tool"],
        "faidx_threads": extract_row["faidx_threads"],
        "faidx_regions": extract_row["faidx_regions"],
    }


@app.command()
def main(
    species: str,
    output_dir: Path,
    genome_fasta: Path | None = None,
    nhmmer_gff: Path | None = None,
    miniprot_evidence_gff: Path | None = None,
    flank_bp: int = typer.Option(50_000, "--flank-bp", min=0),
    faidx_threads: int = typer.Option(70, "--faidx-threads", min=1),
    seqkit_threads: int | None = typer.Option(None, "--seqkit-threads", min=1),
    force: bool = False,
) -> None:
    if seqkit_threads is not None:
        logger.warning("--seqkit-threads is deprecated; using it as --faidx-threads")
        faidx_threads = seqkit_threads

    logger.info(f"Preparing windows for {species}")
    row = prepare_species_windows(
        species,
        output_dir,
        genome_fasta=genome_fasta,
        nhmmer_gff=nhmmer_gff,
        miniprot_evidence_gff=miniprot_evidence_gff,
        flank_bp=flank_bp,
        faidx_threads=faidx_threads,
        force=force,
    )
    for _ in tqdm(range(1), total=1, desc="Window preparation"):
        pass
    logger.success(
        f"{species}: {row['n_merged_windows']} merged windows from {row['n_raw_intervals']} raw intervals"
    )


if __name__ == "__main__":
    app()
