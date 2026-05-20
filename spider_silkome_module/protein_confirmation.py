"""
protein_confirmation.py — parse miniprot GFF on spidroin_full_length.fasta,
translate predicted CDS into protein, and run pairwise alignment vs the
matching NCBI full-length spidroin protein.

Input layout (per species):
  miniprot_raw_dir / <species>.gff           (miniprot full-protein alignment)
  typing_dir       / <species> / spidroin_full_length.fasta
  query_fasta                                  (cd-hit deduplicated NCBI proteins)

Output layout (per species, written under output_dir/<species>/):
  fullprot_miniprot.gff          # cleaned GFF: gene/mRNA/exon/CDS/intron rows
  predicted_proteins.fa          # translated proteins (best hit per spidroin)
  pairwise_vs_ncbi.tsv           # one row per spidroin with pairwise metrics

The "best hit" per spidroin_id is the miniprot mRNA with the highest score on
that target sequence.  Pairwise alignment uses Biopython PairwiseAligner with
the BLOSUM62 matrix in global mode.
"""

from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from loguru import logger
import polars as pl
from tqdm import tqdm
import typer

from spider_silkome_module.config import INTERIM_DATA_DIR, PROCESSED_DATA_DIR

app = typer.Typer()

FASTA_FILE = "spidroin_full_length.fasta"
DEFAULT_TASK_NAME = "protein_confirmation"


# ── data structures ──────────────────────────────────────────────────────────


@dataclass
class MiniprotRecord:
    """One miniprot mRNA hit (a candidate prediction for a spidroin locus)."""

    spidroin_id: str        # target seqid (= entry of spidroin_full_length.fasta)
    query: str              # NCBI protein ID used as query
    rank: int               # miniprot rank (1 = best for this query×target)
    score: int              # miniprot score
    identity: float
    positive: float
    target_from: int        # query coords (1-based, inclusive)
    target_to: int
    strand: str             # almost always '+' (mRNA-oriented input)
    mrna_start: int         # on the spidroin_full_length sequence (1-based)
    mrna_end: int
    cds_blocks: list[tuple[int, int, int]] = field(default_factory=list)
    # list of (start, end, phase) — 1-based inclusive on the target sequence
    has_start: bool = False
    has_stop: bool = False


# ── GFF parsing ──────────────────────────────────────────────────────────────


def _parse_attrs(attr_str: str) -> dict[str, str]:
    """Parse GFF3 attribute field into a dict."""
    attrs: dict[str, str] = {}
    for part in attr_str.strip().split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            attrs[k.strip()] = v.strip()
    return attrs


def parse_miniprot_gff(gff_path: Path) -> dict[str, list[MiniprotRecord]]:
    """
    Parse a miniprot --gff file into MiniprotRecord objects keyed by spidroin_id
    (the target sequence id = a spidroin locus from spidroin_full_length.fasta).
    """
    records: dict[str, MiniprotRecord] = {}
    out: dict[str, list[MiniprotRecord]] = defaultdict(list)

    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, _src, feature, start_s, end_s, score_s, strand, phase_s, attrs_s = parts[:9]
            attrs = _parse_attrs(attrs_s)
            if feature == "mRNA":
                target_field = attrs.get("Target", "")
                target_parts = target_field.split()
                query = target_parts[0] if target_parts else ""
                tfrom = int(target_parts[1]) if len(target_parts) > 1 else 0
                tto = int(target_parts[2]) if len(target_parts) > 2 else 0
                rec = MiniprotRecord(
                    spidroin_id=seqid,
                    query=query,
                    rank=int(attrs.get("Rank", "0")),
                    score=int(float(score_s)) if score_s != "." else 0,
                    identity=float(attrs.get("Identity", "0")),
                    positive=float(attrs.get("Positive", "0")),
                    target_from=tfrom,
                    target_to=tto,
                    strand=strand,
                    mrna_start=int(start_s),
                    mrna_end=int(end_s),
                )
                records[attrs.get("ID", "")] = rec
                out[seqid].append(rec)
            elif feature == "CDS":
                parent = attrs.get("Parent", "")
                if parent in records:
                    phase = int(phase_s) if phase_s.isdigit() else 0
                    records[parent].cds_blocks.append((int(start_s), int(end_s), phase))
            elif feature == "start_codon":
                parent = attrs.get("Parent", "")
                if parent in records:
                    records[parent].has_start = True
            elif feature == "stop_codon":
                parent = attrs.get("Parent", "")
                if parent in records:
                    records[parent].has_stop = True

    # CDS blocks may arrive out of order — sort each record by start.
    for recs in out.values():
        for r in recs:
            r.cds_blocks.sort(key=lambda b: b[0])
    return out


# ── translation ──────────────────────────────────────────────────────────────


def translate_record(record: MiniprotRecord, target_seq: str) -> tuple[str, bool]:
    """
    Concatenate CDS DNA from the target sequence and translate to protein.
    Returns (protein_string, has_premature_stop).

    target_seq is the full spidroin_full_length sequence in mRNA orientation.
    """
    if not record.cds_blocks:
        return "", False

    if record.strand == "+":
        chunks = [target_seq[s - 1:e] for s, e, _ in record.cds_blocks]
    else:
        # In the rare case miniprot reports a '-' hit (target itself is mRNA-
        # oriented but query was matched to its complement), reverse-complement
        # the concatenation in CDS order from highest coordinate.
        chunks = [str(Seq(target_seq[s - 1:e]).reverse_complement())
                  for s, e, _ in reversed(record.cds_blocks)]

    cds = "".join(chunks)
    # Trim to a multiple of 3 to avoid Biopython warnings on partial codons.
    cds = cds[: len(cds) - (len(cds) % 3)]
    protein = str(Seq(cds).translate(to_stop=False))

    has_premature = "*" in protein.rstrip("*")
    # Strip trailing stop codons for downstream alignment / output.
    protein = protein.rstrip("*")
    return protein, has_premature


# ── pairwise alignment ──────────────────────────────────────────────────────


def _make_aligner() -> PairwiseAligner:
    """BLOSUM62 + global alignment with mild gap penalties."""
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.mode = "global"
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -1
    aligner.end_insertion_score = 0
    aligner.end_deletion_score = 0
    return aligner


def pairwise_metrics(
    predicted: str,
    reference: str,
    aligner: PairwiseAligner,
) -> tuple[float, float, float]:
    """
    Run a global pairwise alignment and return:
      identity (fraction of reference length),
      coverage_ref (aligned ref columns / ref length),
      coverage_query (aligned predicted columns / predicted length).

    "Aligned" = position is non-gap in BOTH sequences.  Identity is
    counted only over those columns.
    """
    if not predicted or not reference:
        return 0.0, 0.0, 0.0
    try:
        alignment = next(iter(aligner.align(predicted, reference)))
    except (ValueError, StopIteration):
        return 0.0, 0.0, 0.0

    aligned_p = alignment[0]  # gapped predicted
    aligned_r = alignment[1]  # gapped reference
    same = mismatch = 0
    for a, b in zip(aligned_p, aligned_r):
        if a == "-" or b == "-":
            continue
        if a == b:
            same += 1
        else:
            mismatch += 1
    aligned_cols = same + mismatch
    if aligned_cols == 0:
        return 0.0, 0.0, 0.0

    identity = same / aligned_cols
    coverage_ref = aligned_cols / max(1, len(reference))
    coverage_query = aligned_cols / max(1, len(predicted))
    return identity, coverage_ref, coverage_query


# ── per-species processing ───────────────────────────────────────────────────


def select_best_per_spidroin(
    records_per_id: dict[str, list[MiniprotRecord]],
) -> dict[str, MiniprotRecord]:
    """
    Pick the highest-scoring miniprot record per spidroin_id.
    """
    chosen: dict[str, MiniprotRecord] = {}
    for sid, recs in records_per_id.items():
        if not recs:
            continue
        best = max(recs, key=lambda r: (r.score, r.identity))
        chosen[sid] = best
    return chosen


def write_clean_gff(
    records: dict[str, MiniprotRecord],
    output_path: Path,
) -> None:
    """
    Write a tidy GFF3 (gene → mRNA → exon/CDS/intron) covering the best hit
    of each spidroin locus.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as out:
        out.write("##gff-version 3\n")
        out.write("# Source: protein_confirmation (best miniprot hit per spidroin)\n")
        for sid, r in sorted(records.items()):
            attr_gene = f"ID={sid};Name={sid};query={r.query}"
            out.write(
                f"{sid}\tprotein_confirmation\tgene\t{r.mrna_start}\t{r.mrna_end}"
                f"\t{r.score}\t{r.strand}\t.\t{attr_gene}\n"
            )
            mrna_id = f"{sid}.m1"
            attr_mrna = (
                f"ID={mrna_id};Parent={sid};query={r.query};"
                f"identity={r.identity};positive={r.positive};score={r.score}"
            )
            out.write(
                f"{sid}\tprotein_confirmation\tmRNA\t{r.mrna_start}\t{r.mrna_end}"
                f"\t{r.score}\t{r.strand}\t.\t{attr_mrna}\n"
            )
            for i, (s, e, ph) in enumerate(r.cds_blocks, 1):
                out.write(
                    f"{sid}\tprotein_confirmation\texon\t{s}\t{e}\t.\t"
                    f"{r.strand}\t.\tID={mrna_id}.exon{i};Parent={mrna_id}\n"
                )
                out.write(
                    f"{sid}\tprotein_confirmation\tCDS\t{s}\t{e}\t.\t"
                    f"{r.strand}\t{ph}\tID={mrna_id}.cds{i};Parent={mrna_id}\n"
                )
            for i in range(len(r.cds_blocks) - 1):
                istart = r.cds_blocks[i][1] + 1
                iend = r.cds_blocks[i + 1][0] - 1
                if iend >= istart:
                    out.write(
                        f"{sid}\tprotein_confirmation\tintron\t{istart}\t{iend}\t.\t"
                        f"{r.strand}\t.\tID={mrna_id}.intron{i + 1};Parent={mrna_id}\n"
                    )


def process_species(
    species_name: str,
    miniprot_gff: Path,
    spidroin_fasta: Path,
    query_fasta: Path,
    out_species_dir: Path,
    aligner: PairwiseAligner,
) -> int:
    """
    Run the full confirmation pipeline for one species.
    Returns the number of spidroins processed.
    """
    if not miniprot_gff.exists() or miniprot_gff.stat().st_size == 0:
        logger.warning(f"[{species_name}] no miniprot GFF, skipping")
        return 0
    if not spidroin_fasta.exists():
        logger.warning(f"[{species_name}] no spidroin fasta, skipping")
        return 0

    out_species_dir.mkdir(parents=True, exist_ok=True)

    target_seqs = {
        rec.id: str(rec.seq).upper()
        for rec in SeqIO.parse(spidroin_fasta, "fasta")
    }
    # The cd-hit deduplicated query fasta is loaded once at the caller; we read
    # it once here via the lazy dict on demand, but typically it's tens of MB.
    query_seqs = {
        rec.id: str(rec.seq).upper().rstrip("*")
        for rec in SeqIO.parse(query_fasta, "fasta")
    }

    records_per_id = parse_miniprot_gff(miniprot_gff)
    best = select_best_per_spidroin(records_per_id)

    rows: list[dict] = []
    proteins: list[SeqRecord] = []

    for sid, rec in best.items():
        target = target_seqs.get(sid)
        if target is None:
            logger.warning(f"[{species_name}] {sid} not in fasta, skipping")
            continue
        protein, premature = translate_record(rec, target)
        ref = query_seqs.get(rec.query, "")
        identity = coverage_q = coverage_t = 0.0
        if protein and ref:
            identity, coverage_q, coverage_t = pairwise_metrics(protein, ref, aligner)
        introns = [
            (rec.cds_blocks[i][1] + 1, rec.cds_blocks[i + 1][0] - 1)
            for i in range(len(rec.cds_blocks) - 1)
        ]
        intron_total = sum(max(0, e - s + 1) for s, e in introns)

        rows.append({
            "spidroin_id": sid,
            "species": species_name,
            "best_ncbi_match": rec.query,
            "miniprot_score": rec.score,
            "miniprot_identity": rec.identity,
            "pairwise_identity": round(identity, 4),
            "pairwise_coverage_query": round(coverage_q, 4),
            "pairwise_coverage_target": round(coverage_t, 4),
            "predicted_protein_length": len(protein),
            "ref_protein_length": len(ref),
            "exon_count": len(rec.cds_blocks),
            "intron_count": len(introns),
            "total_intron_len": intron_total,
            "has_premature_stop": premature,
            "miniprot_has_start": rec.has_start,
            "miniprot_has_stop": rec.has_stop,
            "mrna_start": rec.mrna_start,
            "mrna_end": rec.mrna_end,
            "strand": rec.strand,
        })
        if protein:
            proteins.append(SeqRecord(
                Seq(protein),
                id=sid,
                description=f"query={rec.query} score={rec.score} exons={len(rec.cds_blocks)}",
            ))

    if rows:
        df = pl.DataFrame(rows)
        df.write_csv(out_species_dir / "pairwise_vs_ncbi.tsv", separator="\t")
    SeqIO.write(proteins, out_species_dir / "predicted_proteins.fa", "fasta")
    write_clean_gff(best, out_species_dir / "fullprot_miniprot.gff")

    return len(rows)


# ── CLI ──────────────────────────────────────────────────────────────────────


@app.command()
def main(
    miniprot_raw_dir: Path = INTERIM_DATA_DIR / DEFAULT_TASK_NAME / "miniprot_raw",
    typing_dir: Path = PROCESSED_DATA_DIR / "typing_results",
    query_fasta: Path = INTERIM_DATA_DIR / DEFAULT_TASK_NAME / "cdhit" / "cdhit_shortest_seq.fa",
    output_dir: Path = PROCESSED_DATA_DIR / DEFAULT_TASK_NAME / "confirmation",
    species: str | None = None,
):
    """
    Run protein confirmation pipeline for one or all species.
    Pass --species "001.Allagelena_difficilis" to process a single species.
    """
    aligner = _make_aligner()
    output_dir.mkdir(parents=True, exist_ok=True)

    gff_files = sorted(miniprot_raw_dir.glob("*.gff"))
    if species is not None:
        gff_files = [g for g in gff_files if g.stem == species]
    if not gff_files:
        logger.error(f"No miniprot GFF files matching species filter in {miniprot_raw_dir}")
        raise typer.Exit(1)

    logger.info(f"Processing {len(gff_files)} species")
    total = 0
    for gff in tqdm(gff_files, desc="Confirming proteins"):
        species_name = gff.stem
        spidroin_fa = typing_dir / species_name / FASTA_FILE
        n = process_species(
            species_name,
            gff,
            spidroin_fa,
            query_fasta,
            output_dir / species_name,
            aligner,
        )
        logger.info(f"[{species_name}] processed {n} spidroins")
        total += n

    logger.success(f"Done. {total} spidroins confirmed across {len(gff_files)} species.")


if __name__ == "__main__":
    app()
