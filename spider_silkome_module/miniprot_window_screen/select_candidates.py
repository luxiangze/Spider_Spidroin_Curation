from __future__ import annotations

import csv
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any

from loguru import logger
from tqdm import tqdm
import typer

from spider_silkome_module.miniprot_window_screen.common import (
    Window,
    load_windows_bed,
    overlaps,
    parse_attrs,
    parse_float,
    parse_int,
)

app = typer.Typer(help="Select likely full-length spidroin MPIDs from parsed miniprot window models.")


@dataclass(frozen=True)
class DomainHit:
    hit_id: str
    window_id: str
    seqid: str
    feature: str
    start: int
    end: int
    score: float | None
    strand: str
    name: str
    evalue: str


@dataclass(frozen=True)
class DomainPair:
    pair_id: str
    window_id: str
    seqid: str
    strand: str
    ntd: DomainHit
    ctd: DomainHit

    @property
    def start(self) -> int:
        return min(self.ntd.start, self.ctd.start)

    @property
    def end(self) -> int:
        return max(self.ntd.end, self.ctd.end)


@dataclass(frozen=True)
class TypingLocus:
    spidroin_id: str
    window_id: str
    seqid: str
    start: int
    end: int
    strand: str
    spidroin_type: str
    full_length: bool
    hint_type: str
    confidence: str
    needs_review: bool
    ntd_start: int | None = None
    ntd_end: int | None = None
    ctd_start: int | None = None
    ctd_end: int | None = None

    @property
    def pair_id(self) -> str:
        return f"typing:{self.spidroin_id}"

    @property
    def pair_start(self) -> int:
        starts = [self.start]
        if self.ntd_start is not None:
            starts.append(self.ntd_start)
        if self.ctd_start is not None:
            starts.append(self.ctd_start)
        return min(starts)

    @property
    def pair_end(self) -> int:
        ends = [self.end]
        if self.ntd_end is not None:
            ends.append(self.ntd_end)
        if self.ctd_end is not None:
            ends.append(self.ctd_end)
        return max(ends)

    @property
    def terminal_support_count(self) -> int:
        ntd = self.ntd_start is not None and self.ntd_end is not None
        ctd = self.ctd_start is not None and self.ctd_end is not None
        return int(ntd) + int(ctd)


def load_models(models_jsonl: Path) -> list[dict[str, Any]]:
    models: list[dict[str, Any]] = []
    with models_jsonl.open() as handle:
        for line in handle:
            if line.strip():
                models.append(json.loads(line))
    return models


def gff_rows(gff_path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if not gff_path.exists():
        return rows
    with gff_path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            rows.append({
                "seqid": fields[0],
                "source": fields[1],
                "feature": fields[2],
                "start": int(fields[3]),
                "end": int(fields[4]),
                "score": parse_float(fields[5]),
                "strand": fields[6],
                "phase": fields[7],
                "attrs": parse_attrs(fields[8]),
            })
    return rows


def row_window_id(row: dict[str, Any], windows: list[Window]) -> str | None:
    for window in windows:
        if row["seqid"] == window.chrom and overlaps(row["start"], row["end"], window.start1, window.end0):
            return window.window_id
    return None


def parse_bool(value: Any) -> bool:
    return str(value or "").strip().lower() in {"1", "true", "yes", "y"}


def load_typing_loci_from_tsv(typing_tsv: Path, windows: list[Window]) -> list[TypingLocus]:
    loci: list[TypingLocus] = []
    if not typing_tsv.exists() or typing_tsv.stat().st_size == 0:
        return loci
    with typing_tsv.open() as handle:
        for i, row in enumerate(csv.DictReader(handle, delimiter="\t"), 1):
            seqid = row.get("Chr") or row.get("chrom") or row.get("seqid")
            start = parse_int(row.get("Start"))
            end = parse_int(row.get("End"))
            strand = row.get("Stran") or row.get("Strand") or "."
            if not seqid or start is None or end is None:
                continue
            window_id = row_window_id({"seqid": seqid, "start": start, "end": end}, windows)
            if window_id is None:
                continue
            loci.append(
                TypingLocus(
                    spidroin_id=row.get("Spidroin_ID") or f"typing_locus_{i:05d}",
                    window_id=window_id,
                    seqid=seqid,
                    start=start,
                    end=end,
                    strand=strand,
                    spidroin_type=row.get("Spidroin_type") or "",
                    full_length=parse_bool(row.get("Full_length")),
                    hint_type=row.get("Hint_type") or "",
                    confidence=row.get("confidence") or "",
                    needs_review=parse_bool(row.get("needs_review")),
                    ntd_start=parse_int(row.get("ntd_start")),
                    ntd_end=parse_int(row.get("ntd_end")),
                    ctd_start=parse_int(row.get("ctd_start")),
                    ctd_end=parse_int(row.get("ctd_end")),
                )
            )
    return loci


def load_typing_loci_from_gff(typing_gff: Path, windows: list[Window]) -> list[TypingLocus]:
    if not typing_gff.exists() or typing_gff.stat().st_size == 0:
        return []

    genes: dict[str, dict[str, Any]] = {}
    domains: dict[str, dict[str, tuple[int, int]]] = {}
    for row in gff_rows(typing_gff):
        attrs = row["attrs"]
        if row["feature"] == "spidroin_gene":
            spidroin_id = attrs.get("ID")
            if spidroin_id:
                genes[spidroin_id] = row
        elif row["feature"] in {"NTD", "CTD"}:
            parent = attrs.get("Parent")
            if parent:
                domains.setdefault(parent, {})[row["feature"]] = (row["start"], row["end"])

    loci: list[TypingLocus] = []
    for spidroin_id, row in genes.items():
        window_id = row_window_id(row, windows)
        if window_id is None:
            continue
        attrs = row["attrs"]
        ntd = domains.get(spidroin_id, {}).get("NTD")
        ctd = domains.get(spidroin_id, {}).get("CTD")
        hint_type = attrs.get("hint_type") or ""
        loci.append(
            TypingLocus(
                spidroin_id=spidroin_id,
                window_id=window_id,
                seqid=row["seqid"],
                start=row["start"],
                end=row["end"],
                strand=row["strand"],
                spidroin_type=attrs.get("Name") or "",
                full_length=hint_type == "Full_length",
                hint_type=hint_type,
                confidence=attrs.get("confidence") or "",
                needs_review=False,
                ntd_start=ntd[0] if ntd else None,
                ntd_end=ntd[1] if ntd else None,
                ctd_start=ctd[0] if ctd else None,
                ctd_end=ctd[1] if ctd else None,
            )
        )
    return loci


def load_typing_loci(
    typing_tsv: Path | None,
    typing_gff: Path | None,
    windows: list[Window],
) -> list[TypingLocus]:
    if typing_tsv and typing_tsv.exists():
        loci = load_typing_loci_from_tsv(typing_tsv, windows)
        if loci:
            return loci
    if typing_gff and typing_gff.exists():
        return load_typing_loci_from_gff(typing_gff, windows)
    return []


def load_domain_hits(nhmmer_gff: Path, windows: list[Window]) -> list[DomainHit]:
    hits: list[DomainHit] = []
    for i, row in enumerate(gff_rows(nhmmer_gff), 1):
        if row["feature"] not in {"NTD", "CTD"}:
            continue
        window_id = row_window_id(row, windows)
        if window_id is None:
            continue
        attrs = row["attrs"]
        hit_id = attrs.get("ID") or f"domain_{i}"
        hits.append(
            DomainHit(
                hit_id=hit_id,
                window_id=window_id,
                seqid=row["seqid"],
                feature=row["feature"],
                start=row["start"],
                end=row["end"],
                score=row["score"],
                strand=row["strand"],
                name=attrs.get("Name") or row["feature"],
                evalue=attrs.get("E-value") or "",
            )
        )
    return hits


def make_domain_pairs(
    hits: list[DomainHit],
    nearest_ctds_per_ntd: int = 2,
    max_pairs_per_window: int = 200,
) -> list[DomainPair]:
    by_window: dict[str, list[DomainHit]] = {}
    for hit in hits:
        by_window.setdefault(hit.window_id, []).append(hit)

    pairs: list[DomainPair] = []
    for window_id, window_hits in by_window.items():
        ntds = [hit for hit in window_hits if hit.feature == "NTD"]
        ctds = [hit for hit in window_hits if hit.feature == "CTD"]
        ranked_pairs: list[tuple[int, DomainHit, DomainHit]] = []
        for ntd in ntds:
            compatible: list[tuple[int, DomainHit]] = []
            for ctd in ctds:
                if ntd.seqid != ctd.seqid or ntd.strand != ctd.strand:
                    continue
                valid_order = (ntd.end <= ctd.start) if ntd.strand == "+" else (ctd.end <= ntd.start)
                if not valid_order:
                    continue
                distance = ctd.start - ntd.end if ntd.strand == "+" else ntd.start - ctd.end
                compatible.append((distance, ctd))
            for distance, ctd in sorted(compatible, key=lambda item: item[0])[:nearest_ctds_per_ntd]:
                ranked_pairs.append((distance, ntd, ctd))

        rank = 1
        for _distance, ntd, ctd in sorted(ranked_pairs, key=lambda item: item[0])[:max_pairs_per_window]:
            pairs.append(
                DomainPair(
                    pair_id=f"{window_id}|pair{rank:03d}|{ntd.hit_id}|{ctd.hit_id}",
                    window_id=window_id,
                    seqid=ntd.seqid,
                    strand=ntd.strand,
                    ntd=ntd,
                    ctd=ctd,
                )
            )
            rank += 1
    return pairs


def index_domain_hits(hits: list[DomainHit]) -> dict[str, list[DomainHit]]:
    by_window: dict[str, list[DomainHit]] = {}
    for hit in hits:
        by_window.setdefault(hit.window_id, []).append(hit)
    return by_window


def index_evidence_rows(evidence_rows: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    by_seqid: dict[str, list[dict[str, Any]]] = {}
    for row in evidence_rows:
        by_seqid.setdefault(row["seqid"], []).append(row)
    return by_seqid


def evidence_overlap_count(model: dict[str, Any], evidence_rows: list[dict[str, Any]]) -> int:
    n = 0
    for row in evidence_rows:
        if row["feature"] not in {"mRNA", "CDS"}:
            continue
        if row["seqid"] == model["genomic_seqid"] and overlaps(
            row["start"],
            row["end"],
            model["genomic_start"],
            model["genomic_end"],
        ):
            n += 1
    return n


def terminal_support(model: dict[str, Any], hits: list[DomainHit]) -> tuple[int, str]:
    names: list[str] = []
    for hit in hits:
        if (
            hit.window_id == model["window_id"]
            and hit.seqid == model["genomic_seqid"]
            and hit.strand == model["strand"]
            and overlaps(hit.start, hit.end, model["genomic_start"], model["genomic_end"])
        ):
            names.append(f"{hit.feature}:{hit.name}:{hit.hit_id}")
    return len(names), ",".join(sorted(names))


def overlap_len(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def typing_overlap_fraction(model: dict[str, Any], locus: TypingLocus) -> float:
    n = overlap_len(model["genomic_start"], model["genomic_end"], locus.start, locus.end)
    if n == 0:
        return 0.0
    model_len = max(1, model["genomic_end"] - model["genomic_start"] + 1)
    locus_len = max(1, locus.end - locus.start + 1)
    return n / min(model_len, locus_len)


def value(model: dict[str, Any], key: str, default: float = 0.0) -> float:
    parsed = parse_float(model.get(key))
    return default if parsed is None else parsed


def protein_len(model: dict[str, Any]) -> int:
    if model.get("protein_sequence_length") is not None:
        return int(model["protein_sequence_length"])
    seq = model.get("protein_sequence")
    return len(seq) if isinstance(seq, str) else 0


def is_strong(model: dict[str, Any], min_positive: float, min_query_coverage: float, min_aa: int) -> bool:
    return (
        value(model, "positive") >= min_positive
        and value(model, "query_coverage") >= min_query_coverage
        and protein_len(model) >= min_aa
    )


def is_rescue(
    model: dict[str, Any],
    min_positive: float,
    min_query_coverage: float,
    min_aa: int,
    has_domain_support: bool,
) -> bool:
    return (
        has_domain_support
        and value(model, "positive") >= min_positive
        and value(model, "query_coverage") >= min_query_coverage
        and protein_len(model) >= min_aa
    )


def selection_score(
    model: dict[str, Any],
    domain_pair: DomainPair | None,
    terminal_count: int,
    evidence_count: int,
    typing_locus: TypingLocus | None = None,
    typing_overlap_frac: float = 0.0,
) -> float:
    exon_count = int(model.get("exon_count") or 0)
    exon_bonus = 8 if 1 <= exon_count <= 5 else 0
    exon_penalty = 10 if exon_count > 10 else 0
    domain_bonus = 25 if domain_pair else 8 * min(terminal_count, 2)
    if typing_locus:
        domain_bonus = max(domain_bonus, 25 if typing_locus.full_length else 8 * typing_locus.terminal_support_count)
    evidence_bonus = min(evidence_count, 5) * 3
    typing_bonus = typing_overlap_frac * 30 if typing_locus else 0
    confidence_bonus = 5 if typing_locus and typing_locus.confidence == "high" else 0
    review_penalty = 12 if typing_locus and typing_locus.needs_review else 0
    return (
        value(model, "positive") * 100
        + value(model, "query_coverage") * 80
        + min(protein_len(model), 5000) / 100
        + value(model, "score") / 1000
        + domain_bonus
        + evidence_bonus
        + exon_bonus
        + typing_bonus
        + confidence_bonus
        - exon_penalty
        - review_penalty
    )


def typing_support_hits(locus: TypingLocus | None) -> str:
    if locus is None:
        return ""
    hits = []
    if locus.ntd_start is not None and locus.ntd_end is not None:
        hits.append(f"typing:NTD:{locus.spidroin_id}")
    if locus.ctd_start is not None and locus.ctd_end is not None:
        hits.append(f"typing:CTD:{locus.spidroin_id}")
    return ",".join(hits)


def selected_row(
    species: str,
    model: dict[str, Any],
    domain_pair: DomainPair | None,
    terminal_count: int,
    terminal_hits: str,
    evidence_count: int,
    status: str,
    selection_mode: str,
    typing_locus: TypingLocus | None = None,
    typing_overlap_frac: float = 0.0,
) -> dict[str, Any]:
    domain_pair_id = domain_pair.pair_id if domain_pair else ""
    domain_pair_start = domain_pair.start if domain_pair else ""
    domain_pair_end = domain_pair.end if domain_pair else ""
    if typing_locus:
        domain_pair_id = typing_locus.pair_id
        domain_pair_start = typing_locus.pair_start
        domain_pair_end = typing_locus.pair_end
        terminal_count = max(terminal_count, typing_locus.terminal_support_count)
        extra_hits = typing_support_hits(typing_locus)
        terminal_hits = ",".join(item for item in [terminal_hits, extra_hits] if item)

    return {
        "species": species,
        "window_id": model["window_id"],
        "mpid": model["mrna_id"],
        "selection_status": status,
        "selection_mode": selection_mode,
        "typing_spidroin_id": typing_locus.spidroin_id if typing_locus else "",
        "typing_spidroin_type": typing_locus.spidroin_type if typing_locus else "",
        "typing_hint_type": typing_locus.hint_type if typing_locus else "",
        "typing_full_length": typing_locus.full_length if typing_locus else "",
        "typing_confidence": typing_locus.confidence if typing_locus else "",
        "typing_needs_review": typing_locus.needs_review if typing_locus else "",
        "typing_overlap_frac": typing_overlap_frac if typing_locus else "",
        "domain_pair_id": domain_pair_id,
        "domain_pair_start": domain_pair_start,
        "domain_pair_end": domain_pair_end,
        "domain_support_count": terminal_count,
        "domain_support_hits": terminal_hits,
        "miniprot_evidence_overlap_count": evidence_count,
        "genomic_seqid": model["genomic_seqid"],
        "genomic_start": model["genomic_start"],
        "genomic_end": model["genomic_end"],
        "strand": model["strand"],
        "positive": model.get("positive"),
        "identity": model.get("identity"),
        "query_coverage": model.get("query_coverage"),
        "score": model.get("score"),
        "protein_sequence_length": protein_len(model),
        "exon_count": model.get("exon_count"),
        "intron_count": model.get("intron_count"),
        "target_protein": model.get("target_protein"),
        "selection_score": selection_score(
            model,
            domain_pair,
            terminal_count,
            evidence_count,
            typing_locus=typing_locus,
            typing_overlap_frac=typing_overlap_frac,
        ),
    }


def ranked_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return sorted(
        rows,
        key=lambda row: (
            -float(row["selection_score"]),
            -float(row["positive"] or 0),
            -float(row["query_coverage"] or 0),
            -int(row["protein_sequence_length"] or 0),
            str(row["mpid"]),
        ),
    )


def select_with_typing_loci(
    species: str,
    models_by_window: dict[str, list[dict[str, Any]]],
    typing_loci: list[TypingLocus],
    hits_by_window: dict[str, list[DomainHit]],
    evidence_by_seqid: dict[str, list[dict[str, Any]]],
    strong_positive: float,
    strong_query_coverage: float,
    rescue_positive: float,
    rescue_query_coverage: float,
    min_aa: int,
) -> list[dict[str, Any]]:
    selected: list[dict[str, Any]] = []
    selected_mpids: set[str] = set()
    sorted_loci = sorted(typing_loci, key=lambda item: (item.window_id, item.start, item.spidroin_id))

    for locus in tqdm(sorted_loci, desc="Select typing-guided candidates"):
        ranked: list[dict[str, Any]] = []
        for model in models_by_window.get(locus.window_id, []):
            if model["genomic_seqid"] != locus.seqid:
                continue
            if locus.strand in {"+", "-"} and model["strand"] != locus.strand:
                continue
            overlap_frac = typing_overlap_fraction(model, locus)
            if overlap_frac <= 0:
                continue
            terminal_count, terminal_hits = terminal_support(model, hits_by_window.get(model["window_id"], []))
            evidence_count = evidence_overlap_count(model, evidence_by_seqid.get(model["genomic_seqid"], []))
            has_domain_support = terminal_count > 0 or locus.terminal_support_count > 0 or locus.full_length
            if is_strong(model, strong_positive, strong_query_coverage, min_aa):
                status = "selected_strong"
            elif is_rescue(model, rescue_positive, rescue_query_coverage, min_aa, has_domain_support):
                status = "selected_rescue"
            else:
                continue
            ranked.append(
                selected_row(
                    species,
                    model,
                    None,
                    terminal_count,
                    terminal_hits,
                    evidence_count,
                    status,
                    "typing_guided",
                    typing_locus=locus,
                    typing_overlap_frac=overlap_frac,
                )
            )
        for row in ranked_rows(ranked):
            if row["mpid"] not in selected_mpids:
                selected.append(row)
                selected_mpids.add(row["mpid"])
                break

    return selected


def select_with_fallback_pairs(
    species: str,
    windows: list[Window],
    models_by_window: dict[str, list[dict[str, Any]]],
    domain_hits: list[DomainHit],
    domain_pairs: list[DomainPair],
    hits_by_window: dict[str, list[DomainHit]],
    evidence_by_seqid: dict[str, list[dict[str, Any]]],
    strong_positive: float,
    strong_query_coverage: float,
    rescue_positive: float,
    rescue_query_coverage: float,
    min_aa: int,
) -> list[dict[str, Any]]:
    selected: list[dict[str, Any]] = []
    selected_keys: set[tuple[str, str]] = set()

    for pair in tqdm(domain_pairs, desc="Select paired-domain candidates"):
        candidates = [
            model for model in models_by_window.get(pair.window_id, [])
            if model["genomic_seqid"] == pair.seqid
            and model["strand"] == pair.strand
            and overlaps(model["genomic_start"], model["genomic_end"], pair.start, pair.end)
        ]
        ranked: list[dict[str, Any]] = []
        for model in candidates:
            terminal_count, terminal_hits = terminal_support(model, hits_by_window.get(model["window_id"], []))
            evidence_count = evidence_overlap_count(model, evidence_by_seqid.get(model["genomic_seqid"], []))
            if is_strong(model, strong_positive, strong_query_coverage, min_aa):
                status = "selected_strong"
            elif is_rescue(model, rescue_positive, rescue_query_coverage, min_aa, True):
                status = "selected_rescue"
            else:
                continue
            ranked.append(
                selected_row(
                    species,
                    model,
                    pair,
                    terminal_count,
                    terminal_hits,
                    evidence_count,
                    status,
                    "fallback",
                )
            )
        if ranked:
            best = ranked_rows(ranked)[0]
            key = (best["window_id"], best["mpid"])
            if key not in selected_keys:
                selected.append(best)
                selected_keys.add(key)

    paired_windows = {pair.window_id for pair in domain_pairs}
    for window in tqdm(windows, desc="Select unpaired-window candidates"):
        window_models = models_by_window.get(window.window_id, [])
        if not window_models:
            continue
        ranked = []
        for model in window_models:
            terminal_count, terminal_hits = terminal_support(model, hits_by_window.get(model["window_id"], []))
            evidence_count = evidence_overlap_count(model, evidence_by_seqid.get(model["genomic_seqid"], []))
            has_domain = terminal_count > 0
            if is_strong(model, strong_positive, strong_query_coverage, min_aa):
                status = "selected_strong"
            elif is_rescue(model, rescue_positive, rescue_query_coverage, min_aa, has_domain):
                status = "selected_rescue"
            else:
                continue
            ranked.append(
                selected_row(
                    species,
                    model,
                    None,
                    terminal_count,
                    terminal_hits,
                    evidence_count,
                    status,
                    "fallback",
                )
            )
        if not ranked:
            continue
        max_rows = 1 if window.window_id not in paired_windows else 3
        for row in ranked_rows(ranked)[:max_rows]:
            key = (row["window_id"], row["mpid"])
            if key not in selected_keys:
                selected.append(row)
                selected_keys.add(key)

    return selected


def select_candidates(
    species: str,
    models: list[dict[str, Any]],
    windows: list[Window],
    domain_hits: list[DomainHit],
    domain_pairs: list[DomainPair],
    evidence_rows: list[dict[str, Any]],
    typing_loci: list[TypingLocus] | None = None,
    strong_positive: float = 0.60,
    strong_query_coverage: float = 0.80,
    rescue_positive: float = 0.55,
    rescue_query_coverage: float = 0.70,
    min_aa: int = 1000,
) -> list[dict[str, Any]]:
    models_by_window: dict[str, list[dict[str, Any]]] = {}
    for model in models:
        models_by_window.setdefault(model["window_id"], []).append(model)

    hits_by_window = index_domain_hits(domain_hits)
    evidence_by_seqid = index_evidence_rows(evidence_rows)

    if typing_loci:
        selected = select_with_typing_loci(
            species,
            models_by_window,
            typing_loci,
            hits_by_window,
            evidence_by_seqid,
            strong_positive,
            strong_query_coverage,
            rescue_positive,
            rescue_query_coverage,
            min_aa,
        )
    else:
        selected = select_with_fallback_pairs(
            species,
            windows,
            models_by_window,
            domain_hits,
            domain_pairs,
            hits_by_window,
            evidence_by_seqid,
            strong_positive,
            strong_query_coverage,
            rescue_positive,
            rescue_query_coverage,
            min_aa,
        )

    return sorted(selected, key=lambda row: (row["species"], row["window_id"], -float(row["selection_score"])))


def write_rows(rows: list[dict[str, Any]], output_tsv: Path) -> None:
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
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
        "exon_count",
        "intron_count",
        "target_protein",
        "selection_score",
    ]
    with output_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


@app.command()
def main(
    species: str,
    models_jsonl: Path,
    windows_bed: Path,
    output_tsv: Path,
    nhmmer_gff: Path | None = None,
    miniprot_evidence_gff: Path | None = None,
    typing_tsv: Path | None = typer.Option(None, "--typing-tsv"),
    typing_gff: Path | None = typer.Option(None, "--typing-gff"),
    strong_positive: float = typer.Option(0.60, "--strong-positive"),
    strong_query_coverage: float = typer.Option(0.80, "--strong-query-coverage"),
    rescue_positive: float = typer.Option(0.55, "--rescue-positive"),
    rescue_query_coverage: float = typer.Option(0.70, "--rescue-query-coverage"),
    min_aa: int = typer.Option(1000, "--min-aa"),
) -> None:
    logger.info(f"Selecting candidates for {species}")
    models = load_models(models_jsonl)
    windows = load_windows_bed(windows_bed)
    domain_hits = load_domain_hits(nhmmer_gff, windows) if nhmmer_gff and nhmmer_gff.exists() else []
    typing_loci = load_typing_loci(typing_tsv, typing_gff, windows)
    domain_pairs = [] if typing_loci else make_domain_pairs(domain_hits)
    evidence_rows = gff_rows(miniprot_evidence_gff) if miniprot_evidence_gff and miniprot_evidence_gff.exists() else []
    logger.info(
        f"selection_mode={'typing_guided' if typing_loci else 'fallback'}; "
        f"typing_loci={len(typing_loci)}; domain_pairs={len(domain_pairs)}; models={len(models)}"
    )
    rows = select_candidates(
        species,
        models,
        windows,
        domain_hits,
        domain_pairs,
        evidence_rows,
        typing_loci=typing_loci,
        strong_positive=strong_positive,
        strong_query_coverage=strong_query_coverage,
        rescue_positive=rescue_positive,
        rescue_query_coverage=rescue_query_coverage,
        min_aa=min_aa,
    )
    write_rows(rows, output_tsv)
    logger.success(f"Selected {len(rows)} MPIDs for {species}: {output_tsv}")


if __name__ == "__main__":
    app()
