"""Unit + integration tests for the k-mer repeat-profile intron finder."""

import random

from spider_silkome_module.find_introns_kmer import (
    Intron,
    coverage_block,
    find_introns,
    flank_copy,
    kmer_min_copy_profile,
    parse_locus,
    parse_type,
    plot_coverage,
    refine_boundaries,
    rna_coverage,
    single_copy_runs,
    splice_out_introns,
    translate_cds,
)

# ── kmer_min_copy_profile ─────────────────────────────────────────────────────


def test_profile_homopolymer():
    # Every 3-mer of "AAAAAA" is "AAA" (count 4), so every base has copy 4.
    assert kmer_min_copy_profile("AAAAAA", k=3) == [4, 4, 4, 4, 4, 4]


def test_profile_short_sequence_is_single_copy():
    assert kmer_min_copy_profile("ACGT", k=14) == [1, 1, 1, 1]


def test_profile_matches_naive_implementation():
    from collections import Counter

    def naive(seq, k):
        n = len(seq)
        cov = [n + 1] * n
        counts = Counter(seq[i:i + k] for i in range(n - k + 1))
        for i in range(n - k + 1):
            c = counts[seq[i:i + k]]
            for j in range(i, i + k):
                cov[j] = min(cov[j], c)
        return cov

    rng = random.Random(0)
    seq = "".join(rng.choice("ACGT") for _ in range(2000))
    assert kmer_min_copy_profile(seq, k=14) == naive(seq, 14)


# ── single_copy_runs ──────────────────────────────────────────────────────────


def test_single_copy_runs_filters_short_runs():
    profile = [1, 1, 1, 5, 5, 1, 1, 1, 1, 1]
    # 1-based runs of copy==1 with length >= 3
    assert single_copy_runs(profile, min_run=3) == [(1, 3), (6, 10)]


def test_single_copy_runs_ignores_below_threshold():
    profile = [1, 1, 5, 5, 1, 1, 1]
    assert single_copy_runs(profile, min_run=3) == [(5, 7)]


# ── flank_copy ────────────────────────────────────────────────────────────────


def test_flank_copy_median():
    profile = [5, 5, 5, 1, 1, 1, 1, 7, 7, 7]
    # candidate run is the 1s at 1-based 4..7
    left, right = flank_copy(profile, cand_start=4, cand_end=7, window=3)
    assert (left, right) == (5, 7)


# ── refine_boundaries ─────────────────────────────────────────────────────────


def test_refine_boundaries_exact():
    # "C"*30 + GT + "C"*42 + AG + "C"*30 : only one GT and one AG.
    seq = "C" * 30 + "GT" + "C" * 42 + "AG" + "C" * 30
    # true intron 1-based: G at 31 .. G of AG at 76, length 46 (46 % 3 == 1)
    result = refine_boundaries(seq, cand_start=31, cand_end=76)
    assert result == (31, 76, 46, 0)


def test_refine_boundaries_rejects_wrong_frame():
    # Only candidate pair has length 47 (47 % 3 == 2) -> rejected for remainder 1.
    seq = "C" * 30 + "GT" + "C" * 43 + "AG" + "C" * 30
    assert refine_boundaries(seq, cand_start=31, cand_end=77) is None
    # ...but accepted when we ask for remainder 2.
    assert refine_boundaries(seq, cand_start=31, cand_end=77, frame_remainder=2) is not None


# ── _resolve_overlaps via find_introns helper behaviour ───────────────────────


def test_resolve_overlaps_keeps_smaller_span_diff():
    from spider_silkome_module.find_introns_kmer import _resolve_overlaps

    a = Intron(100, 160, 61, "GT", "AG", 100, 160, 5, 5, 8)
    b = Intron(140, 200, 61, "GT", "AG", 140, 200, 5, 5, 2)  # overlaps a, better diff
    c = Intron(300, 360, 61, "GT", "AG", 300, 360, 5, 5, 1)  # disjoint
    kept = _resolve_overlaps([a, b, c])
    assert [(i.start, i.end) for i in kept] == [(140, 200), (300, 360)]


# ── find_introns integration ──────────────────────────────────────────────────


def _overlap(a_start, a_end, b_start, b_end):
    inter = max(0, min(a_end, b_end) - max(a_start, b_start) + 1)
    union = (a_end - a_start + 1) + (b_end - b_start + 1) - inter
    return inter / union


def test_find_introns_recovers_synthetic_intron():
    rng = random.Random(42)

    def uniq(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    ntd = uniq(120)
    unit = uniq(40)
    intron = "GT" + uniq(90) + "AG"  # length 94 (94 % 3 == 1)
    ctd = uniq(120)
    # intron embedded in a high-copy repeat body (8 unit copies -> copy ~8)
    seq = ntd + unit * 4 + intron + unit * 4 + ctd
    true_start = len(ntd) + len(unit * 4) + 1
    true_end = true_start + len(intron) - 1

    # NTD, intron candidate, CTD -> exactly three single-copy runs
    runs = single_copy_runs(kmer_min_copy_profile(seq, 14), min_run=30)
    assert len(runs) == 3

    introns = find_introns(seq, min_run=30, min_flank_cov=3)
    assert len(introns) == 1
    it = introns[0]
    assert it.donor == "GT" and it.acceptor == "AG"
    assert it.length % 3 == 1
    # k-mer edge effect can shift boundaries by ~k; require strong overlap.
    assert _overlap(it.start, it.end, true_start, true_end) > 0.6


def test_find_introns_empty_without_repeat_body():
    rng = random.Random(7)
    seq = "".join(rng.choice("ACGT") for _ in range(2000))  # no repeats -> one run
    assert find_introns(seq, min_run=30) == []


# ── coverage_block ─────────────────────────────────────────────────────────────


def test_coverage_block_format():
    # 3 columns per base: position (1-based), base, coverage; '# id' header line.
    block = coverage_block("seqX", "ACGT", [1, 2, 3, 4])
    assert block == "# seqX\n1\tA\t1\n2\tC\t2\n3\tG\t3\n4\tT\t4\n"


# ── plot_coverage ──────────────────────────────────────────────────────────────


def test_plot_coverage_writes_png(tmp_path):
    profile = [1] * 30 + [8] * 60 + [1] * 30  # NTD valley, repeat peak, CTD valley
    introns = [Intron(50, 110, 61, "GT", "AG", 50, 110, 8, 8, 2)]
    out = tmp_path / "seqX.png"
    plot_coverage("seqX", profile, introns, out, ntd_ctd=((1, 30), (91, 120)))
    assert out.exists() and out.stat().st_size > 0


# ── RNA-seq bigWig ─────────────────────────────────────────────────────────────


def test_parse_locus():
    assert parse_locus("Aldi_1 chr07:100-200(+) type=AcSp") == ("chr07", 100, 200, "+")
    assert parse_locus("x scaf3:5-9(-) foo") == ("scaf3", 5, 9, "-")
    assert parse_locus("no coordinates here") is None


def test_rna_coverage_strand(tmp_path):
    import pyBigWig

    p = tmp_path / "t.bw"
    bw = pyBigWig.open(str(p), "w")
    bw.addHeader([("chr1", 100)])
    bw.addEntries("chr1", 0, values=[float(i) for i in range(20)], span=1, step=1)
    bw.close()

    reader = pyBigWig.open(str(p))
    # chr1:1-10(+) -> genomic 0-based [0, 10) = 0..9
    assert rna_coverage(reader, "s chr1:1-10(+)", 10) == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    # '-' strand reverses to align with the sense sequence
    assert rna_coverage(reader, "s chr1:1-10(-)", 10) == [9, 8, 7, 6, 5, 4, 3, 2, 1, 0]
    assert rna_coverage(reader, "s chr1:1-10(+)", 11) is None  # length mismatch
    assert rna_coverage(reader, "s chrX:1-10(+)", 10) is None  # chromosome absent
    reader.close()


def test_coverage_block_with_rna():
    block = coverage_block("seqX", "ACGT", [1, 5, 5, 1], rna=[0.0, 12.0, float("nan"), 3.0])
    assert block == "# seqX\n1\tA\t1\t0\n2\tC\t5\t12\n3\tG\t5\tNA\n4\tT\t1\t3\n"


def test_plot_coverage_with_rna(tmp_path):
    profile = [1] * 30 + [8] * 60 + [1] * 30
    introns = [Intron(50, 110, 61, "GT", "AG", 50, 110, 8, 8, 2)]
    rna = [0.0] * 30 + [50.0] * 60 + [0.0] * 30  # RNA drops to 0 over the terminal valleys
    out = tmp_path / "rna.png"
    plot_coverage("seqX", profile, introns, out, ntd_ctd=((1, 30), (91, 120)), rna=rna)
    assert out.exists() and out.stat().st_size > 0


# ── parse_type ─────────────────────────────────────────────────────────────────


def test_parse_type():
    assert parse_type("Evsp_spid_00001 chr03:1-2(-) type=MiSp") == "MiSp"
    assert parse_type("x chr1:1-2(+) type=MaSp1") == "MaSp1"
    assert parse_type("no type field here") is None


# ── splice_out_introns ─────────────────────────────────────────────────────────


def test_splice_out_introns_removes_interval():
    # 1-based inclusive 4..6 of "AAACCCGGG" is "CCC" -> "AAAGGG".
    it = Intron(4, 6, 3, "GT", "AG", 4, 6, 5, 5, 0)
    assert splice_out_introns("AAACCCGGG", [it]) == "AAAGGG"


def test_splice_out_introns_multiple_and_empty():
    introns = [
        Intron(4, 6, 3, "GT", "AG", 4, 6, 5, 5, 0),    # removes "CCC"
        Intron(10, 12, 3, "GT", "AG", 10, 12, 5, 5, 0),  # removes "TTT"
    ]
    assert splice_out_introns("AAACCCGGGTTT", introns) == "AAAGGG"
    assert splice_out_introns("AAACCCGGG", []) == "AAACCCGGG"  # no introns -> unchanged


# ── translate_cds ──────────────────────────────────────────────────────────────


def test_translate_cds_clean_orf_has_no_premature_stop():
    # ATG AAA TAA -> "MK*" -> terminal stop stripped, no internal stop.
    protein, premature = translate_cds("ATGAAATAA")
    assert protein == "MK"
    assert premature is False


def test_translate_cds_detects_premature_stop():
    # ATG TAA AAA TAA -> "M*K*" -> internal stop remains after stripping terminal.
    protein, premature = translate_cds("ATGTAAAAATAA")
    assert premature is True
    assert "*" in protein


def test_translate_cds_trims_partial_codon():
    # Trailing "GG" is not a whole codon and is dropped before translation.
    protein, _ = translate_cds("ATGAAAGG")
    assert protein == "MK"
