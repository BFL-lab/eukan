"""Unit tests for eukan.assembly.bam_diagnostic.

Builds tiny in-memory BAMs + small FASTA fixtures to exercise:
- Soft-clip extraction on forward/reverse reads (leading + trailing).
- mRNA-orientation normalization (reverse-complement + side-label swap).
- Anchored-substring cluster key.
- N cigar-op intron walking and dinucleotide lookup against a genome.
- MAPQ + secondary/supplementary filtering.
- diagnose_bam() aggregation of loci and clusters.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pysam
import pytest

from eukan.assembly.bam_diagnostic import (
    DiagnosticReport,
    IntronStats,
    LocusConsistencyStats,
    SoftClipStats,
    _cluster_consensus,
    _cluster_key,
    _dinucleotide,
    _extract_clips,
    _hamming_distance,
    _is_low_complexity,
    _reverse_complement,
    _walk_introns,
    compute_verdict,
    diagnose_bam,
)
from eukan.infra.genome import ContigIndex


def _write_fasta(path: Path, contigs: list[tuple[str, str]]) -> Path:
    with open(path, "w") as f:
        for name, seq in contigs:
            f.write(f">{name}\n{seq}\n")
    return path


def _write_bam(
    path: Path,
    contigs: list[tuple[str, int]],
    reads: list[dict],
) -> Path:
    """Write a minimal BAM file.

    Each read dict requires: query_name, query_sequence, flag,
    reference_id, reference_start, cigartuples, mapping_quality.
    """
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": name, "LN": length} for name, length in contigs],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for r in reads:
            a = pysam.AlignedSegment(out.header)
            a.query_name = r["query_name"]
            a.query_sequence = r["query_sequence"]
            a.query_qualities = pysam.qualitystring_to_array(
                "I" * len(r["query_sequence"])
            )
            a.flag = r["flag"]
            a.reference_id = r["reference_id"]
            a.reference_start = r["reference_start"]
            a.mapping_quality = r["mapping_quality"]
            a.cigartuples = r["cigartuples"]
            out.write(a)
    return path


def _read_single(bam_path: Path, query_name: str) -> pysam.AlignedSegment:
    with pysam.AlignmentFile(str(bam_path), "rb") as f:
        for r in f:
            if r.query_name == query_name:
                return r
    raise KeyError(query_name)


# ------------------------------------------------------------------
# Pure helpers
# ------------------------------------------------------------------


def test_reverse_complement_basic():
    assert _reverse_complement("ACGTN") == "NACGT"
    assert _reverse_complement("aaccgg") == "ccggtt"


def test_cluster_key_anchors_correct_end():
    # 5p clip: anchored at the alignment-touching END → take last K
    assert _cluster_key("5p", "AAAAGGGGGGGGGGGG", k=8) == "GGGGGGGG"
    # 3p clip: anchored at the alignment-touching START → take first K
    assert _cluster_key("3p", "AAAAGGGGGGGGGGGG", k=8) == "AAAAGGGG"
    # Short clip falls back to full sequence
    assert _cluster_key("5p", "ACGT", k=8) == "ACGT"
    # Uppercased
    assert _cluster_key("5p", "aaaaggggggggggg", k=6) == "GGGGGG"


# ------------------------------------------------------------------
# Soft-clip extraction across the 4 strand x side cases
# ------------------------------------------------------------------


@pytest.fixture
def clip_bam(tmp_path):
    contigs = [("chr1", 300)]
    # Construct a clean 5p clip = "AAAACCCCAAAAGGGGGGGG" (20 bp); aligned
    # portion "TTTTTTTTTT" (10 M); for forward reads at ref 50,
    # leading clip is 5p in mRNA orientation.
    fwd_seq_5p = "AAAACCCCAAAAGGGGGGGG" + "TTTTTTTTTT"
    # Same logic, trailing clip on forward = 3p.
    fwd_seq_3p = "TTTTTTTTTT" + "AAAACCCCAAAAGGGGGGGG"
    # For reverse reads, query_sequence in BAM is RC of the original.
    # leading clip on reverse = 3p in mRNA; the mRNA-orient clip seq is
    # RC(BAM clip).
    rev_seq_3p_lead = "AAAACCCCAAAAGGGGGGGG" + "TTTTTTTTTT"
    rev_seq_5p_trail = "TTTTTTTTTT" + "AAAACCCCAAAAGGGGGGGG"

    reads = [
        # Forward, leading clip → 5p, seq as-is
        dict(
            query_name="fwd_lead",
            query_sequence=fwd_seq_5p,
            query_qualities=None,
            flag=0,
            reference_id=0,
            reference_start=50,
            mapping_quality=60,
            cigartuples=[(4, 20), (0, 10)],
        ),
        # Forward, trailing clip → 3p, seq as-is
        dict(
            query_name="fwd_trail",
            query_sequence=fwd_seq_3p,
            query_qualities=None,
            flag=0,
            reference_id=0,
            reference_start=100,
            mapping_quality=60,
            cigartuples=[(0, 10), (4, 20)],
        ),
        # Reverse, leading clip → 3p, seq RC'd
        dict(
            query_name="rev_lead",
            query_sequence=rev_seq_3p_lead,
            query_qualities=None,
            flag=16,
            reference_id=0,
            reference_start=150,
            mapping_quality=60,
            cigartuples=[(4, 20), (0, 10)],
        ),
        # Reverse, trailing clip → 5p, seq RC'd
        dict(
            query_name="rev_trail",
            query_sequence=rev_seq_5p_trail,
            query_qualities=None,
            flag=16,
            reference_id=0,
            reference_start=200,
            mapping_quality=60,
            cigartuples=[(0, 10), (4, 20)],
        ),
    ]
    return _write_bam(tmp_path / "clips.bam", contigs, reads)


def test_extract_clips_forward_leading_is_5p(clip_bam):
    r = _read_single(clip_bam, "fwd_lead")
    clips = list(_extract_clips(r, min_clip_len=8))
    assert len(clips) == 1
    side, seq, anchor = clips[0]
    assert side == "5p"
    assert seq == "AAAACCCCAAAAGGGGGGGG"
    assert anchor == 50


def test_extract_clips_forward_trailing_is_3p(clip_bam):
    r = _read_single(clip_bam, "fwd_trail")
    clips = list(_extract_clips(r, min_clip_len=8))
    assert len(clips) == 1
    side, seq, anchor = clips[0]
    assert side == "3p"
    assert seq == "AAAACCCCAAAAGGGGGGGG"
    # reference_end is exclusive of last aligned base; anchor = end - 1
    assert anchor == 100 + 10 - 1


def test_extract_clips_reverse_leading_is_3p_and_rcd(clip_bam):
    r = _read_single(clip_bam, "rev_lead")
    clips = list(_extract_clips(r, min_clip_len=8))
    assert len(clips) == 1
    side, seq, anchor = clips[0]
    assert side == "3p"
    # mRNA-orient clip = RC(BAM clip "AAAACCCCAAAAGGGGGGGG") = "CCCCCCCCTTTTGGGGTTTT"
    assert seq == _reverse_complement("AAAACCCCAAAAGGGGGGGG")
    assert anchor == 150


def test_extract_clips_reverse_trailing_is_5p_and_rcd(clip_bam):
    r = _read_single(clip_bam, "rev_trail")
    clips = list(_extract_clips(r, min_clip_len=8))
    assert len(clips) == 1
    side, seq, anchor = clips[0]
    assert side == "5p"
    assert seq == _reverse_complement("AAAACCCCAAAAGGGGGGGG")
    assert anchor == 200 + 10 - 1


def test_extract_clips_below_min_length_excluded(clip_bam):
    r = _read_single(clip_bam, "fwd_lead")
    assert list(_extract_clips(r, min_clip_len=21)) == []  # clip is 20 bp


# ------------------------------------------------------------------
# Intron walking and dinucleotide lookup
# ------------------------------------------------------------------


def test_walk_introns_single_N_op(tmp_path):
    # Read: 10M 50N 10M starting at ref 100 → intron at (100+10, 100+10+50)
    bam = _write_bam(
        tmp_path / "intron.bam",
        [("chr1", 500)],
        [
            dict(
                query_name="spliced",
                query_sequence="A" * 20,
                flag=0,
                reference_id=0,
                reference_start=100,
                mapping_quality=60,
                cigartuples=[(0, 10), (3, 50), (0, 10)],
            ),
        ],
    )
    r = _read_single(bam, "spliced")
    introns = list(_walk_introns(r))
    assert introns == [("chr1", 110, 160)]


def test_dinucleotide_forward_canonical(tmp_path):
    # Contig with GT...AG intron at coords 10-19 (10 bp).
    seq = "N" * 10 + "GTAAAAAAAG" + "N" * 30  # intron [10:20), donor "GT", acceptor "AG"
    fa = _write_fasta(tmp_path / "g.fa", [("chr1", seq)])
    with ContigIndex(fa) as contigs:
        assert _dinucleotide(contigs, "chr1", 10, 20) == "GT-AG"


def test_dinucleotide_reverse_strand_canonical(tmp_path):
    # CT...AC on + strand = reverse-strand canonical
    seq = "N" * 10 + "CTAAAAAAAC" + "N" * 30
    fa = _write_fasta(tmp_path / "g.fa", [("chr1", seq)])
    with ContigIndex(fa) as contigs:
        assert _dinucleotide(contigs, "chr1", 10, 20) == "CT-AC"


def test_dinucleotide_out_of_bounds_returns_none(tmp_path):
    seq = "ACGT" * 10
    fa = _write_fasta(tmp_path / "g.fa", [("chr1", seq)])
    with ContigIndex(fa) as contigs:
        assert _dinucleotide(contigs, "chr1", 100, 200) is None
        assert _dinucleotide(contigs, "missing", 0, 10) is None


# ------------------------------------------------------------------
# diagnose_bam aggregation
# ------------------------------------------------------------------


@pytest.fixture
def trans_splice_like_bam(tmp_path):
    """3 forward reads at distinct loci, all sharing the same 5p clip suffix."""
    contigs = [("chr1", 1000)]
    # Anchored 12 bp suffix is "GGGGGGGGGGGG" — should cluster 3 loci into 1
    clip1 = "AAAATTTTGGGGGGGGGGGG"  # 20 bp
    clip2 = "CCCCAAAGGGGGGGGGGGGG"  # 20 bp (different prefix, but last 12 bp = same "GGGGGGGGGGGG"... wait let me check)
    # Actually wait — clip2 last 12: "GGGGGGGGGGGG"? Let me count:
    #   "CCCCAAAGGGGGGGGGGGGG" has 20 chars. Last 12 = chars 8..19 = "GGGGGGGGGGGG" — yes.
    # And a third clip in a different cluster:
    clip3 = "TTTTAAAAACCCCCCCCCCCC"  # 21 bp; last 12 = "CCCCCCCCCCCC"

    reads = []
    for i, clip in enumerate([clip1, clip2, clip3]):
        seq = clip + "AAAAAAAAAA"  # 10 bp aligned
        reads.append(
            dict(
                query_name=f"r{i}",
                query_sequence=seq,
                flag=0,
                reference_id=0,
                reference_start=100 + i * 200,
                mapping_quality=60,
                cigartuples=[(4, len(clip)), (0, 10)],
            )
        )
    return _write_bam(tmp_path / "trans.bam", contigs, reads)


@pytest.fixture
def synthetic_genome(tmp_path):
    fa = _write_fasta(tmp_path / "g.fa", [("chr1", "N" * 1000)])
    return fa


def test_diagnose_bam_trans_splicing_signature(
    trans_splice_like_bam, synthetic_genome,
):
    report = diagnose_bam(
        trans_splice_like_bam, synthetic_genome,
        min_clip_len=8, cluster_key_len=12, min_mapq=20,
    )
    soft, intron = report.softclip, report.intron
    assert soft.n_reads_scanned == 3
    assert soft.n_reads_with_clip == 3
    assert soft.n_clips_total == 3
    assert soft.n_clips_by_side == {"5p": 3}
    assert soft.n_loci == 3
    # 2 unique cluster keys: "GGGGGGGGGGGG" (2 loci) and "CCCCCCCCCCCC" (1 locus)
    assert soft.n_clusters == 2
    assert soft.cluster_to_loci["GGGGGGGGGGGG"] == 2
    assert soft.cluster_to_loci["CCCCCCCCCCCC"] == 1
    # No N ops → no introns
    assert intron.n_introns_total == 0


def test_diagnose_bam_tallies_polya_without_perturbing_sl_verdict(tmp_path, synthetic_genome):
    # Reads whose only soft-clip is a poly-A 3' tail: counted in report.polya, but
    # poly-A is low-complexity so it must NOT register as a trans-splicing (SL) call.
    contigs = [("chr1", 1000)]
    reads = [
        dict(
            query_name=f"r{i}",
            query_sequence="ACGTACGTAC" + "A" * 20,  # 10 M aligned + 20S poly-A tail
            flag=0,
            reference_id=0,
            reference_start=100 + i * 50,
            mapping_quality=60,
            cigartuples=[(0, 10), (4, 20)],
        )
        for i in range(5)
    ]
    bam = _write_bam(tmp_path / "polya.bam", contigs, reads)
    report = diagnose_bam(bam, synthetic_genome, min_clip_len=8, min_mapq=20)

    assert report.polya.n_polya == 5
    assert report.polya.n_polyt == 0
    assert report.polya.n_clips_examined == 5
    verdict = compute_verdict(report)
    assert verdict.trans_splicing.call == "ABSENT"  # poly-A is not an SL signal


def test_diagnose_bam_filters_low_mapq_and_secondary(tmp_path, synthetic_genome):
    contigs = [("chr1", 1000)]
    clip = "AAAACCCCAAAAGGGGGGGG"
    seq = clip + "AAAAAAAAAA"
    reads = [
        # Primary, high MAPQ → kept
        dict(
            query_name="kept",
            query_sequence=seq,
            flag=0,
            reference_id=0,
            reference_start=100,
            mapping_quality=60,
            cigartuples=[(4, 20), (0, 10)],
        ),
        # Primary, low MAPQ → skipped
        dict(
            query_name="low_mapq",
            query_sequence=seq,
            flag=0,
            reference_id=0,
            reference_start=200,
            mapping_quality=5,
            cigartuples=[(4, 20), (0, 10)],
        ),
        # Secondary alignment → skipped (FLAG 0x100 = 256)
        dict(
            query_name="secondary",
            query_sequence=seq,
            flag=256,
            reference_id=0,
            reference_start=300,
            mapping_quality=60,
            cigartuples=[(4, 20), (0, 10)],
        ),
        # Unmapped → skipped (FLAG 0x4 = 4)
        dict(
            query_name="unmapped",
            query_sequence=seq,
            flag=4,
            reference_id=0,
            reference_start=400,
            mapping_quality=60,
            cigartuples=[(4, 20), (0, 10)],
        ),
    ]
    bam = _write_bam(tmp_path / "f.bam", contigs, reads)
    soft = diagnose_bam(bam, synthetic_genome, min_mapq=20).softclip
    assert soft.n_reads_scanned == 1


def test_diagnose_bam_introns_and_dinuc(tmp_path):
    # Contig: filler + GT intron + filler + AG ... AC + filler
    # Intron #1 at [50, 100): donor at 50-51 = "GT", acceptor at 98-99 = "AG"
    seq = (
        "N" * 50
        + "GT" + "A" * 46 + "AG"  # intron 1
        + "N" * 50
        + "CT" + "T" * 46 + "AC"  # intron 2 (reverse canonical)
        + "N" * 50
    )
    fa = _write_fasta(tmp_path / "g.fa", [("chr1", seq)])
    contigs = [("chr1", len(seq))]
    # Read 1: 10M 50N 10M starting at ref 40  → intron at (50, 100)  → "GT-AG"
    # Read 2: 10M 50N 10M starting at ref 140 → intron at (150, 200) → "CT-AC"
    reads = [
        dict(
            query_name="r1",
            query_sequence="A" * 20,
            flag=0,
            reference_id=0,
            reference_start=40,
            mapping_quality=60,
            cigartuples=[(0, 10), (3, 50), (0, 10)],
        ),
        dict(
            query_name="r2",
            query_sequence="A" * 20,
            flag=16,
            reference_id=0,
            reference_start=140,
            mapping_quality=60,
            cigartuples=[(0, 10), (3, 50), (0, 10)],
        ),
    ]
    bam = _write_bam(tmp_path / "i.bam", contigs, reads)
    report = diagnose_bam(bam, fa, min_clip_len=8)
    soft, intron = report.softclip, report.intron
    assert soft.n_clips_total == 0  # no clips on these reads
    assert intron.n_introns_total == 2
    assert intron.by_dinucleotide == {"GT-AG": 1, "CT-AC": 1}
    assert intron.canonical_pct == pytest.approx(100.0)


def test_diagnose_bam_top_clusters_sorted_by_locus_count(
    trans_splice_like_bam, synthetic_genome,
):
    soft = diagnose_bam(
        trans_splice_like_bam, synthetic_genome,
        min_clip_len=8, cluster_key_len=12, top_n=10,
    ).softclip
    keys_in_order = [k for k, _n_loci, _n_reads in soft.top_clusters]
    assert keys_in_order[0] == "GGGGGGGGGGGG"
    assert soft.top_clusters[0][1] == 2  # n_loci
    assert soft.top_clusters[1] == ("CCCCCCCCCCCC", 1, 1)


# ------------------------------------------------------------------
# Per-locus consolidation analysis (LocusConsistencyStats)
# ------------------------------------------------------------------


def _read_with_5p_clip(name: str, ref_start: int, clip: str) -> dict:
    """Build a forward read with a leading soft clip + 10 M aligned bases."""
    return dict(
        query_name=name,
        query_sequence=clip + "AAAAAAAAAA",
        flag=0,
        reference_id=0,
        reference_start=ref_start,
        mapping_quality=60,
        cigartuples=[(4, len(clip)), (0, 10)],
    )


def test_locus_consistent_with_varying_clip_lengths(tmp_path, synthetic_genome):
    """3 reads at one locus, clip lengths 10/14/18 — all share the same anchored suffix."""
    # Anchored 12 bp (last 12 of the longer clips) = "TTAAGGGGGGGG"
    # Shorter 10-bp clip is a suffix of that → consistent
    contigs = [("chr1", 1000)]
    bam = _write_bam(
        tmp_path / "consistent.bam",
        contigs,
        [
            _read_with_5p_clip("a", 100, "AAGGGGGGGG"),           # 10 bp (short)
            _read_with_5p_clip("b", 100, "TTTTAAGGGGGGGG"),       # 14 bp (full key)
            _read_with_5p_clip("c", 100, "CCCCTTTTAAGGGGGGGG"),   # 18 bp (full key)
        ],
    )
    report = diagnose_bam(bam, synthetic_genome, min_clip_len=8, cluster_key_len=12)
    lc = report.locus_consistency
    assert lc.n_loci_total == 1
    assert lc.n_loci_consistent == 1
    assert lc.n_loci_inconsistent == 0
    assert lc.n_loci_singleton == 0
    # 3 clips on this one locus; motif is the K-bp anchored suffix
    assert lc.all_rows[0].n_clips == 3
    assert lc.all_rows[0].status == "consistent"
    assert lc.all_rows[0].motif_key == "TTAAGGGGGGGG"


def test_locus_inconsistent_when_long_keys_differ(tmp_path, synthetic_genome):
    """2 reads at one locus with clearly different anchored substrings."""
    bam = _write_bam(
        tmp_path / "inconsistent.bam",
        [("chr1", 1000)],
        [
            _read_with_5p_clip("a", 100, "AAAATTTTGGGGGGGG"),  # last 12 = TTTTGGGGGGGG
            _read_with_5p_clip("b", 100, "AAAATTTTCCCCCCCC"),  # last 12 = TTTTCCCCCCCC
        ],
    )
    report = diagnose_bam(bam, synthetic_genome, min_clip_len=8, cluster_key_len=12)
    lc = report.locus_consistency
    assert lc.n_loci_total == 1
    assert lc.n_loci_inconsistent == 1
    assert lc.n_loci_consistent == 0


def test_locus_short_only_when_all_clips_below_K(tmp_path, synthetic_genome):
    """2 reads at one locus, both clips < K; counted as short-only, not consistent."""
    bam = _write_bam(
        tmp_path / "short_only.bam",
        [("chr1", 1000)],
        [
            _read_with_5p_clip("a", 100, "AAAAAAAA"),  # 8 bp (< K=12)
            _read_with_5p_clip("b", 100, "AAAAAAAA"),  # same clip
        ],
    )
    report = diagnose_bam(bam, synthetic_genome, min_clip_len=8, cluster_key_len=12)
    lc = report.locus_consistency
    assert lc.n_loci_total == 1
    assert lc.n_loci_short_only == 1
    assert lc.n_loci_consistent == 0
    assert lc.n_loci_inconsistent == 0
    # short-only loci are still emitted to all_rows (status="short_only")
    assert lc.all_rows[0].status == "short_only"
    assert lc.all_rows[0].motif_key == ""
    assert lc.all_rows[0].motif_share == 0


def test_motif_share_histogram_buckets(tmp_path, synthetic_genome):
    """2 loci sharing a motif + 1 locus with a unique motif → buckets '2-10' and '1'."""
    # Locus A (chr1:100) and B (chr1:200): same anchored 12 bp "TTTTGGGGGGGG"
    # Locus C (chr1:300): unique anchored 12 bp "CCCCGGGGGGGG"
    bam = _write_bam(
        tmp_path / "histogram.bam",
        [("chr1", 1000)],
        [
            _read_with_5p_clip("a1", 100, "AAAATTTTGGGGGGGG"),
            _read_with_5p_clip("a2", 100, "GGGGTTTTGGGGGGGG"),
            _read_with_5p_clip("b1", 200, "CCCCTTTTGGGGGGGG"),
            _read_with_5p_clip("b2", 200, "TTTTTTTTGGGGGGGG"),
            _read_with_5p_clip("c1", 300, "AAAACCCCGGGGGGGG"),
            _read_with_5p_clip("c2", 300, "GGGGCCCCGGGGGGGG"),
        ],
    )
    report = diagnose_bam(bam, synthetic_genome, min_clip_len=8, cluster_key_len=12)
    lc = report.locus_consistency
    assert lc.n_loci_total == 3
    assert lc.n_loci_consistent == 3
    # A and B share the same cluster_key (share=2), C is unique (share=1)
    assert lc.motif_share_histogram == {"2-10": 2, "1": 1}


def test_singleton_locus_counted_separately(tmp_path, synthetic_genome):
    """A 1-read locus is consistent by definition, but goes in n_loci_singleton."""
    bam = _write_bam(
        tmp_path / "singleton.bam",
        [("chr1", 1000)],
        [
            _read_with_5p_clip("a", 100, "AAAATTTTGGGGGGGG"),
        ],
    )
    report = diagnose_bam(bam, synthetic_genome, min_clip_len=8, cluster_key_len=12)
    lc = report.locus_consistency
    assert lc.n_loci_total == 1
    assert lc.n_loci_singleton == 1
    assert lc.n_loci_consistent == 0
    # singletons are not emitted to all_rows
    assert lc.all_rows == []


# ------------------------------------------------------------------
# Hamming-tolerant per-locus consistency + cross-locus canonicalization
# ------------------------------------------------------------------


def test_low_complexity_helper():
    """`_is_low_complexity` should reject poly-runs but pass real motifs."""
    assert _is_low_complexity("AAAAAAAAAAAA")
    assert _is_low_complexity("TTTTTTTTTT")
    assert _is_low_complexity("AAAAAAAAAACC")
    assert not _is_low_complexity("CTGTACTTTATT")   # SL motif, T=5/12
    assert not _is_low_complexity("AGAACATGGTCG")
    assert not _is_low_complexity("")               # empty stays False


def test_cluster_consensus_extends_beyond_key_5p():
    """Clips that all share a 24-bp 5p motif yield the full 24-bp consensus."""
    motif = "GAGAGAGAGAGACTGTACTTTATT"  # 24 bp; last 12 = anchor key
    clips = [motif] * 30
    assert _cluster_consensus(clips, "5p") == motif


def test_cluster_consensus_extends_beyond_key_3p():
    """Mirror of the 5p test: clips share a 3p motif starting at the anchor."""
    motif = "AATAAACATGTCAGAGAGAGAGAG"  # 24 bp; first 12 = anchor key
    clips = [motif] * 30
    assert _cluster_consensus(clips, "3p") == motif


def test_cluster_consensus_truncates_at_low_coverage_5p():
    """When only a few clips extend past the anchor, consensus stops at the anchor."""
    long_motif = "GAGAGAGAGAGACTGTACTTTATT"  # 24 bp
    anchor = "CTGTACTTTATT"                   # 12 bp (just the anchor)
    # 30 anchor-only + 3 long. min_coverage=5 → column past anchor has only 3.
    clips = [anchor] * 30 + [long_motif] * 3
    assert _cluster_consensus(clips, "5p") == anchor


def test_cluster_consensus_breaks_on_split_majority():
    """A 50/50 base split at a column terminates the consensus before that column."""
    # 10 clips with 'A' at the extended position, 10 with 'T'; anchor agreed
    # on all. min_majority_fraction=0.6 → the split column fails.
    clips_a = ["ACTGTACTTTATT"] * 10  # 13 bp; col -13 = 'A'
    clips_t = ["TCTGTACTTTATT"] * 10  # 13 bp; col -13 = 'T'
    consensus = _cluster_consensus(clips_a + clips_t, "5p")
    assert consensus == "CTGTACTTTATT"  # just the agreed 12-bp anchor


def test_cluster_consensus_empty_input():
    """Empty input returns empty string, no exception."""
    assert _cluster_consensus([], "5p") == ""
    assert _cluster_consensus([], "3p") == ""


def test_cluster_consensus_min_majority_fraction_relaxed():
    """Lowering ``min_majority_fraction`` lets the consensus extend past a noisy column."""
    # Last 12 bp of both clips = "CTGTACTTTATT" (anchor key, agreed).
    # Differ at idx 11 (col 12): clip_a has 'A', clip_t has 'T'.
    # Cols 13-23 agree again (the "AGAGAGAGAG" prefix).
    clip_a = "AGAGAGAGAGAACTGTACTTTATT"
    clip_t = "AGAGAGAGAGATCTGTACTTTATT"
    clips = [clip_a] * 6 + [clip_t] * 5  # col 12 majority 6/11 = 54.5%

    # Default 0.6 → col 12 fails, consensus stops at the anchor.
    assert _cluster_consensus(clips, "5p") == "CTGTACTTTATT"

    # 0.4 → col 12 passes (54.5% > 40%); consensus continues to the
    # 100%-agreed prefix and recovers the full clip.
    relaxed = _cluster_consensus(clips, "5p", min_majority_fraction=0.4)
    assert relaxed == clip_a


def test_diagnose_bam_consensus_excludes_off_anchor_longest_clip(
    tmp_path, synthetic_genome,
):
    """The cluster consensus ignores per-locus ``longest_clip`` when its
    K-bp anchor doesn't match the cluster's raw keys.

    Setup: 6 SL-only loci with 18-bp clips ending in CTGTACTTTATT, plus
    5 mixed loci where each carries both a 12-bp SL read (just the
    anchor key) and a 30-bp contamination read whose anchor is
    ``A`` x 12. The contamination wins ``longest_clip`` at each mixed
    locus. Without the K-bp anchor filter, contamination bases would
    fight SL bases at column 0 (6 vs 5 → 54.5% < 60%) and the consensus
    would collapse. With the filter the contamination is dropped and
    the consensus recovers the full 18-bp motif.
    """
    sl_clip = "AATCAACTGTACTTTATT"          # 18 bp, suffix = SL key
    sl_short = "CTGTACTTTATT"                # 12 bp, just the key
    contam_clip = "T" * 18 + "A" * 12        # 30 bp, suffix = AAAAA…

    reads: list[dict] = []
    # 6 SL-only loci (one SL read each, well-separated positions)
    for i in range(6):
        reads.append(_read_with_5p_clip(f"sl{i}", 100 + i * 50, sl_clip))
    # 5 mixed loci, each with both an SL anchor-only read and a long
    # contamination read sharing the same alignment start
    for i in range(5):
        pos = 1000 + i * 50
        reads.append(_read_with_5p_clip(f"mix_sl{i}", pos, sl_short))
        reads.append(_read_with_5p_clip(f"mix_co{i}", pos, contam_clip))

    bam = _write_bam(tmp_path / "filter.bam", [("chr1", 2000)], reads)
    report = diagnose_bam(
        bam, synthetic_genome, min_clip_len=8, cluster_key_len=12,
    )
    soft = report.softclip
    assert soft.cluster_consensus.get("CTGTACTTTATT", "") == sl_clip


def _make_report(
    top_clusters: list[tuple[str, int, int]],
    canonical_pct: float,
    by_dinucleotide: dict[str, int] | None = None,
    n_introns_total: int = 1000,
    n_loci_consistent: int = 100,
    motif_share_histogram: dict[str, int] | None = None,
    cluster_consensus: dict[str, str] | None = None,
) -> DiagnosticReport:
    """Build a synthetic DiagnosticReport for verdict-only unit tests."""
    soft = SoftClipStats(
        n_reads_scanned=10_000, n_reads_with_clip=1_000,
        n_clips_total=1_000, n_clips_by_side={"5p": 1_000, "3p": 0},
        n_loci=500, n_clusters=len(top_clusters),
        cluster_to_loci={k: n_loci for k, n_loci, _ in top_clusters},
        cluster_to_reads={k: n_reads for k, _, n_reads in top_clusters},
        cluster_examples={k: k for k, _, _ in top_clusters},
        cluster_consensus=cluster_consensus or {},
        top_clusters=top_clusters,
    )
    intron = IntronStats(
        n_introns_total=n_introns_total,
        by_dinucleotide=by_dinucleotide or {"GT-AG": int(n_introns_total * canonical_pct / 100)},
        canonical_pct=canonical_pct,
    )
    lc = LocusConsistencyStats(
        n_loci_total=500, n_loci_singleton=400,
        n_loci_consistent=n_loci_consistent,
        n_loci_inconsistent=0, n_loci_short_only=0,
        motif_share_histogram=motif_share_histogram or {},
    )
    return DiagnosticReport(softclip=soft, intron=intron, locus_consistency=lc)


def test_compute_verdict_trans_splicing_strong():
    """Top cluster ≥1K loci, ≥10K reads + canonical_pct ≥ 95% → STRONG / ABSENT."""
    report = _make_report(
        top_clusters=[("ATCGATCGATCG", 2000, 20_000)],
        canonical_pct=99.0,
        cluster_consensus={"ATCGATCGATCG": "GGGGGGATCGATCGATCGATCG"},
    )
    verdict = compute_verdict(report)
    assert verdict.trans_splicing.call == "STRONG"
    assert verdict.trans_splicing.top_non_trivial_cluster_key == "ATCGATCGATCG"
    assert (
        verdict.trans_splicing.top_non_trivial_cluster_consensus
        == "GGGGGGATCGATCGATCGATCG"
    )
    assert verdict.non_canonical_splice.call == "ABSENT"


def test_compute_verdict_non_canonical_extensive():
    """Top cluster too small + canonical_pct < 80% → ABSENT / EXTENSIVE."""
    report = _make_report(
        top_clusters=[("ATCGATCGATCG", 20, 500)],  # below MODERATE threshold
        canonical_pct=55.0,
        by_dinucleotide={"GT-AG": 300, "CT-AC": 250, "CT-CG": 200, "CT-CT": 100},
    )
    verdict = compute_verdict(report)
    assert verdict.trans_splicing.call == "ABSENT"
    assert verdict.non_canonical_splice.call == "EXTENSIVE"
    # Top non-canonical dinuc reported with its percentage
    assert "CT-CG" in verdict.non_canonical_splice.top_non_canonical_dinuc


def test_compute_verdict_skips_low_complexity_top_cluster():
    """First non-trivial cluster is picked, even if a poly-run is more populous."""
    report = _make_report(
        top_clusters=[
            ("TTTTTTTTTTTT", 5000, 50_000),    # low-complexity, must be skipped
            ("ATCGATCGATCG", 500, 5_000),      # MODERATE level
        ],
        canonical_pct=99.0,
    )
    verdict = compute_verdict(report)
    assert verdict.trans_splicing.top_non_trivial_cluster_key == "ATCGATCGATCG"
    assert verdict.trans_splicing.call == "MODERATE"


def test_hamming_distance_helper():
    assert _hamming_distance("ACGT", "ACGT") == 0
    assert _hamming_distance("ACGT", "ACGA") == 1
    assert _hamming_distance("ACGT", "AAAA") == 3
    # Different lengths sort as "far"
    assert _hamming_distance("ACGT", "ACG") == 4


def test_locus_consistent_under_hamming_tolerance(tmp_path, synthetic_genome):
    """3+ reads at one locus with K-bp keys differing by 1 base from a dominant key."""
    # Dominant K-bp suffix = "TTTTGGGGGGGG" (3 reads → wins seed selection)
    # Variant 1: differs at position 11 — TTTTGGGGGGGT (1 read)
    # Variant 2: differs at position 9  — TTTTGGGGGAGG (1 read)
    # All within H=1 of dominant. Default hamming_tolerance=1 → consistent.
    reads = []
    for i in range(3):
        reads.append(_read_with_5p_clip(f"dom{i}", 100, "AAAATTTTGGGGGGGG"))
    reads.append(_read_with_5p_clip("v1", 100, "AAAATTTTGGGGGGGT"))
    reads.append(_read_with_5p_clip("v2", 100, "AAAATTTTGGGGGAGG"))
    bam = _write_bam(tmp_path / "tolerant.bam", [("chr1", 1000)], reads)
    report = diagnose_bam(bam, synthetic_genome, min_clip_len=8, cluster_key_len=12)
    lc = report.locus_consistency
    assert lc.n_loci_total == 1
    assert lc.n_loci_consistent == 1
    assert lc.n_loci_inconsistent == 0
    # Dominant K-bp key wins (highest read count → first seed)
    assert lc.all_rows[0].motif_key == "TTTTGGGGGGGG"


def test_locus_strict_inconsistent_under_tolerance_zero(tmp_path, synthetic_genome):
    """Same 3 H=1 variants — with tolerance=0 they're 3 distinct keys → inconsistent."""
    reads = []
    for i in range(3):
        reads.append(_read_with_5p_clip(f"dom{i}", 100, "AAAATTTTGGGGGGGG"))
    reads.append(_read_with_5p_clip("v1", 100, "AAAATTTTGGGGGGGT"))
    reads.append(_read_with_5p_clip("v2", 100, "AAAATTTTGGGGGAGG"))
    bam = _write_bam(tmp_path / "strict.bam", [("chr1", 1000)], reads)
    report = diagnose_bam(
        bam, synthetic_genome, min_clip_len=8, cluster_key_len=12,
        hamming_tolerance=0,
    )
    lc = report.locus_consistency
    assert lc.n_loci_inconsistent == 1
    assert lc.n_loci_consistent == 0


def test_canonical_merging_across_loci(tmp_path, synthetic_genome):
    """3 loci with K-bp keys all within H=1 of a shared seed → canonical share = 3."""
    # Seed (4 reads, processed first by descending count): key = "AAAACCCCGGGG"
    # Locus B (2 reads): key = "AAAACCCCGGGT" — H=1 from seed (last base flipped G→T)
    # Locus C (2 reads): key = "AAAACCCCGGCG" — H=1 from seed (position 10 flipped G→C)
    # Pairwise B-vs-C is H=2, but both are H=1 from seed → both absorb into seed.
    reads = []
    # Locus A: 4 identical reads
    for i in range(4):
        reads.append(_read_with_5p_clip(f"a{i}", 100, "TTTTAAAACCCCGGGG"))
    # Locus B: 2 identical reads
    for i in range(2):
        reads.append(_read_with_5p_clip(f"b{i}", 200, "TTTTAAAACCCCGGGT"))
    # Locus C: 2 identical reads
    for i in range(2):
        reads.append(_read_with_5p_clip(f"c{i}", 300, "TTTTAAAACCCCGGCG"))

    bam = _write_bam(tmp_path / "canonical.bam", [("chr1", 1000)], reads)
    report = diagnose_bam(bam, synthetic_genome, min_clip_len=8, cluster_key_len=12)
    lc = report.locus_consistency
    soft = report.softclip
    # All 3 loci collapse to a single canonical cluster → share = 3 each
    assert lc.n_loci_total == 3
    assert lc.n_loci_consistent == 3
    assert lc.motif_share_histogram == {"2-10": 3}
    # The canonical seed is the highest-coverage raw key
    assert soft.n_clusters == 1
    assert "AAAACCCCGGGG" in soft.cluster_to_loci
    assert soft.cluster_to_loci["AAAACCCCGGGG"] == 3
    # Every emitted row's motif_key is the canonical seed
    motif_keys = {r.motif_key for r in lc.all_rows}
    assert motif_keys == {"AAAACCCCGGGG"}


def test_canonical_no_merging_when_tolerance_zero(tmp_path, synthetic_genome):
    """Same 3-locus setup but cluster_hamming_tolerance=0 leaves all 3 raw keys distinct."""
    reads = []
    for i in range(4):
        reads.append(_read_with_5p_clip(f"a{i}", 100, "TTTTAAAACCCCGGGG"))
    for i in range(2):
        reads.append(_read_with_5p_clip(f"b{i}", 200, "TTTTAAAACCCCGGGT"))
    for i in range(2):
        reads.append(_read_with_5p_clip(f"c{i}", 300, "TTTTAAAACCCCGGCG"))

    bam = _write_bam(tmp_path / "no_merge.bam", [("chr1", 1000)], reads)
    report = diagnose_bam(
        bam, synthetic_genome, min_clip_len=8, cluster_key_len=12,
        hamming_tolerance=0, cluster_hamming_tolerance=0,
    )
    soft = report.softclip
    assert soft.n_clusters == 3
    # Each locus has its own unique key with share=1
    assert all(soft.cluster_to_loci[k] == 1 for k in soft.cluster_to_loci)


def test_locus_consistent_when_outlier_fraction_below_threshold(
    tmp_path, synthetic_genome,
):
    """19 dominant + 1 outlier (5% noise) consolidates at default 0.95 threshold."""
    reads = []
    # 19 reads with the dominant key
    for i in range(19):
        reads.append(_read_with_5p_clip(f"d{i}", 100, "AAAATTTTGGGGGGGG"))
    # 1 outlier with a key Hamming-8 from dominant (not within H=1)
    reads.append(_read_with_5p_clip("o", 100, "AAAACCCCCCCCCCCC"))

    bam = _write_bam(tmp_path / "fraction_pass.bam", [("chr1", 1000)], reads)
    report = diagnose_bam(bam, synthetic_genome, min_clip_len=8, cluster_key_len=12)
    lc = report.locus_consistency
    assert lc.n_loci_consistent == 1
    assert lc.n_loci_inconsistent == 0
    # 19/20 = 95% matched → consistent at default threshold


def test_locus_inconsistent_when_outlier_fraction_above_threshold(
    tmp_path, synthetic_genome,
):
    """9 dominant + 1 outlier (10% noise) fails default 0.95 threshold."""
    reads = []
    for i in range(9):
        reads.append(_read_with_5p_clip(f"d{i}", 100, "AAAATTTTGGGGGGGG"))
    reads.append(_read_with_5p_clip("o", 100, "AAAACCCCCCCCCCCC"))

    bam = _write_bam(tmp_path / "fraction_fail.bam", [("chr1", 1000)], reads)
    report = diagnose_bam(bam, synthetic_genome, min_clip_len=8, cluster_key_len=12)
    lc = report.locus_consistency
    assert lc.n_loci_inconsistent == 1
    assert lc.n_loci_consistent == 0
    # 9/10 = 90% matched, below default 0.95 → inconsistent


def test_locus_strict_when_min_consistency_fraction_is_1(tmp_path, synthetic_genome):
    """With min_consistency_fraction=1.0 the old all-must-match rule comes back."""
    reads = []
    for i in range(19):
        reads.append(_read_with_5p_clip(f"d{i}", 100, "AAAATTTTGGGGGGGG"))
    reads.append(_read_with_5p_clip("o", 100, "AAAACCCCCCCCCCCC"))

    bam = _write_bam(tmp_path / "strict_frac.bam", [("chr1", 1000)], reads)
    report = diagnose_bam(
        bam, synthetic_genome, min_clip_len=8, cluster_key_len=12,
        min_consistency_fraction=1.0,
    )
    lc = report.locus_consistency
    assert lc.n_loci_inconsistent == 1
    assert lc.n_loci_consistent == 0


def test_split_tolerance_per_locus_tolerant_cross_locus_strict(
    tmp_path, synthetic_genome,
):
    """Per-locus de-noising on, cross-locus merging off.

    A locus with H=1 variants is consolidated (per-locus tolerant) but
    a parallel locus with the same dominant key is NOT merged with
    other near-neighbor canonical clusters (cross-locus strict).
    """
    reads = []
    # Locus A: 3 dominant + 2 H=1 variants. Per-locus check should see 1
    # consistent locus with motif "TTTTGGGGGGGG".
    for i in range(3):
        reads.append(_read_with_5p_clip(f"a{i}", 100, "AAAATTTTGGGGGGGG"))
    reads.append(_read_with_5p_clip("a3", 100, "AAAATTTTGGGGGGGT"))   # H=1
    reads.append(_read_with_5p_clip("a4", 100, "AAAATTTTGGGGGAGG"))   # H=1
    # Locus B: 2 reads with key "TTTTGGGGGGGT" (H=1 from A's dominant)
    for i in range(2):
        reads.append(_read_with_5p_clip(f"b{i}", 200, "AAAATTTTGGGGGGGT"))

    bam = _write_bam(tmp_path / "split.bam", [("chr1", 1000)], reads)
    report = diagnose_bam(
        bam, synthetic_genome, min_clip_len=8, cluster_key_len=12,
        hamming_tolerance=1, cluster_hamming_tolerance=0,
    )
    lc = report.locus_consistency
    soft = report.softclip
    # Both loci consolidate per-locus (H=1 de-noising) → consistent
    assert lc.n_loci_consistent == 2
    assert lc.n_loci_inconsistent == 0
    # But cross-locus merging is OFF, so A and B keep distinct cluster keys
    # A's dominant is "TTTTGGGGGGGG" (3 reads); B's is "TTTTGGGGGGGT" (2 reads).
    assert soft.n_clusters == 3
    assert soft.cluster_to_loci["TTTTGGGGGGGG"] == 1
    assert soft.cluster_to_loci["TTTTGGGGGGGT"] == 2  # locus A's H=1 variant + locus B


# --- adapter WARNING in the soft-clip analysis -----------------------------


def test_log_trans_splicing_verdict_flags_adapter(caplog):
    """A dominant soft-clip cluster that is Illumina adapter read-through is
    reported as adapter contamination, not as a spliced leader."""
    from eukan.assembly.align_hints import _log_trans_splicing_verdict
    report = _make_report(
        top_clusters=[("AGATCGGAAGAG", 2000, 20_000)],
        canonical_pct=99.0,
        cluster_consensus={"AGATCGGAAGAG": "AGATCGGAAGAGCACACGTCTG"},
    )
    verdict = compute_verdict(report)
    assert verdict.trans_splicing.call == "STRONG"  # adapter looks like trans-splicing
    with caplog.at_level(logging.WARNING):
        _log_trans_splicing_verdict(verdict.trans_splicing)
    text = caplog.text.lower()
    assert "adapter" in text
    assert "spliced leader" in text
    # The misleading "trans-splicing STRONG" line must not be emitted on its own.
    assert "top motif" not in text


def test_log_trans_splicing_verdict_flags_adapter_even_when_absent(caplog):
    """Adapter is surfaced even when its support is below the trans-splicing
    thresholds (call ABSENT)."""
    from eukan.assembly.align_hints import _log_trans_splicing_verdict
    report = _make_report(
        top_clusters=[("AGATCGGAAGAG", 20, 500)],  # below MODERATE
        canonical_pct=99.0,
        cluster_consensus={"AGATCGGAAGAG": "AGATCGGAAGAGCACACGTCTG"},
    )
    verdict = compute_verdict(report)
    assert verdict.trans_splicing.call == "ABSENT"
    with caplog.at_level(logging.WARNING):
        _log_trans_splicing_verdict(verdict.trans_splicing)
    assert "adapter" in caplog.text.lower()


def test_log_trans_splicing_verdict_genuine_sl_not_flagged(caplog):
    """A genuine (non-adapter) spliced-leader signal logs the trans-splicing
    warning and says nothing about adapters."""
    from eukan.assembly.align_hints import _log_trans_splicing_verdict
    report = _make_report(
        top_clusters=[("ATCGATCGATCG", 2000, 20_000)],
        canonical_pct=99.0,
        cluster_consensus={"ATCGATCGATCG": "GGGGGGATCGATCGATCGATCG"},
    )
    verdict = compute_verdict(report)
    with caplog.at_level(logging.WARNING):
        _log_trans_splicing_verdict(verdict.trans_splicing)
    text = caplog.text.lower()
    assert "adapter" not in text
    assert "trans-splicing signal strong" in text
