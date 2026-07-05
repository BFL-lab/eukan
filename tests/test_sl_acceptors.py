"""Unit tests for eukan.assembly.sl_acceptors (SL trans-splice acceptor detection)."""

from __future__ import annotations

import json

import pysam

from eukan.assembly.sl_acceptors import (
    AcceptorSite,
    _iter_sl_ops,
    build_joint_consensus,
    detect_sl_acceptors,
    load_sl_acceptors,
)
from eukan.assembly.sl_depletion import _revcomp, is_adapter
from eukan.infra.artifacts import Artifact
from eukan.settings import AssemblyConfig

SL = "GTACTTTATT"            # 10 bp, meets _MIN_MOTIF_LEN
SL16 = "GTACTTTATTCCGGAA"    # 16 bp, meets _CONSENSUS_LEN
ADAPT = "AGATCGGAAGAGC"      # Illumina universal / TruSeq read-through seed
ADAPT16 = "AGATCGGAAGAGCACA"  # 16 bp read-through window (PE1/PE2 rc start; holds the seed)


def _header(ref="chr1", ref_len=100_000):
    return pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}, "SQ": [{"SN": ref, "LN": ref_len}]})


def _seg(flag, start, cigar, seq, header=None):
    s = pysam.AlignedSegment(header or _header())
    s.query_name = "r"
    s.flag = flag
    s.reference_id = 0
    s.reference_start = start
    s.mapping_quality = 60
    s.cigartuples = cigar
    s.query_sequence = seq
    return s


def _make_bam(path, reads, ref="chr1", ref_len=100_000):
    """reads = [(name, flag, start, cigar, seq), ...]."""
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": ref, "LN": ref_len}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for name, flag, start, cigar, seq in reads:
            s = pysam.AlignedSegment(out.header)
            s.query_name = name
            s.flag = flag
            s.reference_id = 0
            s.reference_start = start
            s.mapping_quality = 60
            s.cigartuples = cigar
            s.query_sequence = seq
            out.write(s)


def _config(tmp_path, **kw):
    return AssemblyConfig(
        genome=tmp_path / "g.fa", work_dir=tmp_path, num_cpu=1, **kw
    )


# --- _iter_sl_ops geometry -------------------------------------------------


def test_iter_sl_ops_forward_5p_clip():
    seg = _seg(0, 100, [(4, 10), (0, 20)], SL + "A" * 20)
    ops = list(_iter_sl_ops(seg, min_clip_len=8, min_ins_len=10, scan_insertions=False))
    assert ops == [(101, "+", SL)]  # acceptor = reference_start + 1


def test_iter_sl_ops_reverse_5p_clip():
    seg = _seg(16, 100, [(0, 20), (4, 10)], "A" * 20 + SL)
    ops = list(_iter_sl_ops(seg, min_clip_len=8, min_ins_len=10, scan_insertions=False))
    assert ops == [(120, "-", SL)]  # acceptor = reference_end (1-based last aligned base)


def test_iter_sl_ops_internal_insertion():
    seg = _seg(0, 100, [(0, 10), (1, 10), (0, 10)], "A" * 10 + SL + "A" * 10)
    ops = list(_iter_sl_ops(seg, min_clip_len=8, min_ins_len=10, scan_insertions=True))
    assert ops == [(111, "+", SL)]  # acceptor = first aligned base after the insertion


def test_iter_sl_ops_terminal_insertion_ignored():
    # A leading insertion is not internal → never an acceptor.
    seg = _seg(0, 100, [(1, 10), (0, 20)], SL + "A" * 20)
    ops = list(_iter_sl_ops(seg, min_clip_len=8, min_ins_len=10, scan_insertions=True))
    assert ops == []


def test_iter_sl_ops_reverse_read_leading_clip():
    # NODE_574: an antisense-assembled contig maps reverse yet carries its forward
    # SL in the *leading* clip. Geometry (not is_reverse) makes it a '+' acceptor.
    seg = _seg(16, 100, [(4, 16), (0, 20)], SL16 + "A" * 20)
    ops = list(_iter_sl_ops(seg, min_clip_len=8, min_ins_len=10, scan_insertions=False))
    assert ops == [(101, "+", SL16)]


def test_iter_sl_ops_forward_read_trailing_clip():
    # A forward read's trailing clip is now also examined → a '-' acceptor.
    seg = _seg(0, 100, [(0, 20), (4, 16)], "A" * 20 + SL16)
    ops = list(_iter_sl_ops(seg, min_clip_len=8, min_ins_len=10, scan_insertions=False))
    assert ops == [(120, "-", SL16)]  # acceptor = reference_end


def test_iter_sl_ops_both_terminal_clips():
    seg = _seg(0, 100, [(4, 16), (0, 20), (4, 16)], SL16 + "A" * 20 + SL16)
    ops = list(_iter_sl_ops(seg, min_clip_len=8, min_ins_len=10, scan_insertions=False))
    assert ops == [(101, "+", SL16), (120, "-", SL16)]


# --- build_joint_consensus -------------------------------------------------


def _write_verdict(tmp_path, call, consensus=SL):
    summary = {"verdict": {"trans_splicing": {
        "call": call, "top_non_trivial_cluster_consensus": consensus,
        "top_non_trivial_cluster_key": consensus,
    }}}
    (tmp_path / Artifact.SOFTCLIP_DIAGNOSTIC.value).write_text(json.dumps(summary))


def test_consensus_override_wins(tmp_path):
    _write_verdict(tmp_path, "ABSENT")
    cfg = _config(tmp_path, sl_sequence="acgtacgtac")
    assert build_joint_consensus(cfg, []) == "ACGTACGTAC"


def test_consensus_from_strong_verdict(tmp_path):
    _write_verdict(tmp_path, "STRONG")
    assert build_joint_consensus(_config(tmp_path), []) == SL


def test_consensus_from_de_novo_fallback(tmp_path):
    # No override, no verdict file → fall back to the dominant de novo insertion.
    bam = tmp_path / "trinity-denovo.genome.bam"
    _make_bam(bam, [
        ("t1", 0, 100, [(0, 10), (1, 16), (0, 10)], "A" * 10 + SL16 + "A" * 10),
        ("t2", 0, 200, [(0, 10), (1, 16), (0, 10)], "A" * 10 + SL16 + "A" * 10),
    ])
    assert build_joint_consensus(_config(tmp_path), [bam]) == SL16


def test_consensus_none_without_signal(tmp_path):
    assert build_joint_consensus(_config(tmp_path), []) is None


# --- detect_sl_acceptors end-to-end ----------------------------------------


def test_detect_pools_reads_and_de_novo(tmp_path):
    cfg = _config(tmp_path, sl_sequence=SL)
    _make_bam(tmp_path / "minimap2_Aligned.sortedByCoord.out.bam",
              [("r1", 0, 100, [(4, 10), (0, 30)], SL + "A" * 30)])           # acceptor 101 (+)
    _make_bam(tmp_path / "trinity-denovo.genome.bam",
              [("t1", 0, 100, [(0, 10), (1, 10), (0, 10)], "A" * 10 + SL + "A" * 10)])  # 111 (+)

    detect_sl_acceptors(cfg)

    sites = {(s.pos, s.strand): s for s in load_sl_acceptors(tmp_path / "sl_acceptors.gff3")}
    assert (101, "+") in sites and sites[(101, "+")].sources == ("reads",)
    assert (111, "+") in sites and sites[(111, "+")].sources == ("trinity-denovo",)


def test_detect_clusters_and_unions_sources(tmp_path):
    cfg = _config(tmp_path, sl_sequence=SL, sl_cluster_window=5)
    # Two acceptors 2 bp apart (within the window) from different sources → one site.
    _make_bam(tmp_path / "minimap2_Aligned.sortedByCoord.out.bam",
              [("r1", 0, 100, [(4, 10), (0, 30)], SL + "A" * 30)])                       # 101
    _make_bam(tmp_path / "trinity-denovo.genome.bam",
              [("t1", 0, 102, [(4, 10), (0, 30)], SL + "A" * 30)])                       # 103
    detect_sl_acceptors(cfg)

    sites = load_sl_acceptors(tmp_path / "sl_acceptors.gff3")
    plus = [s for s in sites if s.strand == "+"]
    assert len(plus) == 1
    assert plus[0].support == 2
    assert set(plus[0].sources) == {"reads", "trinity-denovo"}


def test_detect_recovers_antisense_leading_clip(tmp_path):
    """Fix A: a reverse-mapped de novo contig with its forward SL in the leading
    clip (NODE_574) now yields a '+' acceptor — previously dropped because only the
    trailing clip was inspected for reverse reads."""
    cfg = _config(tmp_path, sl_sequence=SL16)
    _make_bam(tmp_path / "trinity-denovo.genome.bam",
              [("t1", 16, 100, [(4, 16), (0, 20)], SL16 + "A" * 20)])  # reverse, leading SL
    detect_sl_acceptors(cfg)

    sites = {(s.pos, s.strand): s for s in load_sl_acceptors(tmp_path / "sl_acceptors.gff3")}
    assert (101, "+") in sites and sites[(101, "+")].sources == ("trinity-denovo",)


def test_detect_matches_core_of_overlong_consensus(tmp_path):
    """Fix B: a 25 nt read-verdict consensus is matched by its conserved 16 nt 3'
    core, so a 17 nt captured leader carrying only the core is still detected."""
    cfg = _config(tmp_path, sl_sequence="AAAGCTACAGTTTCTGTACTTTATT")  # core = GTTTCTGTACTTTATT
    _make_bam(tmp_path / "trinity-denovo.genome.bam",
              [("t1", 0, 100, [(4, 17), (0, 30)], "AGTTTCTGTACTTTATT" + "A" * 30)])
    detect_sl_acceptors(cfg)

    sites = {(s.pos, s.strand) for s in load_sl_acceptors(tmp_path / "sl_acceptors.gff3")}
    assert (101, "+") in sites


def test_detect_is_orientation_aware(tmp_path):
    """A clip whose motif orientation contradicts its geometry strand is rejected;
    a genuinely reverse-complement leader at a trailing clip is kept as '-'."""
    cfg = _config(tmp_path, sl_sequence=SL16)
    _make_bam(tmp_path / "trinity-denovo.genome.bam", [
        ("wrong", 0, 100, [(0, 20), (4, 16)], "A" * 20 + SL16),            # fwd SL, '-' geom → drop
        ("right", 0, 200, [(0, 20), (4, 16)], "A" * 20 + _revcomp(SL16)),  # RC SL, '-' geom → keep
    ])
    detect_sl_acceptors(cfg)

    sites = {(s.pos, s.strand) for s in load_sl_acceptors(tmp_path / "sl_acceptors.gff3")}
    assert (220, "-") in sites      # reference_end of the RC-leader read
    assert (120, "-") not in sites  # forward SL at a trailing clip is not a '-' leader


def test_detect_noop_without_signal(tmp_path):
    cfg = _config(tmp_path)  # no override, no verdict, no de novo signal
    _make_bam(tmp_path / "minimap2_Aligned.sortedByCoord.out.bam",
              [("r1", 0, 100, [(0, 30)], "A" * 30)])  # clean, no clip
    detect_sl_acceptors(cfg)

    out = tmp_path / "sl_acceptors.gff3"
    assert out.exists()
    assert load_sl_acceptors(out) == []


def test_acceptor_site_roundtrip(tmp_path):
    out = tmp_path / "acc.gff3"
    from eukan.assembly.sl_acceptors import _write_acceptors
    sites = [
        AcceptorSite("chr1", 500, "+", 12, ("reads", "trinity-denovo")),
        AcceptorSite("chr2", 30, "-", 3, ("trinity-denovo",)),
    ]
    _write_acceptors(sites, out)
    assert load_sl_acceptors(out) == sites


# --- is_adapter (adapter-vs-SL discrimination) -----------------------------


def test_is_adapter_illumina_read_through():
    assert is_adapter("AGATCGGAAGAGCACACGTCTGAAC")  # full Illumina read-through


def test_is_adapter_reverse_complement():
    # Adapter on the opposite strand (RC of a read-through window) is still caught.
    assert is_adapter(_revcomp(ADAPT16))


def test_is_adapter_eleven_bp_window_flagged():
    # 11 contiguous bp of the seed embedded in a longer clip → adapter.
    assert is_adapter("GGGG" + ADAPT[:11] + "TTTT")


def test_is_adapter_ten_bp_share_not_flagged():
    # Shares only 10 bp with the seed (breaks at bp 11) → below the 11 bp floor.
    assert not is_adapter("TTTT" + ADAPT[:10] + "TTTT")


def test_is_adapter_nextera_flagged():
    assert is_adapter("GGG" + "CTGTCTCTTATACACATCT" + "GGG")


def test_is_adapter_real_sl_not_flagged():
    assert not is_adapter(SL16)


def test_is_adapter_empty_list_disables():
    assert not is_adapter(ADAPT16, [])


# --- adapter discrimination in consensus selection -------------------------


def _denovo_ins_reads(n, motif16, prefix):
    """n de novo reads each carrying *motif16* as an internal 16 bp insertion."""
    return [
        (f"{prefix}{i}", 0, 100 + i * 50, [(0, 10), (1, 16), (0, 10)],
         "C" * 10 + motif16 + "C" * 10)
        for i in range(n)
    ]


def test_consensus_adapter_only_denovo_is_none(tmp_path):
    bam = tmp_path / "trinity-denovo.genome.bam"
    _make_bam(bam, _denovo_ins_reads(3, ADAPT16, "a"))
    assert build_joint_consensus(_config(tmp_path), [bam]) is None


def test_consensus_real_sl_wins_over_dominant_adapter(tmp_path):
    # Adapter is the most common window (3 reads) but a real SL ranks behind it
    # (2 reads, still >= _MIN_DENOVO_SUPPORT once head+tail are pooled). The
    # skip-and-continue must surface the SL rather than lock onto the adapter.
    bam = tmp_path / "trinity-denovo.genome.bam"
    _make_bam(bam, _denovo_ins_reads(3, ADAPT16, "a") + _denovo_ins_reads(2, SL16, "s"))
    assert build_joint_consensus(_config(tmp_path), [bam]) == SL16


def test_consensus_adapter_verdict_rejected(tmp_path):
    _write_verdict(tmp_path, "STRONG", consensus=ADAPT16)
    assert build_joint_consensus(_config(tmp_path), []) is None


def test_consensus_adapter_override_not_filtered(tmp_path):
    # An explicit --sl-sequence override is trusted as-is, even if adapter-like.
    assert build_joint_consensus(_config(tmp_path, sl_sequence=ADAPT16), []) == ADAPT16


def test_consensus_adapter_filter_off_lets_adapter_through(tmp_path):
    bam = tmp_path / "trinity-denovo.genome.bam"
    _make_bam(bam, _denovo_ins_reads(3, ADAPT16, "a"))
    cfg = _config(tmp_path, sl_adapter_filter=False)
    assert build_joint_consensus(cfg, [bam]) == ADAPT16


def test_detect_adapter_only_is_noop(tmp_path):
    """S. pombe in miniature: the only soft-clip / insertion signal is Illumina
    adapter read-through, so no SL consensus is recovered and no acceptors written."""
    cfg = _config(tmp_path)
    _make_bam(tmp_path / "minimap2_Aligned.sortedByCoord.out.bam",
              [("r1", 0, 100, [(0, 30), (4, 16)], "A" * 30 + ADAPT16)])  # trailing adapter clip
    _make_bam(tmp_path / "trinity-denovo.genome.bam", _denovo_ins_reads(3, ADAPT16, "a"))
    detect_sl_acceptors(cfg)
    out = tmp_path / "sl_acceptors.gff3"
    assert out.exists()
    assert load_sl_acceptors(out) == []


def test_detect_scan_skips_adapter_clip_under_adapter_override(tmp_path):
    # Even when an adapter is the (override-trusted) consensus, the defensive scan
    # skips adapter ops, so adapter read-through seeds no acceptor sites.
    cfg = _config(tmp_path, sl_sequence=ADAPT16)
    _make_bam(tmp_path / "trinity-denovo.genome.bam",
              [("t1", 0, 100, [(4, 16), (0, 20)], ADAPT16 + "A" * 20)])  # leading adapter clip
    detect_sl_acceptors(cfg)
    assert load_sl_acceptors(tmp_path / "sl_acceptors.gff3") == []


def test_detect_scan_filter_off_keeps_adapter_clip(tmp_path):
    # With the filter off the same adapter clip is matched against the adapter
    # consensus and recorded — confirming the toggle gates the scan, not just
    # consensus selection.
    cfg = _config(tmp_path, sl_sequence=ADAPT16, sl_adapter_filter=False)
    _make_bam(tmp_path / "trinity-denovo.genome.bam",
              [("t1", 0, 100, [(4, 16), (0, 20)], ADAPT16 + "A" * 20)])
    detect_sl_acceptors(cfg)
    sites = {(s.pos, s.strand) for s in load_sl_acceptors(tmp_path / "sl_acceptors.gff3")}
    assert (101, "+") in sites
