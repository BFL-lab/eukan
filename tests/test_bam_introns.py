"""Unit tests for eukan.assembly.bam_introns (max-intron split of a BAM)."""

from __future__ import annotations

import pysam

from eukan.assembly import bam_introns
from eukan.assembly.bam_introns import _split_read, split_long_introns

_HEADER = pysam.AlignmentHeader.from_dict(
    {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 100_000}]}
)


def _aln(name, flag, start, cigar, seq):
    a = pysam.AlignedSegment(_HEADER)
    a.query_name = name
    a.flag = flag
    a.reference_id = 0
    a.reference_start = start
    a.mapping_quality = 60
    a.cigartuples = cigar
    a.query_sequence = seq
    a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    return a


# --- _split_read -----------------------------------------------------------


def test_split_one_long_intron():
    read = _aln("q", 0, 100, [(0, 10), (3, 6000), (0, 10)], "A" * 20)
    pieces = _split_read(read, _HEADER, 5000)
    assert pieces is not None and len(pieces) == 2
    assert pieces[0].reference_start == 100
    assert pieces[0].cigartuples == [(0, 10), (4, 10)]   # M then soft-clip the 3' half
    assert pieces[1].reference_start == 6110             # past the 6000-nt intron
    assert pieces[1].cigartuples == [(4, 10), (0, 10)]   # soft-clip the 5' half then M
    assert [p.query_name for p in pieces] == ["q.0", "q.1"]
    for p in pieces:
        assert p.query_sequence == "A" * 20             # full SEQ retained
        # query-consuming CIGAR lengths still sum to len(SEQ)
        q = sum(length for op, length in p.cigartuples if op in (0, 1, 4, 7, 8))
        assert q == 20


def test_split_short_intron_untouched():
    read = _aln("q", 0, 100, [(0, 10), (3, 100), (0, 10)], "A" * 20)
    assert _split_read(read, _HEADER, 5000) is None


def test_split_two_long_introns():
    read = _aln("q", 0, 100, [(0, 10), (3, 6000), (0, 10), (3, 7000), (0, 10)], "A" * 30)
    pieces = _split_read(read, _HEADER, 5000)
    assert pieces is not None and len(pieces) == 3
    assert pieces[0].reference_start == 100
    assert pieces[1].reference_start == 6110            # 100 + 10 + 6000
    assert pieces[2].reference_start == 13120           # 6110 + 10 + 7000
    assert [p.query_name for p in pieces] == ["q.0", "q.1", "q.2"]


def test_split_clears_pairing_keeps_strand():
    # paired (0x1) + reverse (0x10) + first-in-pair (0x40)
    read = _aln("q", 0x1 | 0x10 | 0x40, 100, [(0, 10), (3, 6000), (0, 10)], "A" * 20)
    pieces = _split_read(read, _HEADER, 5000)
    assert pieces is not None
    for p in pieces:
        assert not p.is_paired           # single-end now
        assert p.is_reverse              # strand preserved
        assert not p.is_secondary and not p.is_supplementary


def test_secondary_passthrough():
    read = _aln("q", 0x100, 100, [(0, 10), (3, 6000), (0, 10)], "A" * 20)
    assert _split_read(read, _HEADER, 5000) is None


def test_hard_clipped_passthrough():
    read = _aln("q", 0, 100, [(5, 4), (0, 10), (3, 6000), (0, 10)], "A" * 20)
    assert _split_read(read, _HEADER, 5000) is None


def test_disabled_when_limit_zero():
    read = _aln("q", 0, 100, [(0, 10), (3, 6000), (0, 10)], "A" * 20)
    assert _split_read(read, _HEADER, 0) is None


# --- split_long_introns (end to end) ---------------------------------------


def _sort_with_pysam(monkeypatch):
    """Patch run_cmd so the samtools sort runs via pysam's bundled htslib."""

    def fake(cmd, **kw):
        cwd = kw["cwd"]
        out_name = cmd[cmd.index("-o") + 1]
        tmp_name = cmd[-1]
        pysam.sort("-o", str(cwd / out_name), str(cwd / tmp_name))

    monkeypatch.setattr(bam_introns, "run_cmd", fake)


def test_split_long_introns_end_to_end(tmp_path, monkeypatch):
    _sort_with_pysam(monkeypatch)
    in_bam = tmp_path / "in.bam"
    with pysam.AlignmentFile(str(in_bam), "wb", header=_HEADER) as out:
        out.write(_aln("normal", 0, 50, [(0, 30)], "C" * 30))
        out.write(_aln("spliced", 0, 100, [(0, 10), (3, 6000), (0, 10)], "A" * 20))

    out_bam = tmp_path / "out.bam"
    n = split_long_introns(out_bam.parent / "in.bam", out_bam, max_intron_len=5000)
    assert n == 1
    assert not (tmp_path / "out.bam.unsorted.bam").exists()

    with pysam.AlignmentFile(str(out_bam), "rb") as bam:
        names = sorted(r.query_name for r in bam)
    assert names == ["normal", "spliced.0", "spliced.1"]
