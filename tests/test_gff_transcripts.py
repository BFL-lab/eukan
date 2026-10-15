"""Unit tests for eukan.gff.transcripts (transcript records, reader, writer, split, stitch)."""

from __future__ import annotations

from eukan.gff.transcripts import (
    Transcript,
    iter_transcripts,
    read_transcripts,
    split_at,
    stitch,
    write_transcripts,
)
from eukan.infra.genome import ContigIndex


def test_transcript_properties():
    tx = Transcript("t", "c", "-", "s", [(1, 10), (21, 30)])
    assert (tx.start, tx.end, tx.length) == (1, 30, 20)
    assert tx.exons_5to3() == [(21, 30), (1, 10)]
    assert Transcript("u", "c", ".", "s", [(1, 10), (21, 30)]).exons_5to3() == [(1, 10), (21, 30)]


def test_iter_transcripts_gff3_and_gtf(tmp_path):
    gff = tmp_path / "a.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "c1\tsrc\tgene\t1\t90\t.\t-\t.\tID=g1\n"
        "c1\tsrc\tmRNA\t1\t90\t.\t-\t.\tID=t1;Parent=g1\n"
        "c1\tsrc\texon\t51\t90\t.\t-\t.\tID=t1.e2;Parent=t1\n"
        "c1\tsrc\texon\t1\t40\t.\t-\t.\tID=t1.e1;Parent=t1\n"
        "c1\tsrc\tmRNA\t100\t200\t.\t+\t.\tID=empty\n"
    )
    (t1,) = read_transcripts(gff)  # exon-less "empty" is dropped
    assert (t1.tid, t1.chrom, t1.strand, t1.source) == ("t1", "c1", "-", "src")
    assert t1.exons == [(1, 40), (51, 90)]

    gtf = tmp_path / "a.gtf"
    gtf.write_text(
        'c2\tStringTie\texon\t5\t9\t.\t+\t.\tgene_id "g"; transcript_id "p";\n'
        'c2\tStringTie\texon\t1\t3\t.\t+\t.\tgene_id "g"; transcript_id "p";\n'
    )
    (p,) = iter_transcripts(gtf)
    assert (p.tid, p.source, p.exons) == ("p", "StringTie", [(1, 3), (5, 9)])


def test_iter_transcripts_requires_contiguous_rows(tmp_path):
    gff = tmp_path / "a.gff3"
    gff.write_text(
        "c\ts\texon\t1\t9\t.\t+\t.\tParent=a\n"
        "c\ts\texon\t20\t29\t.\t+\t.\tParent=b\n"
        "c\ts\texon\t40\t49\t.\t+\t.\tParent=a\n"
    )
    assert [t.tid for t in iter_transcripts(gff)] == ["a", "b", "a"]


def test_write_then_read_round_trips_sorted(tmp_path):
    models = [
        Transcript("b", "c2", "+", "src", [(5, 9)]),
        Transcript("a", "c1", "-", "src", [(1, 40), (51, 90)]),
    ]
    out = tmp_path / "out.gff3"
    write_transcripts(models, out)
    lines = out.read_text().splitlines()
    assert lines[0] == "##gff-version 3"
    assert lines[1:4] == [
        "c1\tsrc\tgene\t1\t90\t.\t-\t.\tID=a.gene",
        "c1\tsrc\tmRNA\t1\t90\t.\t-\t.\tID=a;Parent=a.gene",
        "c1\tsrc\texon\t1\t40\t.\t-\t.\tID=a.exon1;Parent=a",
    ]
    assert read_transcripts(out) == sorted(models, key=lambda t: t.chrom)


def test_split_at_suffix_floor_and_strand_hook():
    tx = Transcript("t", "c", "+", "s", [(1, 40), (51, 90)])  # spliced length 80
    pieces = split_at(tx, [45, 0, 80], 25)  # out-of-range clips are ignored
    assert [(p.tid, p.exons) for p in pieces] == [
        ("t.j1", [(1, 40), (51, 55)]),
        ("t.j2", [(56, 90)]),
    ]
    assert split_at(tx, [0, 80], 25) is None
    assert [p.tid for p in split_at(tx, [10], 25)] == ["t.j1"]  # 10 bp head dropped

    spans = []

    def strand_of(lo, hi):
        spans.append((lo, hi))
        return "-" if lo > 1 else "+"

    pieces = split_at(tx, [45], 1, suffix="d", strand_of=strand_of)
    assert [(p.tid, p.strand) for p in pieces] == [("t.d1", "+"), ("t.d2", "-")]
    assert spans == [(1, 45), (46, 80)]


def test_stitch_joins_exons_and_reverse_complements(tmp_path):
    fa = tmp_path / "g.fa"
    fa.write_text(">c\nAAAACCCCGGGGTTTTACGT\n")
    with ContigIndex(fa) as contigs:
        assert stitch(Transcript("p", "c", "+", "s", [(1, 4), (9, 12)]), contigs) == "AAAAGGGG"
        assert stitch(Transcript("m", "c", "-", "s", [(1, 4), (9, 12)]), contigs) == "CCCCTTTT"
