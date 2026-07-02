"""Unit tests for eukan.assembly.sl_cut (genomic SL cut of transcript models).

The over-long-intron split moved to :mod:`eukan.assembly.max_intron`
(see ``tests/test_max_intron.py``); this module now covers the SL cut only.
"""

from __future__ import annotations

import pysam
import pytest

from eukan.assembly.jaccard import _parse_transcript_models
from eukan.assembly.sl_acceptors import AcceptorSite
from eukan.assembly.sl_cut import (
    _project_genomic_to_spliced,
    bam_to_transcript_gff3,
    cut_models_at_sl,
    run_sl_cut,
)
from eukan.settings import AssemblyConfig


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


def _transcript_gff(tid, strand, exons):
    lo, hi = exons[0][0], exons[-1][1]
    rows = [
        "##gff-version 3",
        f"chr1\ttest\tgene\t{lo}\t{hi}\t.\t{strand}\t.\tID={tid}.g",
        f"chr1\ttest\tmRNA\t{lo}\t{hi}\t.\t{strand}\t.\tID={tid};Parent={tid}.g",
    ]
    for k, (s, e) in enumerate(exons, 1):
        rows.append(f"chr1\ttest\texon\t{s}\t{e}\t.\t{strand}\t.\tID={tid}.e{k};Parent={tid}")
    return "\n".join(rows) + "\n"


# --- _project_genomic_to_spliced -------------------------------------------


def test_project_plus_strand():
    exons = [(1, 40), (51, 90)]  # 5'->3' for '+'
    assert _project_genomic_to_spliced(exons, "+", 60) == 50
    assert _project_genomic_to_spliced(exons, "+", 1) == 1
    assert _project_genomic_to_spliced(exons, "+", 45) is None  # intronic (41-50)


def test_project_minus_strand():
    exons_5to3 = [(51, 90), (1, 40)]  # descending genomic for '-'
    assert _project_genomic_to_spliced(exons_5to3, "-", 90) == 1  # 5'-most base
    assert _project_genomic_to_spliced(exons_5to3, "-", 60) == 31


# --- bam_to_transcript_gff3 ------------------------------------------------


def test_bam_to_gff3_splits_at_intron(tmp_path):
    bam = tmp_path / "tx.bam"
    _make_bam(bam, [("q1", 0, 100, [(0, 10), (3, 50), (0, 10)], "A" * 20)])
    out = tmp_path / "tx.gff3"
    assert bam_to_transcript_gff3(bam, out, "trinity-denovo.genome") == 1
    (m,) = _parse_transcript_models(out)
    assert m.exons == [(101, 110), (161, 170)]
    assert m.strand == "+"


def test_bam_to_gff3_minus_strand(tmp_path):
    bam = tmp_path / "tx.bam"
    _make_bam(bam, [("q1", 16, 100, [(0, 10), (3, 50), (0, 10)], "A" * 20)])
    out = tmp_path / "tx.gff3"
    bam_to_transcript_gff3(bam, out, "trinity-denovo.genome")
    (m,) = _parse_transcript_models(out)
    assert m.strand == "-" and m.exons == [(101, 110), (161, 170)]


def test_bam_to_gff3_unique_ids_for_multimaps(tmp_path):
    bam = tmp_path / "tx.bam"
    _make_bam(bam, [
        ("q1", 0, 100, [(0, 20)], "A" * 20),
        ("q1", 256, 500, [(0, 20)], "A" * 20),  # secondary: a second locus
    ])
    out = tmp_path / "tx.gff3"
    assert bam_to_transcript_gff3(bam, out, "trinity-gg.genome") == 2
    ids = {m.tid for m in _parse_transcript_models(out)}
    assert ids == {"q1.m1", "q1.m2"}  # no Parent collision → loci stay distinct


# --- cut_models_at_sl ------------------------------------------------------


@pytest.mark.parametrize("strand,expected", [
    ("+", {"t1.j1": [(1, 40), (51, 59)], "t1.j2": [(60, 90)]}),
    ("-", {"t1.j1": [(61, 90)], "t1.j2": [(1, 40), (51, 60)]}),
])
def test_cut_both_strands(tmp_path, strand, expected):
    gff = tmp_path / "t.gff3"
    gff.write_text(_transcript_gff("t1", strand, [(1, 40), (51, 90)]))
    sites = [AcceptorSite("chr1", 60, strand, 5, ("reads",))]
    out = tmp_path / "cut.gff3"
    assert cut_models_at_sl(gff, sites, out, min_segment=25) == 1
    models = {m.tid: m.exons for m in _parse_transcript_models(out)}
    assert models == expected


def test_cut_dot_strand_oriented_by_acceptor(tmp_path):
    gff = tmp_path / "t.gff3"
    gff.write_text(_transcript_gff("t1", ".", [(1, 40), (51, 90)]))
    sites = [AcceptorSite("chr1", 60, "+", 5, ("reads",))]
    out = tmp_path / "cut.gff3"
    assert cut_models_at_sl(gff, sites, out, min_segment=25) == 1
    models = {m.tid: (m.strand, m.exons) for m in _parse_transcript_models(out)}
    assert all(strand == "+" for strand, _ in models.values())  # SL imposed strand
    assert models["t1.j2"][1] == [(60, 90)]


def test_cut_ignores_opposite_strand_site(tmp_path):
    gff = tmp_path / "t.gff3"
    gff.write_text(_transcript_gff("t1", "+", [(1, 40), (51, 90)]))
    sites = [AcceptorSite("chr1", 60, "-", 5, ("reads",))]  # wrong strand
    out = tmp_path / "cut.gff3"
    assert cut_models_at_sl(gff, sites, out, min_segment=25) == 0
    models = {m.tid: m.exons for m in _parse_transcript_models(out)}
    assert models == {"t1": [(1, 40), (51, 90)]}


def test_cut_passthrough_without_sites(tmp_path):
    gff = tmp_path / "t.gff3"
    gff.write_text(_transcript_gff("t1", "+", [(1, 40), (51, 90)]))
    out = tmp_path / "cut.gff3"
    assert cut_models_at_sl(gff, [], out, min_segment=25) == 0
    models = {m.tid: m.exons for m in _parse_transcript_models(out)}
    assert models == {"t1": [(1, 40), (51, 90)]}


# --- run_sl_cut ------------------------------------------------------------


def _track_models_gff3(tag):
    """Two-model gene>mRNA>exon GFF3 whose source column encodes the variant *tag*."""
    return (
        "##gff-version 3\n"
        f"chr1\t{tag}\tgene\t1\t40\t.\t+\t.\tID=A.gene\n"
        f"chr1\t{tag}\tmRNA\t1\t40\t.\t+\t.\tID=A;Parent=A.gene\n"
        f"chr1\t{tag}\texon\t1\t40\t.\t+\t.\tID=A.e1;Parent=A\n"
        f"chr1\t{tag}\tgene\t60\t100\t.\t+\t.\tID=B.gene\n"
        f"chr1\t{tag}\tmRNA\t60\t100\t.\t+\t.\tID=B;Parent=B.gene\n"
        f"chr1\t{tag}\texon\t60\t100\t.\t+\t.\tID=B.e1;Parent=B\n"
    )


def _source_of(gff3_text):
    """The source (column 2) of the first exon row — identifies which variant was read."""
    for line in gff3_text.splitlines():
        cols = line.split("\t")
        if len(cols) >= 9 and cols[2] == "exon":
            return cols[1]
    return None


def test_run_sl_cut_reads_maxintron_models(tmp_path):
    # run_sl_cut reads each track's max-intron-split models ({stem}.maxintron.gff3,
    # always written by the max_intron_split step) and writes {stem}.cut.gff3.
    denovo, gg = "trinity-denovo.genome", "trinity-gg.genome"
    (tmp_path / f"{denovo}.maxintron.gff3").write_text(_track_models_gff3("maxintron"))
    (tmp_path / f"{gg}.maxintron.gff3").write_text(_track_models_gff3("maxintron"))

    config = AssemblyConfig(genome=tmp_path / "g.fa", work_dir=tmp_path, num_cpu=1)
    run_sl_cut(config)  # no SL acceptors -> passthrough copy

    for stem in (denovo, gg):
        out = (tmp_path / f"{stem}.cut.gff3").read_text()
        assert _source_of(out) == "maxintron"
        assert out.count("\tmRNA\t") == 2
        assert "ID=A;" in out and "ID=B;" in out


def test_run_sl_cut_skips_track_with_no_maxintron(tmp_path):
    # A track without a {stem}.maxintron.gff3 (e.g. that Trinity mode produced no
    # models, so max_intron_split skipped it) is skipped here too.
    denovo, gg = "trinity-denovo.genome", "trinity-gg.genome"
    (tmp_path / f"{denovo}.maxintron.gff3").write_text(_track_models_gff3("maxintron"))

    config = AssemblyConfig(genome=tmp_path / "g.fa", work_dir=tmp_path, num_cpu=1)
    run_sl_cut(config)

    assert (tmp_path / f"{denovo}.cut.gff3").exists()
    assert not (tmp_path / f"{gg}.cut.gff3").exists()
