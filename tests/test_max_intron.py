"""Unit tests for eukan.assembly.max_intron (max-intron split of transcript models)."""

from __future__ import annotations

from eukan.assembly.jaccard import _parse_transcript_models, _Tx
from eukan.assembly.max_intron import (
    _count_long_introns,
    _long_intron_cut_offsets,
    cut_models_at_max_intron,
    run_max_intron_split,
)
from eukan.assembly.sl_acceptors import AcceptorSite
from eukan.assembly.sl_cut import cut_models_at_sl
from eukan.settings import AssemblyConfig


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


# --- cut_models_at_max_intron ----------------------------------------------


def test_splits_plus(tmp_path):
    gff = tmp_path / "t.gff3"
    gff.write_text(_transcript_gff("t1", "+", [(1, 40), (5100, 5140)]))  # intron 5059 nt
    out = tmp_path / "cut.gff3"
    assert cut_models_at_max_intron(gff, out, min_segment=25, max_intron_len=5000) == (1, 1)
    models = {m.tid: m.exons for m in _parse_transcript_models(out)}
    assert models == {"t1.j1": [(1, 40)], "t1.j2": [(5100, 5140)]}


def test_short_intron_untouched(tmp_path):
    gff = tmp_path / "t.gff3"
    gff.write_text(_transcript_gff("t1", "+", [(1, 40), (51, 90)]))  # intron 10 nt
    out = tmp_path / "cut.gff3"
    assert cut_models_at_max_intron(gff, out, min_segment=25, max_intron_len=5000) == (0, 0)
    models = {m.tid: m.exons for m in _parse_transcript_models(out)}
    assert models == {"t1": [(1, 40), (51, 90)]}


def test_min_fragment_drops_tiny_tail(tmp_path):
    gff = tmp_path / "t.gff3"
    gff.write_text(_transcript_gff("t1", "+", [(1, 40), (6000, 6010)]))  # 11-nt 3' tail
    out = tmp_path / "cut.gff3"
    assert cut_models_at_max_intron(gff, out, min_segment=25, max_intron_len=5000) == (1, 1)
    models = {m.tid: m.exons for m in _parse_transcript_models(out)}
    assert models == {"t1.j1": [(1, 40)]}  # tiny tail past the long intron dropped


def test_splits_minus(tmp_path):
    gff = tmp_path / "t.gff3"
    gff.write_text(_transcript_gff("t1", "-", [(1, 40), (6000, 6040)]))
    out = tmp_path / "cut.gff3"
    assert cut_models_at_max_intron(gff, out, min_segment=25, max_intron_len=5000) == (1, 1)
    models = {m.tid: m.exons for m in _parse_transcript_models(out)}
    # 5'->3' for '-' starts at the high-coordinate exon
    assert models == {"t1.j1": [(6000, 6040)], "t1.j2": [(1, 40)]}


def test_dot_strand_split(tmp_path):
    gff = tmp_path / "t.gff3"
    gff.write_text(_transcript_gff("t1", ".", [(1, 40), (6000, 6040)]))
    out = tmp_path / "cut.gff3"
    assert cut_models_at_max_intron(gff, out, min_segment=25, max_intron_len=5000) == (1, 1)
    models = {m.tid: (m.strand, m.exons) for m in _parse_transcript_models(out)}
    assert models == {"t1.j1": (".", [(1, 40)]), "t1.j2": (".", [(6000, 6040)])}


def test_mono_exon_untouched(tmp_path):
    gff = tmp_path / "t.gff3"
    gff.write_text(_transcript_gff("t1", "+", [(1, 200)]))
    out = tmp_path / "cut.gff3"
    assert cut_models_at_max_intron(gff, out, min_segment=25, max_intron_len=5000) == (0, 0)
    models = {m.tid: m.exons for m in _parse_transcript_models(out)}
    assert models == {"t1": [(1, 200)]}


def test_idempotent(tmp_path):
    gff = tmp_path / "t.gff3"
    gff.write_text(_transcript_gff("t1", "+", [(1, 40), (5100, 5140)]))
    out1 = tmp_path / "cut1.gff3"
    cut_models_at_max_intron(gff, out1, min_segment=25, max_intron_len=5000)
    out2 = tmp_path / "cut2.gff3"
    assert cut_models_at_max_intron(out1, out2, min_segment=25, max_intron_len=5000) == (0, 0)
    a = {m.tid: m.exons for m in _parse_transcript_models(out1)}
    b = {m.tid: m.exons for m in _parse_transcript_models(out2)}
    assert a == b


def test_disabled_with_max_intron_zero(tmp_path):
    # -M 0 disables the split: the over-long model passes through unchanged...
    gff = tmp_path / "t.gff3"
    gff.write_text(_transcript_gff("t1", "+", [(1, 40), (6000, 6040)]))
    out = tmp_path / "cut.gff3"
    assert cut_models_at_max_intron(gff, out, min_segment=25, max_intron_len=0) == (0, 0)
    assert out.exists()  # ...but the file is still written (copy-through)
    models = {m.tid: m.exons for m in _parse_transcript_models(out)}
    assert models == {"t1": [(1, 40), (6000, 6040)]}


def test_passthrough_always_writes_output(tmp_path):
    # Even with nothing to cut, the sentinel output is written so the SL cut has a
    # stable input to read directly.
    gff = tmp_path / "t.gff3"
    gff.write_text(_transcript_gff("t1", "+", [(1, 40), (51, 90)]))
    out = tmp_path / "cut.gff3"
    assert cut_models_at_max_intron(gff, out, min_segment=25, max_intron_len=5000) == (0, 0)
    assert out.exists()


# --- offset primitives ------------------------------------------------------


def test_long_intron_cut_offsets_plus():
    exons = [(1, 40), (5100, 5140), (5150, 5190)]  # gap1 5059 > 5000; gap2 9 short
    assert _long_intron_cut_offsets(exons, "+", 5000) == {40}


def test_long_intron_cut_offsets_minus():
    exons_5to3 = [(5100, 5140), (1, 40)]  # descending genomic; intron 5059 nt
    assert _long_intron_cut_offsets(exons_5to3, "-", 5000) == {41}


def test_long_intron_cut_offsets_none():
    assert _long_intron_cut_offsets([(1, 40), (51, 90)], "+", 5000) == set()
    assert _long_intron_cut_offsets([(1, 40)], "+", 5000) == set()
    assert _long_intron_cut_offsets([(1, 40), (6000, 6040)], "+", 0) == set()  # disabled


def test_count_long_introns():
    tx = _Tx("t", "chr1", "+", "src", [(1, 40), (5100, 5140), (5150, 5190)])
    assert _count_long_introns(tx, 5000) == 1
    tx2 = _Tx("t", "chr1", "+", "src", [(1, 40), (51, 90)])
    assert _count_long_introns(tx2, 5000) == 0


# --- two-pass geometry equivalence -----------------------------------------


def test_two_pass_geometry_matches_union(tmp_path):
    # max_intron_split THEN sl_cut must reproduce the exon geometry the old single
    # union pass produced (the synthetic .jN ids may differ — combinr re-ids
    # everything downstream — so compare the set of exon-block lists, not ids).
    gff = tmp_path / "t.gff3"
    gff.write_text(_transcript_gff("t1", "+", [(1, 40), (51, 90), (6000, 6040)]))
    sites = [AcceptorSite("chr1", 60, "+", 5, ("reads",))]
    maxi = tmp_path / "t.maxintron.gff3"
    cut_models_at_max_intron(gff, maxi, min_segment=25, max_intron_len=5000)
    cut = tmp_path / "t.cut.gff3"
    cut_models_at_sl(maxi, sites, cut, min_segment=25)
    geometry = {tuple(m.exons) for m in _parse_transcript_models(cut)}
    assert geometry == {
        ((1, 40), (51, 59)),
        ((60, 90),),
        ((6000, 6040),),
    }


# --- run_max_intron_split ---------------------------------------------------


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


def test_run_prefers_defuse_over_stranded_over_raw(tmp_path):
    # Per track, max_intron_split picks the most-processed model variant via
    # tracks.resolve_model_source: {stem}.defuse > {stem}.stranded > {stem}.gff3, and
    # writes {stem}.maxintron.gff3 (copy-through here, since the models are mono-exon).
    denovo, gg = "trinity-denovo.genome", "trinity-gg.genome"
    (tmp_path / f"{denovo}.gff3").write_text(_track_models_gff3("raw"))
    (tmp_path / f"{denovo}.stranded.gff3").write_text(_track_models_gff3("stranded"))
    (tmp_path / f"{denovo}.defuse.gff3").write_text(_track_models_gff3("defuse"))
    (tmp_path / f"{gg}.gff3").write_text(_track_models_gff3("raw"))
    (tmp_path / f"{gg}.stranded.gff3").write_text(_track_models_gff3("stranded"))

    config = AssemblyConfig(genome=tmp_path / "g.fa", work_dir=tmp_path, num_cpu=1)
    run_max_intron_split(config)

    denovo_out = (tmp_path / f"{denovo}.maxintron.gff3").read_text()
    gg_out = (tmp_path / f"{gg}.maxintron.gff3").read_text()
    assert _source_of(denovo_out) == "defuse"  # defuse beats stranded + raw
    assert _source_of(gg_out) == "stranded"  # stranded beats raw (no defuse)
    for out in (denovo_out, gg_out):
        assert out.count("\tmRNA\t") == 2
        assert "ID=A;" in out and "ID=B;" in out


def test_run_skips_track_with_no_models(tmp_path):
    # A Trinity mode that produced nothing leaves resolve_model_source -> None;
    # that track is skipped and no {stem}.maxintron.gff3 is written for it.
    denovo, gg = "trinity-denovo.genome", "trinity-gg.genome"
    (tmp_path / f"{denovo}.gff3").write_text(_track_models_gff3("raw"))

    config = AssemblyConfig(genome=tmp_path / "g.fa", work_dir=tmp_path, num_cpu=1)
    run_max_intron_split(config)

    assert (tmp_path / f"{denovo}.maxintron.gff3").exists()
    assert not (tmp_path / f"{gg}.maxintron.gff3").exists()
