"""Unit tests for eukan.assembly.defuse (homology-grounded transcript de-fusion)."""

from __future__ import annotations

from eukan.assembly import defuse
from eukan.assembly.defuse import (
    _Hit,
    _overlap_frac,
    distinct_nonoverlapping,
    parse_positional_hits,
    run_defuse,
    split_fused,
)
from eukan.assembly.jaccard import _Tx
from eukan.settings import AssemblyConfig


def _hit(sseqid, lo, hi, bitscore=100.0, qframe=1):
    return _Hit(sseqid, lo, hi, bitscore, qframe)


# --- hit parsing -----------------------------------------------------------


def test_parse_positional_hits_normalizes_and_groups(tmp_path):
    tsv = tmp_path / "h.tsv"
    tsv.write_text(
        "trinity-denovo.genome:t1\tP1\t10\t80\t120.0\t1\n"
        "trinity-denovo.genome:t1\tP2\t300\t120\t90.0\t-2\n"  # qstart > qend (minus frame)
        "trinity-gg.genome:t2\tP3\t1\t50\t60.0\t1\n"
        "# comment\n"
    )
    hits = parse_positional_hits(tsv)
    assert set(hits) == {"trinity-denovo.genome:t1", "trinity-gg.genome:t2"}
    a, b = hits["trinity-denovo.genome:t1"]
    assert (a.sseqid, a.qlo, a.qhi, a.qframe) == ("P1", 10, 80, 1)
    assert (b.sseqid, b.qlo, b.qhi, b.qframe) == ("P2", 120, 300, -2)  # normalized


def test_overlap_frac():
    # 21 bp overlap of a 31-bp shorter hit
    assert _overlap_frac(_hit("A", 100, 200), _hit("B", 180, 210)) == 21 / 31
    assert _overlap_frac(_hit("A", 1, 100), _hit("B", 200, 300)) == 0.0


# --- distinct, non-overlapping hit selection -------------------------------


def test_distinct_nonoverlapping_two_genes():
    hits = [_hit("P1", 1, 80), _hit("P2", 120, 200)]
    chosen = distinct_nonoverlapping(hits, 0.10)
    assert chosen is not None
    assert [h.sseqid for h in chosen] == ["P1", "P2"]  # ordered by qlo


def test_distinct_nonoverlapping_same_subject_is_not_a_fusion():
    # Same protein hitting two places (repeat/paralog domain) is one gene -> no split.
    hits = [_hit("P1", 1, 80), _hit("P1", 120, 200)]
    assert distinct_nonoverlapping(hits, 0.10) is None


def test_distinct_nonoverlapping_overlapping_hits_rejected():
    # Two distinct subjects but their query ranges overlap beyond tolerance.
    hits = [_hit("P1", 1, 100), _hit("P2", 50, 150)]  # 51/100 overlap
    assert distinct_nonoverlapping(hits, 0.10) is None


def test_distinct_nonoverlapping_within_tolerance_kept():
    # 6 bp overlap of a 95-bp shorter hit ~ 6.3% <= 10% tolerance -> still distinct.
    hits = [_hit("P1", 1, 100), _hit("P2", 95, 200)]
    assert distinct_nonoverlapping(hits, 0.10) is not None
    assert distinct_nonoverlapping(hits, 0.05) is None  # tighter tolerance rejects


def test_distinct_nonoverlapping_keeps_best_per_subject():
    hits = [_hit("P1", 1, 80, bitscore=50.0), _hit("P1", 1, 80, bitscore=200.0)]
    assert distinct_nonoverlapping(hits, 0.10) is None  # collapses to one subject


# --- the split itself ------------------------------------------------------


def test_split_fused_two_exon_genes_at_boundary():
    tx = _Tx("t1", "c", "+", "src", [(1, 100), (201, 300)])  # spliced len 200
    chosen = [_hit("P1", 1, 80), _hit("P2", 120, 200)]  # gap midpoint = 100
    pieces = split_fused(tx, chosen, min_segment=25)
    assert pieces is not None and len(pieces) == 2
    assert pieces[0].exons == [(1, 100)] and pieces[1].exons == [(201, 300)]
    assert [p.tid for p in pieces] == ["t1.d1", "t1.d2"]
    assert all(p.strand == "+" for p in pieces)


def test_split_fused_opposite_strand_pieces():
    tx = _Tx("t1", "c", "+", "src", [(1, 100), (201, 300)])
    chosen = [_hit("P1", 1, 80, qframe=1), _hit("P2", 120, 200, qframe=-1)]
    pieces = split_fused(tx, chosen, min_segment=25)
    assert pieces is not None
    assert pieces[0].strand == "+" and pieces[1].strand == "-"


def test_split_fused_three_genes():
    tx = _Tx("t", "c", "+", "s", [(1, 100), (201, 300), (401, 500)])  # spliced 300
    chosen = [_hit("P1", 1, 80), _hit("P2", 120, 180), _hit("P3", 220, 300)]
    pieces = split_fused(tx, chosen, min_segment=25)
    assert pieces is not None and len(pieces) == 3
    assert [p.exons for p in pieces] == [[(1, 100)], [(201, 300)], [(401, 500)]]


def test_split_fused_minus_strand():
    tx = _Tx("m", "c", "-", "s", [(1, 100), (201, 300)])
    chosen = [_hit("P1", 1, 80), _hit("P2", 120, 200)]  # query is 5'->3' (genomic 3'->5')
    pieces = split_fused(tx, chosen, min_segment=25)
    assert pieces is not None and len(pieces) == 2
    assert {tuple(p.exons) for p in pieces} == {((1, 100),), ((201, 300),)}


def test_split_fused_drops_to_single_piece_is_no_split():
    tx = _Tx("t", "c", "+", "s", [(1, 200)])  # single exon, spliced 200
    chosen = [_hit("P1", 1, 5), _hit("P2", 40, 200)]  # midpoint = 22
    # piece1 (spliced 1..22 = 22 bp) is below min_segment -> only one survivor -> no split
    assert split_fused(tx, chosen, min_segment=50) is None


# --- step command construction ---------------------------------------------


def _config(tmp_path, **kw):
    (tmp_path / "genome.fa").write_text(">c\n" + "ACGT" * 100 + "\n")  # 400 bp
    return AssemblyConfig(
        genome=tmp_path / "genome.fa", work_dir=tmp_path, num_cpu=2,
        defuse=True, uniprot_db=tmp_path / "db.dmnd", **kw,
    )


def _exon_gff3(tid: str) -> str:
    return (
        f"##gff-version 3\n"
        f"c\tTrinity\texon\t1\t100\t.\t+\t.\tParent={tid}\n"
        f"c\tTrinity\texon\t201\t300\t.\t+\t.\tParent={tid}\n"
    )


def test_run_defuse_diamond_command_and_copy_through(tmp_path, monkeypatch):
    # Two mapped Trinity tracks, each with a raw genome GFF3 (no .stranded.gff3).
    (tmp_path / "trinity-denovo.genome.gff3").write_text(_exon_gff3("t1"))
    (tmp_path / "trinity-gg.genome.gff3").write_text(_exon_gff3("t2"))
    cmds: list[list[str]] = []
    monkeypatch.setattr(defuse, "run_cmd", lambda cmd, **kw: cmds.append(cmd))
    monkeypatch.setattr(defuse, "_resolve_diamond_db", lambda config: "db")
    monkeypatch.setattr(defuse, "parse_positional_hits", lambda path: {})  # no hits
    monkeypatch.setattr(defuse, "_open_indexed_bam", lambda path: None)

    run_defuse(_config(tmp_path))

    (cmd,) = cmds
    assert cmd[:2] == ["diamond", "blastx"]
    assert "--ultra-sensitive" in cmd
    # positional output columns are requested (needed to locate the split point)
    assert "qstart" in cmd and "qend" in cmd
    assert cmd[cmd.index("--outfmt") + 1] == "6"
    assert cmd[cmd.index("--db") + 1] == "db"
    # no fusion -> copy-through of each track's input model
    out_denovo = tmp_path / "trinity-denovo.genome.defuse.gff3"
    out_gg = tmp_path / "trinity-gg.genome.defuse.gff3"
    assert out_denovo.exists() and "t1" in out_denovo.read_text()
    assert out_gg.exists() and "t2" in out_gg.read_text()


def test_run_defuse_splits_a_fused_transcript(tmp_path, monkeypatch):
    # Homology-stranded models from strand_correct are preferred over the raw GFF3.
    (tmp_path / "trinity-denovo.genome.stranded.gff3").write_text(_exon_gff3("t1"))
    monkeypatch.setattr(defuse, "run_cmd", lambda cmd, **kw: None)
    monkeypatch.setattr(defuse, "_resolve_diamond_db", lambda config: "db")
    monkeypatch.setattr(defuse, "_open_indexed_bam", lambda path: None)
    monkeypatch.setattr(
        defuse, "parse_positional_hits",
        lambda path: {
            "trinity-denovo.genome:t1": [_hit("P1", 1, 80), _hit("P2", 120, 200)]
        },
    )

    run_defuse(_config(tmp_path))

    text = (tmp_path / "trinity-denovo.genome.defuse.gff3").read_text()
    assert "t1.d1" in text and "t1.d2" in text  # split into two
    audit = (tmp_path / "defuse.tsv").read_text().splitlines()
    assert audit[0].startswith("set\ttid")
    assert any("t1" in row and "P1,P2" in row for row in audit[1:])
