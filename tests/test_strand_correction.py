"""Unit + integration tests for eukan.assembly.strand_correction."""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import pysam

from eukan.assembly import strand_correction as sc
from eukan.assembly.jaccard import _parse_transcript_models, _Tx
from eukan.infra.genome import ContigIndex
from eukan.settings import AssemblyConfig

# --------------------------------------------------------------------------- #
# Fixtures / builders
# --------------------------------------------------------------------------- #


def _seq_with(length: int, **pos_bases: str) -> str:
    """An ``A``-filled sequence with specific 1-based positions overridden."""
    s = ["A"] * length
    for pos, base in pos_bases.items():
        s[int(pos[1:]) - 1] = base  # keys like p41 -> 1-based 41
    return "".join(s)


# A canonical GT-AG intron for exons (1,40),(51,90): donor at 41/42, acceptor 49/50.
_GT_AG = {"p41": "G", "p42": "T", "p49": "A", "p50": "G"}
# Its plus-genome reverse-complement twin CT-AC (a mislabelled minus gene).
_CT_AC = {"p41": "C", "p42": "T", "p49": "A", "p50": "C"}


def _write_fasta(path: Path, contigs: dict[str, str]) -> None:
    with open(path, "w") as fh:
        for name, seq in contigs.items():
            fh.write(f">{name}\n{seq}\n")


def _gtf_tx(chrom: str, tid: str, strand: str, exons: list[tuple[int, int]]) -> str:
    lo, hi = exons[0][0], exons[-1][1]
    rows = [f'{chrom}\tStringTie\ttranscript\t{lo}\t{hi}\t.\t{strand}\t.\ttranscript_id "{tid}";']
    rows += [
        f'{chrom}\tStringTie\texon\t{s}\t{e}\t.\t{strand}\t.\ttranscript_id "{tid}";'
        for s, e in exons
    ]
    return "\n".join(rows) + "\n"


def _make_bam(path: Path, reads, ref="chr1", ref_len=100_000) -> None:
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": ref, "LN": ref_len}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for name, flag, start, cigar, seq in reads:
            s = pysam.AlignedSegment(out.header)
            s.query_name, s.flag, s.reference_id = name, flag, 0
            s.reference_start, s.mapping_quality = start, 60
            s.cigartuples, s.query_sequence = cigar, seq
            out.write(s)


# --------------------------------------------------------------------------- #
# Pure helpers
# --------------------------------------------------------------------------- #


def test_parse_hits(tmp_path):
    p = tmp_path / "hits.tsv"
    p.write_text(
        "# comment\n"
        "trinity-gg.genome:T1\tsp|P1\t120.5\t2\n"     # best for T1 (forward frame)
        "trinity-gg.genome:T1\tsp|P9\t90.0\t-1\n"     # lower bitscore, same query -> ignored
        "trinity-denovo.genome:N2\tsp|P2\t77.0\t-3\n"  # reverse-frame hit
        "shortrow\n"
    )
    assert sc.parse_hits(p) == {"trinity-gg.genome:T1": 2, "trinity-denovo.genome:N2": -3}


def test_coding_strand():
    assert sc._coding_strand("+", 2) == "+"
    assert sc._coding_strand("+", -1) == "-"
    assert sc._coding_strand("-", 1) == "-"   # reverse-labelled query, forward frame
    assert sc._coding_strand("-", -2) == "+"  # antisense-assembled rescue
    assert sc._coding_strand(".", 3) == "+"   # '.' stitched forward, like '+'
    assert sc._coding_strand(".", -3) == "-"


def test_introns_of():
    assert sc.introns_of([(1, 40), (51, 90)]) == [(40, 50)]
    assert sc.introns_of([(1, 40), (51, 90), (101, 140)]) == [(40, 50), (90, 100)]
    assert sc.introns_of([(1, 90)]) == []


def test_rc_swap():
    assert sc._rc_swap("GT-AG") == "CT-AC"
    assert sc._rc_swap("CT-AC") == "GT-AG"
    assert sc._rc_swap("CG-AG") == "CT-CG"  # non-canonical diplonemid twin


def test_consensus_on_strand():
    assert sc.consensus_on_strand("GT-AG", "+") == "GT-AG"
    assert sc.consensus_on_strand("GT-AG", ".") == "GT-AG"   # extracted forward
    assert sc.consensus_on_strand("CT-AC", "-") == "GT-AG"   # minus reads as coding
    assert sc.consensus_on_strand(None, "+") is None


def test_pick_consensus_dominant():
    assert sc._pick_consensus(Counter({"GT-AG": 60, "GC-AG": 2}), 50) == ("GT-AG", "CT-AC")


def test_pick_consensus_noncanonical_dominant():
    assert sc._pick_consensus(Counter({"CG-AG": 60}), 50) == ("CG-AG", "CT-CG")


def test_pick_consensus_thin_falls_back_to_canonical():
    assert sc._pick_consensus(Counter({"CG-AG": 3}), 50) == ("GT-AG", "CT-AC")
    assert sc._pick_consensus(Counter(), 1) == ("GT-AG", "CT-AC")


# --------------------------------------------------------------------------- #
# _decide
# --------------------------------------------------------------------------- #


def _decide_contigs(tmp_path) -> ContigIndex:
    _write_fasta(tmp_path / "g.fa", {
        "plus": _seq_with(140, **_GT_AG),
        "twin": _seq_with(140, **_CT_AC),
        "mixed": _seq_with(  # intron1 GT-AG (votes +), intron2 CT-AC (votes -)
            160, **_GT_AG, **{"p91": "C", "p92": "T", "p99": "A", "p100": "C"}
        ),
    })
    return ContigIndex(tmp_path / "g.fa")


def test_decide_hit_keeps_label(tmp_path):
    with _decide_contigs(tmp_path) as c:
        tx = _Tx("a", "plus", "+", "src", [(1, 40), (51, 90)])
        assert sc._decide(tx, "+", "GT-AG", "CT-AC", c) == ("+", "keep")


def test_decide_hit_resolves_dot_to_plus(tmp_path):
    with _decide_contigs(tmp_path) as c:
        tx = _Tx("a", "plus", ".", "src", [(1, 40), (51, 90)])
        assert sc._decide(tx, "+", "GT-AG", "CT-AC", c) == ("+", "assign")


def test_decide_hit_flips_antisense(tmp_path):
    # Fix A: a frame-derived coding strand opposite the label flips the model
    # (the NODE_574 case: antisense-assembled contig, label '-', coding '+').
    with _decide_contigs(tmp_path) as c:
        tx = _Tx("a", "plus", "-", "src", [(1, 40), (51, 90)])
        assert sc._decide(tx, "+", "GT-AG", "CT-AC", c) == ("+", "flip")


def test_decide_flips_reverse_twin(tmp_path):
    with _decide_contigs(tmp_path) as c:
        tx = _Tx("a", "twin", "+", "src", [(1, 40), (51, 90)])
        assert sc._decide(tx, None, "GT-AG", "CT-AC", c) == ("-", "flip")


def test_decide_canonical_votes_when_dominant_noncanonical(tmp_path):
    # Fix B: canonical GT-AG introns still vote even when the learned dominant is a
    # non-canonical motif (CG-AG) they don't match -- previously this was ambiguous.
    with _decide_contigs(tmp_path) as c:
        tx = _Tx("a", "plus", "-", "src", [(1, 40), (51, 90)])  # GT-AG intron
        assert sc._decide(tx, None, "CG-AG", "CT-CG", c) == ("+", "flip")


def test_decide_canonical_twin_votes_minus(tmp_path):
    # Fix B: the +genome CT-AC twin votes '-' even under a non-canonical dominant.
    with _decide_contigs(tmp_path) as c:
        tx = _Tx("a", "twin", "+", "src", [(1, 40), (51, 90)])  # CT-AC intron
        assert sc._decide(tx, None, "CG-AG", "CT-CG", c) == ("-", "flip")


def test_decide_assigns_dot_from_consensus(tmp_path):
    with _decide_contigs(tmp_path) as c:
        tx = _Tx("a", "plus", ".", "src", [(1, 40), (51, 90)])
        assert sc._decide(tx, None, "GT-AG", "CT-AC", c) == ("+", "assign")


def test_decide_monoexon_left_alone(tmp_path):
    with _decide_contigs(tmp_path) as c:
        tx = _Tx("a", "plus", "+", "src", [(1, 90)])
        assert sc._decide(tx, None, "GT-AG", "CT-AC", c) == ("+", "mono-exon")


def test_decide_ambiguous_left_alone(tmp_path):
    with _decide_contigs(tmp_path) as c:
        tx = _Tx("a", "mixed", "+", "src", [(1, 40), (51, 90), (101, 140)])
        assert sc._decide(tx, None, "GT-AG", "CT-AC", c) == ("+", "ambiguous")


# --------------------------------------------------------------------------- #
# run_strand_correction
# --------------------------------------------------------------------------- #


def _config(tmp_path, **kw):
    _write_fasta(tmp_path / "genome.fa", {
        "chrA": _seq_with(140, **_GT_AG),
        "chrB": _seq_with(140, **_CT_AC),
        "chrC": _seq_with(140, **_GT_AG),
        "chrD": "A" * 140,
    })
    return AssemblyConfig(genome=tmp_path / "genome.fa", work_dir=tmp_path, num_cpu=2, **kw)


def test_run_corrects_strands(tmp_path, monkeypatch):
    (tmp_path / "db.dmnd").write_text("x")  # prebuilt .dmnd → no makedb
    # The mapped Trinity genome track's raw GFF3 (built from the BAM in step 1);
    # here written directly so step 3 strand-corrects it into *.stranded.gff3.
    (tmp_path / "trinity-gg.genome.gff3").write_text(
        _gtf_tx("chrA", "T_confirmed_plus", "+", [(1, 40), (51, 90)])
        + _gtf_tx("chrB", "T_wrong", "+", [(1, 40), (51, 90)])
        + _gtf_tx("chrC", "T_dot", ".", [(1, 40), (51, 90)])
        + _gtf_tx("chrA", "T_antisense", "-", [(1, 40), (51, 90)])
        + _gtf_tx("chrD", "T_mono", "+", [(1, 90)])
    )

    def fake_run_cmd(cmd, **kw):
        if "blastx" in cmd:
            out = cmd[cmd.index("--out") + 1]
            Path(out).write_text(
                "trinity-gg.genome:T_confirmed_plus\tsp|P1\t150.0\t1\n"  # forward frame -> '+'
                "trinity-gg.genome:T_antisense\tsp|P2\t140.0\t-1\n"      # reverse on '-' -> flip '+'
            )

    monkeypatch.setattr(sc, "run_cmd", fake_run_cmd)

    config = _config(tmp_path, uniprot_db=tmp_path / "db.dmnd", min_strand_consensus=1)
    sc.run_strand_correction(config)

    out = tmp_path / "trinity-gg.genome.stranded.gff3"
    assert out.exists()
    strands = {m.tid: m.strand for m in _parse_transcript_models(out)}
    assert strands == {
        "T_confirmed_plus": "+",  # forward-frame hit → kept
        "T_wrong": "-",           # CT-AC twin → flipped by vote
        "T_dot": "+",             # consensus assigns +
        "T_antisense": "+",       # reverse-frame hit → homology flip
        "T_mono": "+",            # single exon → untouched
    }

    audit = (tmp_path / "strand_correction.tsv").read_text().splitlines()
    decisions = {row.split("\t")[1]: row.split("\t")[5] for row in audit[1:]}
    assert decisions == {
        "T_confirmed_plus": "keep", "T_wrong": "flip",
        "T_dot": "assign", "T_antisense": "flip", "T_mono": "mono-exon",
    }


def test_run_converts_denovo_bam_then_gates_on_no_uniprot(tmp_path, monkeypatch):
    """Even when disabled, each Trinity transcript->genome BAM is converted to GFF3
    for the SL cut."""
    _make_bam(
        tmp_path / "trinity-denovo.genome.bam",
        [("q1", 0, 100, [(0, 10), (3, 50), (0, 10)], "A" * 20)],
    )
    called: list[list[str]] = []
    monkeypatch.setattr(sc, "run_cmd", lambda cmd, **kw: called.append(cmd))

    config = _config(tmp_path, uniprot_db=None)  # disabled
    sc.run_strand_correction(config)

    assert (tmp_path / "trinity-denovo.genome.gff3").exists()         # always produced
    assert not (tmp_path / "trinity-gg.genome.stranded.gff3").exists()  # gated off
    assert not (tmp_path / "trinity-denovo.genome.stranded.gff3").exists()
    assert not any("blastx" in c for c in called)                    # diamond not run


def test_run_clears_stale_stranded_on_noop(tmp_path, monkeypatch):
    """A no-op run removes prior *.stranded.gff3 for BOTH Trinity tracks so sl_cut
    falls back to the fresh raw genome GFF3 instead of a stale stranded file from a
    run that had --uniprot."""
    monkeypatch.setattr(sc, "run_cmd", lambda cmd, **kw: None)
    (tmp_path / "trinity-gg.genome.stranded.gff3").write_text("##gff-version 3\n")
    (tmp_path / "trinity-denovo.genome.stranded.gff3").write_text("##gff-version 3\n")

    config = _config(tmp_path, uniprot_db=None)  # correction now disabled
    sc.run_strand_correction(config)

    assert not (tmp_path / "trinity-gg.genome.stranded.gff3").exists()
    assert not (tmp_path / "trinity-denovo.genome.stranded.gff3").exists()


def test_run_gates_on_strand_specific(tmp_path, monkeypatch):
    (tmp_path / "db.dmnd").write_text("x")
    (tmp_path / "trinity-gg.genome.gff3").write_text(
        _gtf_tx("chrA", "T1", "+", [(1, 40), (51, 90)])
    )
    called: list[list[str]] = []
    monkeypatch.setattr(sc, "run_cmd", lambda cmd, **kw: called.append(cmd))

    config = _config(tmp_path, uniprot_db=tmp_path / "db.dmnd", strand_specific="RF")
    sc.run_strand_correction(config)

    assert not (tmp_path / "trinity-gg.genome.stranded.gff3").exists()
    assert not any("blastx" in c for c in called)
