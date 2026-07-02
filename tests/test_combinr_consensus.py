"""Unit tests for eukan.annotation.combinr_consensus (the consensus engine)."""

from __future__ import annotations

import os
import shutil
from collections import Counter
from pathlib import Path

import pytest

from eukan.annotation import combinr_consensus as cc
from eukan.annotation import consensus as cons
from eukan.annotation.combinr_consensus import (
    _chains_to_match,
    _stage_combinr_inputs,
    run_combinr_consensus,
)
from eukan.settings import PipelineConfig


def _resolve_combinr() -> Path | None:
    """Locate a combinr binary: $COMBINR_BIN, else PATH. None → skip e2e test."""
    env = os.environ.get("COMBINR_BIN")
    if env and Path(env).exists():
        return Path(env)
    on_path = shutil.which("combinr")
    return Path(on_path) if on_path else None


_COMBINR = _resolve_combinr()


def _config(tmp_path: Path, **kw) -> PipelineConfig:
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\nACGT\n")
    prot = tmp_path / "proteins.faa"
    prot.write_text(">p\nM\n")
    return PipelineConfig(
        genome=genome, proteins=[prot], work_dir=tmp_path, num_cpu=2, **kw
    )


def _with_transcripts(tmp_path: Path, nr_gff: Path, **kw) -> PipelineConfig:
    """Config whose has_transcripts is True (all three artifacts set explicitly)."""
    fa = tmp_path / "nr_transcripts.fasta"
    fa.write_text(">m1\nACGT\n")
    hints = tmp_path / "hints_rnaseq.gff"
    hints.write_text("")
    return _config(
        tmp_path,
        transcripts_fasta=fa, transcripts_gff=nr_gff, rnaseq_hints=hints, **kw
    )


def _read_rows(path: Path) -> list[list[str]]:
    return [
        ln.split("\t")
        for ln in path.read_text().splitlines()
        if ln.strip() and not ln.startswith("#")
    ]


# --- _chains_to_match ------------------------------------------------------


def test_chains_to_match_exon_groups_and_targets(tmp_path):
    src = tmp_path / "nr.gff3"
    src.write_text(
        "##gff-version 3\n"
        "chr1\tcombinr-assembly\texon\t300\t400\t.\t+\t.\tID=m1:exon:2;Parent=m1\n"
        "chr1\tcombinr-assembly\texon\t100\t200\t.\t+\t.\tID=m1:exon:1;Parent=m1\n"
        "chr2\tcombinr-assembly\texon\t10\t19\t.\t-\t.\tID=m2:exon:1;Parent=m2\n"
    )
    n = _chains_to_match(src, tmp_path / "out.gff3", feature_type="exon", match_type="cDNA_match")
    assert n == 2
    rows = _read_rows(tmp_path / "out.gff3")
    # every row is a Target-bearing cDNA_match keeping the source token
    assert all(r[1] == "combinr-assembly" and r[2] == "cDNA_match" for r in rows)
    assert all("Target=" in r[8] for r in rows)
    # m1's exons emitted in genomic order with cumulative target coords
    m1 = [r for r in rows if "Parent" not in r[8] and r[8].startswith("ID=m1")]
    assert [r[3:5] for r in m1] == [["100", "200"], ["300", "400"]]
    assert "Target=m1 1 101" in m1[0][8] and "Target=m1 102 202" in m1[1][8]
    # strand preserved per chain
    m2 = [r for r in rows if r[8].startswith("ID=m2")]
    assert m2[0][6] == "-"


def test_chains_to_match_cds_for_proteins(tmp_path):
    src = tmp_path / "prot.gff3"
    src.write_text(
        "##gff-version 3\n"
        "chr1\tprot_align\tgene\t1\t300\t.\t+\t.\tID=p1\n"
        "chr1\tprot_align\tmRNA\t1\t300\t.\t+\t.\tID=p1.t1;Parent=p1\n"
        "chr1\tprot_align\tCDS\t1\t100\t.\t+\t0\tID=c1;Parent=p1.t1\n"
        "chr1\tprot_align\tCDS\t200\t300\t.\t+\t0\tID=c2;Parent=p1.t1\n"
    )
    n = _chains_to_match(
        src, tmp_path / "p.match.gff3",
        feature_type="CDS", match_type="nucleotide_to_protein_match",
    )
    assert n == 1  # both CDS grouped under one Parent chain
    rows = _read_rows(tmp_path / "p.match.gff3")
    assert len(rows) == 2
    assert all(r[2] == "nucleotide_to_protein_match" and "Target=p1.t1" in r[8] for r in rows)


# --- _stage_combinr_inputs -------------------------------------------------


def _prot_gff(path: Path) -> Path:
    path.write_text(
        "##gff-version 3\n"
        "chr1\tprot_align\tmRNA\t1\t300\t.\t+\t.\tID=p1;Parent=pg1\n"
        "chr1\tprot_align\tCDS\t1\t300\t.\t+\t0\tID=c1;Parent=p1\n"
    )
    return path


def _pred_gff(path: Path, source: str) -> Path:
    path.write_text(
        "##gff-version 3\n"
        f"chr1\t{source}\tgene\t1\t300\t.\t+\t.\tID=g1\n"
        f"chr1\t{source}\tmRNA\t1\t300\t.\t+\t.\tID=g1.t1;Parent=g1\n"
        f"chr1\t{source}\tCDS\t1\t300\t.\t+\t0\tID=cds1;Parent=g1.t1\n"
    )
    return path


def _nr_gff(path: Path) -> Path:
    path.write_text(
        "##gff-version 3\n"
        "chr1\tcombinr-assembly\texon\t100\t200\t.\t+\t.\tID=m1:exon:1;Parent=m1\n"
    )
    return path


def test_stage_combinr_inputs_full(tmp_path):
    sdir = tmp_path / "consensus_models"
    sdir.mkdir()
    src = tmp_path / "src"
    src.mkdir()
    prot = _prot_gff(src / "prot.gff3")
    aug = _pred_gff(src / "augustus.gff3", "augustus")
    nr = _nr_gff(src / "nr.gff3")

    cfg = _with_transcripts(tmp_path, nr)
    have_t = _stage_combinr_inputs(cfg, sdir, [prot, aug], nr)
    assert have_t is True

    weights = (sdir / "weights.txt").read_text()
    assert "PROTEIN\tprot_align\t2" in weights
    assert "ABINITIO_PREDICTION\taugustus\t1" in weights
    assert "TRANSCRIPT\tcombinr-assembly\t3" in weights

    # proteins are converted (not concatenated into gene_predictions)
    assert (sdir / "prot.match.gff3").exists()
    assert (sdir / "transcripts.match.gff3").exists()
    preds = (sdir / "gene_predictions.gff3").read_text()
    assert "augustus\tCDS" in preds
    assert "prot_align" not in preds


def test_stage_combinr_inputs_without_transcripts(tmp_path):
    sdir = tmp_path / "consensus_models"
    sdir.mkdir()
    src = tmp_path / "src"
    src.mkdir()
    prot = _prot_gff(src / "prot.gff3")
    aug = _pred_gff(src / "augustus.gff3", "augustus")

    cfg = _config(tmp_path)  # has_transcripts is False
    have_t = _stage_combinr_inputs(cfg, sdir, [prot, aug], None)
    assert have_t is False
    assert not (sdir / "transcripts.match.gff3").exists()
    assert "TRANSCRIPT" not in (sdir / "weights.txt").read_text()


# --- run_combinr_consensus command construction ----------------------------


def _capture_run_cmd(monkeypatch) -> list[tuple[list[str], dict]]:
    calls: list[tuple[list[str], dict]] = []
    monkeypatch.setattr(cc, "run_cmd", lambda cmd, **kw: calls.append((cmd, kw)))
    return calls


def test_run_combinr_consensus_command_with_transcripts(tmp_path, monkeypatch):
    sdir = tmp_path / "consensus_models"
    sdir.mkdir()
    src = tmp_path / "src"
    src.mkdir()
    prot = _prot_gff(src / "prot.gff3")
    aug = _pred_gff(src / "augustus.gff3", "augustus")
    nr = _nr_gff(src / "nr.gff3")
    cfg = _with_transcripts(tmp_path, nr, genetic_code="6")

    calls = _capture_run_cmd(monkeypatch)
    run_combinr_consensus(cfg, sdir, [prot, aug], transcripts=nr)

    cmd, kw = calls[0]
    assert cmd[:2] == ["combinr", "consensus"]
    assert cmd[cmd.index("--genetic-code") + 1] == "6"
    assert cmd[cmd.index("--gene-predictions") + 1] == "gene_predictions.gff3"
    assert cmd[cmd.index("--protein-alignments") + 1] == "prot.match.gff3"
    assert cmd[cmd.index("--transcript-alignments") + 1] == "transcripts.match.gff3"
    assert "--alt-splice" in cmd
    # stringent-overlap accompanies the alt-splice isoform grouping; 0.0 = off default.
    assert cmd[cmd.index("--stringent-overlap") + 1] == "0.0"
    assert kw["out_file"] == "consensus_models.gff3"


def test_run_combinr_consensus_passes_stringent_overlap(tmp_path, monkeypatch):
    sdir = tmp_path / "consensus_models"
    sdir.mkdir()
    src = tmp_path / "src"
    src.mkdir()
    prot = _prot_gff(src / "prot.gff3")
    aug = _pred_gff(src / "augustus.gff3", "augustus")
    nr = _nr_gff(src / "nr.gff3")
    cfg = _with_transcripts(tmp_path, nr, combinr_stringent_overlap=30.0)

    calls = _capture_run_cmd(monkeypatch)
    run_combinr_consensus(cfg, sdir, [prot, aug], transcripts=nr)

    cmd, _ = calls[0]
    assert cmd[cmd.index("--stringent-overlap") + 1] == "30.0"


def test_run_combinr_consensus_command_without_transcripts(tmp_path, monkeypatch):
    sdir = tmp_path / "consensus_models"
    sdir.mkdir()
    src = tmp_path / "src"
    src.mkdir()
    prot = _prot_gff(src / "prot.gff3")
    aug = _pred_gff(src / "augustus.gff3", "augustus")
    cfg = _config(tmp_path)

    calls = _capture_run_cmd(monkeypatch)
    run_combinr_consensus(cfg, sdir, [prot, aug], transcripts=None)

    cmd, _ = calls[0]
    assert "--transcript-alignments" not in cmd
    assert "--alt-splice" not in cmd
    assert "--stringent-overlap" not in cmd  # rides with the alt-splice grouping only
    assert "--protein-alignments" in cmd  # protein evidence still passed


def test_combinr_path_override(tmp_path, monkeypatch):
    sdir = tmp_path / "consensus_models"
    sdir.mkdir()
    src = tmp_path / "src"
    src.mkdir()
    prot = _prot_gff(src / "prot.gff3")
    binpath = tmp_path / "bin" / "combinr"
    binpath.parent.mkdir()
    binpath.write_text("")
    cfg = _config(tmp_path, combinr_path=binpath)

    calls = _capture_run_cmd(monkeypatch)
    run_combinr_consensus(cfg, sdir, [prot], transcripts=None)
    assert calls[0][0][0] == str(binpath)


# --- integration: combinr output through eukan's prettify tail --------------

# A compact slice of real `combinr consensus --alt-splice` output: one locus
# with a consensus mRNA plus a wider transcript isoform carrying 5'/3' UTR, and
# combinr's mRNA attributes (support/score/low_support/contains).
_COMBINR_ALTSPLICE_GFF = """\
chr1\tcombinr\tgene\t60\t260\t.\t+\t.\tID=consensus.chr1.g1
chr1\tcombinr\tmRNA\t100\t198\t.\t+\t.\tID=consensus.chr1.g1.consensus;Parent=consensus.chr1.g1;support=consensus;score=1234.5;score_ratio=10.0;coding_length=99;low_support=false;partial5=false;partial3=false
chr1\tcombinr\texon\t100\t198\t.\t+\t.\tID=consensus.chr1.g1.consensus.exon1;Parent=consensus.chr1.g1.consensus
chr1\tcombinr\tCDS\t100\t198\t.\t+\t0\tID=consensus.chr1.g1.consensus.cds1;Parent=consensus.chr1.g1.consensus
chr1\tcombinr\tmRNA\t60\t260\t.\t+\t.\tID=consensus.chr1.g1.iso1;Parent=consensus.chr1.g1;support=transcript_isoform;coding_altered=false;sources=transcript_alignments.gff3;contains=t1,t2
chr1\tcombinr\tfive_prime_UTR\t60\t99\t.\t+\t.\tID=consensus.chr1.g1.iso1.utr5;Parent=consensus.chr1.g1.iso1
chr1\tcombinr\texon\t60\t260\t.\t+\t.\tID=consensus.chr1.g1.iso1.exon1;Parent=consensus.chr1.g1.iso1
chr1\tcombinr\tCDS\t100\t198\t.\t+\t0\tID=consensus.chr1.g1.iso1.cds1;Parent=consensus.chr1.g1.iso1
chr1\tcombinr\tthree_prime_UTR\t199\t260\t.\t+\t.\tID=consensus.chr1.g1.iso1.utr3;Parent=consensus.chr1.g1.iso1
"""


def _feature_spans(path: Path) -> tuple[dict, dict, dict]:
    """Return (mRNA->gene, mRNA->exon spans, mRNA->CDS spans) from a GFF3."""
    parent, exons, cds = {}, {}, {}
    for ln in path.read_text().splitlines():
        if ln.startswith("#") or not ln.strip():
            continue
        c = ln.split("\t")
        a = dict(p.split("=", 1) for p in c[8].split(";") if "=" in p)
        if c[2] == "mRNA":
            parent[a["ID"]] = a.get("Parent")
        elif c[2] == "exon":
            exons.setdefault(a.get("Parent"), []).append((int(c[3]), int(c[4])))
        elif c[2] == "CDS":
            cds.setdefault(a.get("Parent"), []).append((int(c[3]), int(c[4])))
    return parent, exons, cds


def test_combinr_handoff_preserves_isoforms_and_implicit_utr(tmp_path, monkeypatch):
    """The prettify tail must keep multi-isoform genes and retain UTR as exon
    overhang beyond CDS (combinr's explicit UTR features collapse into exons,
    the canonical NCBI representation)."""
    def fake(config, sdir, evidence, *, transcripts=None):
        (config.work_dir / "consensus_models" / "consensus_models.gff3").write_text(
            _COMBINR_ALTSPLICE_GFF
        )
        return sdir / "consensus_models.gff3"

    monkeypatch.setattr(cons, "run_combinr_consensus", fake)
    cfg = _config(tmp_path)
    ev = tmp_path / "prot.gff3"
    ev.touch()

    final = cons.build_consensus_models(cfg, ev, transcripts=None)

    parent, exons, cds = _feature_spans(final)
    # one gene, two isoforms preserved
    assert len(parent) == 2
    assert len(set(parent.values())) == 1
    # at least one isoform keeps exon extent beyond its CDS (implicit UTR)
    overhang = [
        m for m in cds
        if m in exons and (
            min(s for s, _ in exons[m]) < min(s for s, _ in cds[m])
            or max(e for _, e in exons[m]) > max(e for _, e in cds[m])
        )
    ]
    assert overhang, "UTR lost: no mRNA has exons extending past its CDS"
    # combinr's isoform provenance attribute survives prettification
    assert "contains=" in final.read_text()


@pytest.mark.skipif(_COMBINR is None, reason="combinr binary not found (set COMBINR_BIN or add to PATH)")
def test_combinr_consensus_end_to_end(tmp_path):
    """Real combinr binary through build_consensus_models on minimal inputs.

    A transcript wider than the predicted CDS must yield an isoform with UTR;
    asserts the staged + converted evidence drives combinr to a valid final
    GFF3 (the path CI cannot exercise without the Rust binary)."""
    seq = list("A" * 600)
    orf = "ATG" + "AAG" * 31 + "TAA"  # clean ORF, no internal in-frame stop
    for i, ch in enumerate(orf):
        seq[99 + i] = ch
    (tmp_path / "genome.fa").write_text(">chr1\n" + "".join(seq) + "\n")
    (tmp_path / "proteins.faa").write_text(">p\nMK\n")
    (tmp_path / "augustus.gff3").write_text(
        "##gff-version 3\n"
        "chr1\taugustus\tgene\t100\t198\t.\t+\t.\tID=g1\n"
        "chr1\taugustus\tmRNA\t100\t198\t.\t+\t.\tID=m1;Parent=g1\n"
        "chr1\taugustus\texon\t100\t198\t.\t+\t.\tID=e1;Parent=m1\n"
        "chr1\taugustus\tCDS\t100\t198\t.\t+\t0\tID=c1;Parent=m1\n"
    )
    (tmp_path / "prot.gff3").write_text(
        "##gff-version 3\n"
        "chr1\tprot_align\tmRNA\t100\t198\t.\t+\t.\tID=pm1;Parent=pg1\n"
        "chr1\tprot_align\tCDS\t100\t198\t.\t+\t0\tID=pc1;Parent=pm1\n"
    )
    nr = tmp_path / "nr_transcripts.gff3"
    nr.write_text(
        "##gff-version 3\n"
        "chr1\tcombinr-assembly\texon\t60\t260\t.\t+\t.\tID=t1:exon:1;Parent=t1\n"
    )
    (tmp_path / "nr_transcripts.fasta").write_text(">t1\nACGT\n")
    (tmp_path / "hints.gff").write_text("")

    cfg = PipelineConfig(
        genome=tmp_path / "genome.fa", proteins=[tmp_path / "proteins.faa"],
        work_dir=tmp_path, num_cpu=2, combinr_path=_COMBINR,
        transcripts_fasta=tmp_path / "nr_transcripts.fasta", transcripts_gff=nr,
        rnaseq_hints=tmp_path / "hints.gff",
    )

    final = cons.build_consensus_models(cfg, tmp_path / "prot.gff3", tmp_path / "augustus.gff3", transcripts=nr)

    assert final.exists()
    types = Counter(
        ln.split("\t")[2]
        for ln in final.read_text().splitlines()
        if ln.strip() and not ln.startswith("#")
    )
    assert types["gene"] >= 1 and types["CDS"] >= 1 and types["exon"] >= 1
    # the wider transcript must have produced an isoform whose exon overhangs CDS
    _parent, exons, cds = _feature_spans(final)
    assert any(
        m in exons and (
            min(s for s, _ in exons[m]) < min(s for s, _ in cds[m])
            or max(e for _, e in exons[m]) > max(e for _, e in cds[m])
        )
        for m in cds
    ), "expected a UTR-bearing isoform from the wider transcript"
