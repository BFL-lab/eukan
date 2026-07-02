"""Unit tests for eukan.assembly.combinr (transcript consolidation)."""

from __future__ import annotations

from pathlib import Path

from eukan.assembly import combinr
from eukan.assembly.combinr import (
    _parse_attrs,
    _parse_combinr_gff3,
    _Transcript,
    _write_transcript_alignments_and_hints,
    _write_transcript_fasta,
    run_combinr,
)
from eukan.infra.artifacts import Artifact
from eukan.settings import AssemblyConfig

COMBINR_GFF = """##gff-version 3
chr1\tcombinr\tgene\t100\t400\t.\t+\t.\tID=g1
chr1\tcombinr\tmRNA\t100\t400\t.\t+\t.\tID=m1;Parent=g1;contains=a,b
chr1\tcombinr\texon\t300\t400\t.\t+\t.\tID=m1.exon2;Parent=m1
chr1\tcombinr\texon\t100\t200\t.\t+\t.\tID=m1.exon1;Parent=m1
"""


def _config(tmp_path, **kw):
    return AssemblyConfig(
        genome=tmp_path / "genome.fa", work_dir=tmp_path, num_cpu=2,
        max_intron_len=5000, **kw,
    )


# --- parsing ---------------------------------------------------------------


def test_parse_attrs():
    assert _parse_attrs("ID=m1;Parent=g1;contains=a,b") == {
        "ID": "m1", "Parent": "g1", "contains": "a,b"
    }


def test_parse_combinr_gff3_groups_and_sorts_exons(tmp_path):
    gff = tmp_path / "c.gff3"
    gff.write_text(COMBINR_GFF)
    (tx,) = _parse_combinr_gff3(gff)
    assert tx.tid == "m1" and tx.chrom == "chr1" and tx.strand == "+"
    assert tx.exons == [(100, 200), (300, 400)]  # exon-sorted
    assert tx.start == 100 and tx.end == 400


# --- artifact writers ------------------------------------------------------


def test_write_transcript_alignments_and_hints(tmp_path):
    txs = [_Transcript("m1", "chr1", "+", [(100, 200), (300, 400)])]
    gff = tmp_path / "nr.gff3"
    hints = tmp_path / "hints.gff"
    _write_transcript_alignments_and_hints(txs, gff, hints)

    rows = [ln.split("\t") for ln in gff.read_text().splitlines()]
    assert all(r[1] == "combinr-assembly" and r[2] == "exon" for r in rows)
    assert [r[3:5] for r in rows] == [["100", "200"], ["300", "400"]]
    assert all("Parent=m1" in r[8] for r in rows)
    assert "ID=m1:exon:1" in rows[0][8] and "ID=m1:exon:2" in rows[1][8]

    hint_rows = hints.read_text().splitlines()
    assert all("pri=3;src=E;group=m1" in h for h in hint_rows)


def test_write_transcript_fasta_splices_and_revcomps(tmp_path):
    (tmp_path / "genome.fa").write_text(">chr1\nAAAACCCCGGGGTTTTACGT\n")
    plus = _Transcript("p", "chr1", "+", [(1, 4), (9, 12)])   # AAAA + GGGG
    minus = _Transcript("m", "chr1", "-", [(1, 4), (9, 12)])  # revcomp(AAAAGGGG)
    out = tmp_path / "tx.fasta"
    _write_transcript_fasta([plus, minus], tmp_path / "genome.fa", out)

    recs = dict(_read_fasta(out))
    assert recs["p"] == "AAAAGGGG"
    assert recs["m"] == "CCCCTTTT"


def _read_fasta(path):
    name, seq = None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(seq)
            name, seq = line[1:].strip(), []
        else:
            seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


# --- command construction --------------------------------------------------


def test_run_combinr_assemble_command(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(combinr, "run_cmd", lambda cmd, **kw: calls.append((cmd, kw)))
    combinr._run_combinr_assemble(
        _config(tmp_path), [tmp_path / "a.bam", tmp_path / "b.bam"], tmp_path / "out.gff3"
    )
    cmd, kw = calls[0]
    assert cmd[:2] == ["combinr", "assemble"]
    assert cmd.count("-i") == 2
    assert cmd[cmd.index("--format") + 1] == "gff3"
    assert cmd[cmd.index("--max-intron") + 1] == "5000"
    # stringent-overlap defaults to 0.0 (off) unless configured.
    assert cmd[cmd.index("--stringent-overlap") + 1] == "0.0"
    assert kw["out_file"] == "out.gff3"


def test_run_combinr_assemble_passes_stringent_overlap(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(combinr, "run_cmd", lambda cmd, **kw: calls.append(cmd))
    cfg = _config(tmp_path, combinr_stringent_overlap=30.0)
    combinr._run_combinr_assemble(cfg, [tmp_path / "a.bam"], tmp_path / "o.gff3")
    assert calls[0][calls[0].index("--stringent-overlap") + 1] == "30.0"


def test_combinr_path_override(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(combinr, "run_cmd", lambda cmd, **kw: calls.append(cmd))
    cfg = _config(tmp_path, combinr_path=tmp_path / "bin" / "combinr")
    combinr._run_combinr_assemble(cfg, [tmp_path / "a.bam"], tmp_path / "o.gff3")
    assert calls[0][0] == str(tmp_path / "bin" / "combinr")


# --- run_combinr integration (combinr mocked) ------------------------------

_CUT_MODELS = (
    "trinity-denovo.genome.cut.gff3",
    "trinity-gg.genome.cut.gff3",
)


def _setup_run(tmp_path):
    (tmp_path / "genome.fa").write_text(">chr1\n" + "ACGT" * 250 + "\n")  # 1000 bp
    for f in _CUT_MODELS:
        (tmp_path / f).write_text("##gff-version 3\n")  # non-empty; combinr mocked


def _fake_assemble(calls):
    def fake(config, inputs, out_gff):
        calls.append({Path(p).name for p in inputs})
        out_gff.write_text(COMBINR_GFF)  # the consolidated combinr output
    return fake


def test_run_combinr_consolidates_cut_models(tmp_path, monkeypatch):
    _setup_run(tmp_path)
    calls: list[set] = []
    monkeypatch.setattr(combinr, "_run_combinr_assemble", _fake_assemble(calls))

    run_combinr(_config(tmp_path))

    # single combinr run over the cut models (both Trinity tracks)
    assert len(calls) == 1
    assert calls[0] == set(_CUT_MODELS)
    assert (tmp_path / Artifact.NR_TRANSCRIPTS_GFF).exists()
    assert (tmp_path / Artifact.NR_TRANSCRIPTS_FASTA).exists()
    assert (tmp_path / Artifact.RNASEQ_HINTS).exists()


def test_run_combinr_skips_empty_cut_models(tmp_path, monkeypatch):
    _setup_run(tmp_path)
    # An empty (zero-byte) cut model is ignored, not fed to combinr.
    (tmp_path / "trinity-gg.genome.cut.gff3").write_text("")
    calls: list[set] = []
    monkeypatch.setattr(combinr, "_run_combinr_assemble", _fake_assemble(calls))

    run_combinr(_config(tmp_path))

    assert calls[0] == {"trinity-denovo.genome.cut.gff3"}


def test_run_combinr_errors_without_inputs(tmp_path, monkeypatch):
    (tmp_path / "genome.fa").write_text(">chr1\nACGT\n")
    monkeypatch.setattr(combinr, "_run_combinr_assemble", _fake_assemble([]))
    import pytest
    with pytest.raises(FileNotFoundError):
        run_combinr(_config(tmp_path))
