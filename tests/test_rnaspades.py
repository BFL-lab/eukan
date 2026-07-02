"""Unit tests for eukan.assembly.rnaspades (command construction, no real run)."""

from __future__ import annotations

import pytest

from eukan.assembly import rnaspades
from eukan.settings import AssemblyConfig


def _config(tmp_path, **kw):
    return AssemblyConfig(
        genome=tmp_path / "genome.fa",
        work_dir=tmp_path,
        num_cpu=8,
        memory_gb=16,
        phred_quality=33,
        **kw,
    )


def _capture(monkeypatch):
    """Replace run_cmd in the rnaspades module with a recorder."""
    calls: list[tuple[list[str], dict]] = []
    monkeypatch.setattr(rnaspades, "run_cmd", lambda cmd, **kw: calls.append((cmd, kw)))
    return calls


def test_paired_reads_command(tmp_path, monkeypatch):
    calls = _capture(monkeypatch)
    cfg = _config(
        tmp_path,
        left_reads=tmp_path / "R1.fq.gz",
        right_reads=tmp_path / "R2.fq.gz",
    )
    rnaspades.run_rnaspades(cfg)

    assert len(calls) == 1
    cmd, kw = calls[0]
    assert cmd[0] == "rnaspades.py"
    assert cmd[cmd.index("-1") + 1] == str(tmp_path / "R1.fq.gz")
    assert cmd[cmd.index("-2") + 1] == str(tmp_path / "R2.fq.gz")
    assert cmd[cmd.index("-t") + 1] == "8"
    assert cmd[cmd.index("-m") + 1] == "16"
    assert cmd[cmd.index("--phred-offset") + 1] == "33"
    assert cmd[cmd.index("-o") + 1] == str(tmp_path / "rnaspades_out")
    assert kw["cwd"] == tmp_path


def test_single_reads_command(tmp_path, monkeypatch):
    calls = _capture(monkeypatch)
    cfg = _config(tmp_path, single_reads=tmp_path / "reads.fq.gz")
    rnaspades.run_rnaspades(cfg)

    cmd, _ = calls[0]
    assert cmd[cmd.index("-s") + 1] == str(tmp_path / "reads.fq.gz")
    assert "-1" not in cmd and "-2" not in cmd


def test_skips_when_output_exists(tmp_path, monkeypatch):
    calls = _capture(monkeypatch)
    (tmp_path / "rnaspades.fasta").write_text(">t\nACGT\n")
    cfg = _config(tmp_path, single_reads=tmp_path / "reads.fq.gz")
    rnaspades.run_rnaspades(cfg)

    assert calls == []  # cached output → no invocation


def test_normalizes_transcripts_fasta(tmp_path, monkeypatch):
    """transcripts.fasta in the -o dir is renamed to rnaspades.fasta."""

    def fake_run(cmd, **kw):
        out_dir = tmp_path / "rnaspades_out"
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "transcripts.fasta").write_text(">t1\nACGTACGT\n")

    monkeypatch.setattr(rnaspades, "run_cmd", fake_run)
    cfg = _config(tmp_path, single_reads=tmp_path / "reads.fq.gz")
    rnaspades.run_rnaspades(cfg)

    final = tmp_path / "rnaspades.fasta"
    assert final.exists() and final.read_text().startswith(">t1")
    assert not (tmp_path / "rnaspades_out").exists()  # intermediates cleaned up


def test_no_reads_raises(tmp_path, monkeypatch):
    _capture(monkeypatch)
    cfg = _config(tmp_path)
    with pytest.raises(ValueError):
        rnaspades.run_rnaspades(cfg)
