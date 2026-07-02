"""Unit tests for eukan.assembly.stringtie (command construction)."""

from __future__ import annotations

from eukan.assembly import stringtie
from eukan.settings import AssemblyConfig


def _config(tmp_path, **kw):
    kw.setdefault("aligner", "segemehl")
    return AssemblyConfig(
        genome=tmp_path / "genome.fa", work_dir=tmp_path, num_cpu=4, **kw,
    )


def test_run_stringtie_bounds_segemehl_bam(tmp_path, monkeypatch):
    cmds: list[list[str]] = []
    monkeypatch.setattr(stringtie, "run_cmd", lambda cmd, **kw: cmds.append(cmd))
    split_calls: list[tuple] = []

    def fake_split(in_bam, out_bam, *, max_intron_len, num_cpu=1):
        split_calls.append((in_bam, out_bam, max_intron_len))
        out_bam.write_text("bounded")  # so the post-run unlink is exercised
        return 3

    monkeypatch.setattr(stringtie, "split_long_introns", fake_split)

    stringtie.run_stringtie(_config(tmp_path))

    in_bam, out_bam, mi = split_calls[0]
    assert in_bam == tmp_path / "segemehl_Aligned.sortedByCoord.out.bam"
    assert out_bam == tmp_path / "stringtie_input.bam"
    assert mi == 5000
    (cmd,) = cmds
    assert cmd[0] == "stringtie"
    assert cmd[1] == str(tmp_path / "stringtie_input.bam")  # reads the bounded copy
    assert cmd[cmd.index("-p") + 1] == "4"
    assert cmd[cmd.index("-o") + 1] == "stringtie.gtf"
    # stringency knobs (defaults raised above StringTie's -c 1 / -f 0.01; -j at its
    # default of 1 but passed explicitly so it can be tuned)
    assert cmd[cmd.index("-c") + 1] == "1.5"
    assert cmd[cmd.index("-f") + 1] == "0.1"
    assert cmd[cmd.index("-j") + 1] == "1.0"
    # unstranded: no strand flag
    assert "--rf" not in cmd and "--fr" not in cmd
    assert not (tmp_path / "stringtie_input.bam").exists()  # disposable copy removed


def test_run_stringtie_stringency_from_config(tmp_path, monkeypatch):
    cmds: list[list[str]] = []
    monkeypatch.setattr(stringtie, "run_cmd", lambda cmd, **kw: cmds.append(cmd))
    monkeypatch.setattr(stringtie, "split_long_introns", lambda *a, **k: 0)

    stringtie.run_stringtie(
        _config(
            tmp_path,
            stringtie_min_coverage=3.0,
            stringtie_min_isoform_fraction=0.25,
            stringtie_min_junction_coverage=3.0,
        )
    )

    (cmd,) = cmds
    assert cmd[cmd.index("-c") + 1] == "3.0"
    assert cmd[cmd.index("-f") + 1] == "0.25"
    assert cmd[cmd.index("-j") + 1] == "3.0"


def test_run_stringtie_star_bam_not_bounded(tmp_path, monkeypatch):
    cmds: list[list[str]] = []
    monkeypatch.setattr(stringtie, "run_cmd", lambda cmd, **kw: cmds.append(cmd))
    called: list[int] = []
    monkeypatch.setattr(stringtie, "split_long_introns", lambda *a, **k: called.append(1))

    stringtie.run_stringtie(_config(tmp_path, aligner="star"))

    assert called == []  # STAR BAM is already --alignIntronMax-bounded
    (cmd,) = cmds
    assert cmd[1] == str(tmp_path / "STAR_Aligned.sortedByCoord.out.bam")


def test_run_stringtie_skips_when_output_exists(tmp_path, monkeypatch):
    (tmp_path / "stringtie.gtf").write_text("# already done\n")
    cmds: list[list[str]] = []
    monkeypatch.setattr(stringtie, "run_cmd", lambda cmd, **kw: cmds.append(cmd))
    monkeypatch.setattr(stringtie, "split_long_introns", lambda *a, **k: 0)

    stringtie.run_stringtie(_config(tmp_path))

    assert cmds == []
