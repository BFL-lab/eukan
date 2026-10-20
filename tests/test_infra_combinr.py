"""Tests for eukan.infra.combinr — binary resolution and the version gate."""

from __future__ import annotations

import pytest

from eukan.exceptions import ConfigurationError
from eukan.infra.combinr import combinr_bin, run_combinr
from eukan.settings import AssemblyConfig, PipelineConfig


def _fake_combinr(tmp_path, banner: str):
    """An executable printing *banner* for --version; other calls echo argv to stdout."""
    exe = tmp_path / "bin" / "combinr"
    exe.parent.mkdir()
    exe.write_text(
        f'#!/bin/sh\n[ "$1" = --version ] && {{ echo "{banner}"; exit 0; }}\necho "$@"\n'
    )
    exe.chmod(0o755)
    return exe


def _config(tmp_path, **kw):
    return AssemblyConfig(genome=tmp_path / "g.fa", work_dir=tmp_path, num_cpu=1, **kw)


def test_combinr_bin_defaults_to_path_lookup(tmp_path):
    assert combinr_bin(_config(tmp_path)) == "combinr"
    cfg = PipelineConfig(
        genome=tmp_path / "g.fa", proteins=[tmp_path / "p.faa"], work_dir=tmp_path,
        num_cpu=1, combinr_path=tmp_path / "bin" / "combinr",
    )
    assert combinr_bin(cfg) == str(tmp_path / "bin" / "combinr")


@pytest.mark.parametrize("banner", ["combinr 0.3.0", "combinr 0.3.1", "combinr 0.10.0"])
def test_supported_release_runs_the_command(tmp_path, banner):
    cfg = _config(tmp_path, combinr_path=_fake_combinr(tmp_path, banner))
    run_combinr(cfg, ["assemble", "-i", "x.gff3"], cwd=tmp_path, out_file="out.txt")
    assert (tmp_path / "out.txt").read_text() == "assemble -i x.gff3\n"
    assert not (tmp_path / "combinr.version").exists()  # the probe leaves nothing behind


@pytest.mark.parametrize("banner", ["combinr 0.1.1", "combinr 0.2.9", "combinr dev build"])
def test_old_or_unknown_release_is_refused(tmp_path, banner):
    cfg = _config(tmp_path, combinr_path=_fake_combinr(tmp_path, banner))
    with pytest.raises(ConfigurationError, match=r"combinr >= 0\.3\.0") as exc:
        run_combinr(cfg, ["assemble"], cwd=tmp_path, out_file="out.txt")
    assert banner in str(exc.value)
    assert not (tmp_path / "out.txt").exists()
