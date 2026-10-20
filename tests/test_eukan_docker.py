"""Tests for the ``eukan-docker`` wrapper's ``docker run`` argument list.

A stub ``docker`` first on PATH prints the argv the wrapper would run.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

WRAPPER = Path(__file__).resolve().parents[1] / "eukan-docker"
KEY_MOUNT = "type=bind,source={home}/.gm_key,target=/home/eukan/.gm_key,readonly"


def _docker_argv(tmp_path: Path, home: Path) -> list[str]:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    stub = bin_dir / "docker"
    stub.write_text('#!/bin/sh\nprintf "%s\\n" "$@"\n')
    stub.chmod(0o755)
    out = subprocess.run(
        ["bash", str(WRAPPER), "annotate", "-g", "genome.fasta"],
        cwd=tmp_path, stdin=subprocess.DEVNULL, capture_output=True, text=True,
        env={"PATH": f"{bin_dir}:/usr/bin:/bin", "HOME": str(home)}, check=True,
    )
    return out.stdout.splitlines()


def test_wrapper_parses() -> None:
    subprocess.run(["bash", "-n", str(WRAPPER)], check=True)


def test_host_gm_key_is_mounted_read_only(tmp_path: Path) -> None:
    home = tmp_path / "home"
    home.mkdir()
    (home / ".gm_key").write_text("key\n")
    argv = _docker_argv(tmp_path, home)
    mount = KEY_MOUNT.format(home=home)
    assert mount in argv
    assert argv[argv.index(mount) - 1] == "--mount"
    # The mount is a docker option, so it precedes the image name.
    assert argv.index(mount) < argv.index("eukan")
    assert argv[-3:] == ["annotate", "-g", "genome.fasta"]


def test_no_host_gm_key_keeps_the_baked_one(tmp_path: Path) -> None:
    home = tmp_path / "home"
    home.mkdir()
    argv = _docker_argv(tmp_path, home)
    assert not any(".gm_key" in a for a in argv)
    assert argv[-4:] == ["eukan", "annotate", "-g", "genome.fasta"]
