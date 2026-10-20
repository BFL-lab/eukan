"""conda-activate.sh exports every env var tools.toml sets under conda.

The eukan CLI sets these itself (eukan.infra.environ); the script is for
running the tools directly, so it must not drift from tools.toml.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

from eukan.infra.tools_registry import load_tools

SCRIPT = Path(__file__).resolve().parents[1] / "conda-activate.sh"


def test_exports_every_conda_path_env_var(tmp_path: Path) -> None:
    specs = [s for t in load_tools() for s in t.env_vars if s.path]
    for spec in specs:
        (tmp_path / spec.path).mkdir(parents=True, exist_ok=True)
    names = " ".join(s.var for s in specs)
    out = subprocess.run(
        ["bash", "-c", f'. "$1"; for v in {names}; do printf "%s=%s\\n" "$v" "${{!v}}"; done',
         "bash", str(SCRIPT)],
        env={"CONDA_PREFIX": str(tmp_path), "PATH": "/usr/bin:/bin"},
        capture_output=True, text=True, check=True,
    )
    exported = dict(line.split("=", 1) for line in out.stdout.splitlines())
    for spec in specs:
        assert exported[spec.var].split(":")[0] == str(tmp_path / spec.path), spec.var
