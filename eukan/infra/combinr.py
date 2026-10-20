"""The external ``combinr`` binary, shared by the two pipelines that run it.

``eukan assemble`` consolidates transcript models with ``combinr assemble``
(:mod:`eukan.assembly.combinr`) and ``eukan annotate`` builds consensus gene
models with ``combinr consensus`` (:mod:`eukan.annotation.combinr_consensus`).
Both resolve the binary the same way and refuse a release older than
:data:`MIN_VERSION`, the one eukan pins: the consensus command line passes
``--evidence-report``, which combinr 0.1.x rejects with a bare clap usage error
("unexpected argument '--evidence-report' found").
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import TYPE_CHECKING

from eukan.exceptions import ConfigurationError
from eukan.infra.runner import run_cmd

if TYPE_CHECKING:
    from eukan.settings import AssemblyConfig, PipelineConfig

MIN_VERSION = (0, 3, 0)


def combinr_bin(config: AssemblyConfig | PipelineConfig) -> str:
    """The combinr executable: explicit ``combinr_path`` or ``combinr`` on PATH."""
    return str(config.combinr_path) if config.combinr_path else "combinr"


def run_combinr(
    config: AssemblyConfig | PipelineConfig, args: list[str], *, cwd: Path, out_file: str,
) -> None:
    """Run ``combinr <args>`` in *cwd* with stdout to *out_file* (within *cwd*).

    Raises:
        ConfigurationError: ``combinr --version`` reports a release older than
            :data:`MIN_VERSION`, or no version at all.
    """
    binary = combinr_bin(config)
    # Same cwd as the real call, so a relative combinr_path resolves the same way.
    banner_file = cwd / "combinr.version"
    run_cmd([binary, "--version"], cwd=cwd, out_file=banner_file.name)
    banner = banner_file.read_text().strip()
    banner_file.unlink()
    found = re.search(r"(\d+)\.(\d+)\.(\d+)", banner)
    if found is None or tuple(map(int, found.groups())) < MIN_VERSION:
        need = ".".join(map(str, MIN_VERSION))
        raise ConfigurationError(
            f"eukan needs combinr >= {need}, but '{binary} --version' reported "
            f"{banner or 'nothing'!r}.",
            hint=(
                f"Install the combinr v{need} release binary from "
                "https://github.com/BFL-lab/combinr/releases (or run "
                "scripts/install-extras.sh) and put it on PATH, or point "
                "--combinr-path at it."
            ),
        )
    run_cmd([binary, *args], cwd=cwd, out_file=out_file)
