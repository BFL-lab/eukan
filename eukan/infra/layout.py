"""Per-step output directory layout.

Each ``eukan`` subcommand writes into its own subdirectory under a
shared "run dir" so a multi-step run stays organized:

    run-dir/
    ├── repeats/        eukan mask-repeats
    ├── assemble/       eukan assemble
    ├── annotate/       eukan annotate
    ├── func-annot/     eukan func-annot
    └── submission/     eukan prep-submission

The CLI for each subcommand resolves its own ``work_dir`` to
``<cwd>/<step-subdir>/``. Cross-pipeline artifact lookups
(:mod:`eukan.infra.artifacts`) then fall back from the step's own
``work_dir`` to the sibling step that produced the artifact, so a
downstream subcommand finds upstream outputs without explicit paths.
"""

from __future__ import annotations

from pathlib import Path

# Step name (matches CLI subcommand) → subdirectory under the run dir.
PIPELINE_SUBDIRS: dict[str, str] = {
    "annotate":        "annotate",
    "assemble":        "assemble",
    "mask-repeats":    "repeats",
    "func-annot":      "func-annot",
    "prep-submission": "submission",
}


def step_work_dir(step: str, base: Path | None = None) -> Path:
    """Return the per-step work_dir under *base* (defaults to cwd).

    Raises ``KeyError`` if *step* is not a registered pipeline step.
    """
    if base is None:
        base = Path.cwd()
    return base / PIPELINE_SUBDIRS[step]


def sibling_step_dir(work_dir: Path, step: str) -> Path:
    """Return ``<work_dir.parent>/<step-subdir>`` (the sibling step dir).

    Used by auto-discovery: given the current step's ``work_dir``,
    locate where another step would have written its outputs.
    """
    return work_dir.parent / PIPELINE_SUBDIRS[step]
