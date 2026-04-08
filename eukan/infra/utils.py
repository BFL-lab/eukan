"""Shared utility functions for the eukan pipeline."""

from __future__ import annotations

import os
from pathlib import Path


def symlink(target: Path, link: Path) -> None:
    """Create a symlink, replacing any existing one."""
    link.unlink(missing_ok=True)
    os.symlink(target, link)


def step_done(wd: Path, marker_files: list[str]) -> bool:
    """Check if a step already completed by verifying expected output files exist.

    This checks for the presence of multiple marker files within a working
    directory.  For single-output step caching used by the annotation pipeline,
    see :func:`eukan.infra.steps.step_complete`.
    """
    return all((wd / f).exists() for f in marker_files)
