"""Shared utility functions for the eukan pipeline."""

from __future__ import annotations

import os
import shutil
from collections.abc import Iterable
from pathlib import Path


def symlink(target: Path, link: Path) -> None:
    """Create a symlink, replacing any existing one."""
    link.unlink(missing_ok=True)
    os.symlink(target, link)


def concat_files(srcs: Iterable[Path], dest: Path) -> None:
    """Stream-concatenate ``srcs`` into ``dest`` (binary, no whole-file reads).

    Used to merge GFF3 / FASTA / hint files without loading large inputs
    into memory.  Output mode is binary so the call works regardless of
    line endings or encoding declared by the source files.
    """
    with open(dest, "wb") as outfile:
        for src in srcs:
            with open(src, "rb") as infile:
                shutil.copyfileobj(infile, outfile)
