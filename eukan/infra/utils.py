"""Shared utility functions for the eukan pipeline."""

from __future__ import annotations

import os
import shutil
from collections.abc import Iterable
from importlib.resources import files
from pathlib import Path


def symlink(target: Path, link: Path) -> None:
    """Create a symlink, replacing any existing one."""
    link.unlink(missing_ok=True)
    os.symlink(target, link)


def package_resource(relpath: str) -> Path | None:
    """Locate a file shipped under ``eukan/data/`` in the installed package.

    Returns the filesystem path or None if the resource doesn't exist.
    Works regardless of install layout (editable, wheel, conda package).
    """
    parts = relpath.split("/")
    res = files("eukan").joinpath("data", *parts)
    p = Path(str(res))
    return p if p.is_file() else None


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
