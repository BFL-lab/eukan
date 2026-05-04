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


def find_resource(relpath: str) -> Path | None:
    """Locate a file shipped under the repo root (e.g. ``configs/x.cfg``).

    Searches, in order: ``$EUKAN_ROOT/<relpath>`` (set by the Dockerfile),
    ``<package_parent>/<relpath>`` (working from a source checkout), and
    ``<cwd>/<relpath>``.  The package-parent fallback is correct only when
    eukan is run from its source tree — once installed via ``pip install
    .`` the package lives in site-packages and has no ``configs/`` sibling,
    which is why the env var is the primary lookup.
    """
    candidates: list[Path] = []
    eukan_root = os.environ.get("EUKAN_ROOT")
    if eukan_root:
        candidates.append(Path(eukan_root) / relpath)
    candidates.append(Path(__file__).resolve().parent.parent.parent / relpath)
    candidates.append(Path.cwd() / relpath)
    return next((c for c in candidates if c.exists()), None)


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
