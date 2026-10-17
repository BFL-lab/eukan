"""Shared utility functions for the eukan pipeline."""

from __future__ import annotations

import hashlib
import os
import shutil
import time
from collections.abc import Iterable
from importlib.resources import files
from pathlib import Path

_MD5_CHUNK_SIZE = 8 * 1024 * 1024  # 8 MiB

# Digests already computed in this process, keyed by (absolute path, size,
# mtime_ns). The kernel stamps mtime from a coarse clock (4 ms ticks here,
# whole seconds on some filesystems), so a same-size rewrite within one tick
# keeps the key: a digest is stored only for a file last modified at least
# this long before it was hashed.
_MD5_CACHE_MIN_AGE_NS = 2 * 10**9
_md5_cache: dict[tuple[str, int, int], str] = {}


def md5_file(path: Path) -> str:
    """Compute MD5 of a file in chunks (never loads entire file into memory).

    Within one process an unchanged file is hashed once: fingerprinting an
    input, recording a step output and verifying it on resume often hash
    the same multi-GB file.
    """
    now = time.time_ns()
    st = os.stat(path)
    key = (os.path.abspath(path), st.st_size, st.st_mtime_ns)
    if key in _md5_cache:
        return _md5_cache[key]
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(_MD5_CHUNK_SIZE), b""):
            h.update(chunk)
    digest = h.hexdigest()
    if now - st.st_mtime_ns >= _MD5_CACHE_MIN_AGE_NS:
        _md5_cache[key] = digest
    return digest


def available_cpus() -> int:
    """Number of CPUs this process may run on.

    That is its affinity mask, which is how SLURM cpusets and
    ``docker --cpuset-cpus`` limit a job, or the machine's count where the
    platform has no affinity call (macOS).
    """
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        return os.cpu_count() or 1


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
