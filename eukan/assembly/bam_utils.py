"""Shared BAM helpers for the mapping steps.

Aligner-neutral utilities used by :mod:`eukan.assembly.minimap2` and
:mod:`eukan.assembly.jaccard`: a quickcheck-based BAM resume guard, an
unmapped-record extractor, a FASTA record counter, and the coordinate-sort/filter
recipe with its full-disk error translation. Living here (rather than in an
aligner module) keeps them importable without pulling in the aligner. The
mapped-transcript track table lives in :mod:`eukan.assembly.tracks`.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pysam

from eukan.exceptions import ExternalToolError
from eukan.infra.logging import get_logger
from eukan.infra.runner import run_cmd, run_piped

log = get_logger(__name__)

# The coordinate-sort needs transient room for its temp spill (~= the sorted
# size) plus the final BAM, both landing on the run-dir filesystem — roughly
# twice the unsorted BAM. Warn before sorting when free space is below this.
_SORT_DISK_HEADROOM = 2
# samtools' multithreaded BGZF writer reports a failed write on a full disk
# with a misleading errno (e.g. "...: Illegal seek"); match these stderr
# fragments so we can translate the cryptic failure into a clear out-of-space
# message instead of letting the user chase a phantom seek bug.
_DISK_FULL_MARKERS = ("illegal seek", "no space", "failed writing", "write failed")


def _bam_is_complete(path: Path) -> bool:
    """True if *path* is a non-empty BAM with a valid BGZF EOF (``quickcheck``).

    Lets a re-run reuse a BAM from a previous attempt that failed only in a
    downstream step (e.g. the sort on a full disk), rather than re-mapping (the
    expensive step). A truncated BAM from an interrupted or killed mapper fails
    quickcheck and is re-mapped.
    """
    if not path.exists() or path.stat().st_size == 0:
        return False
    try:
        run_cmd(["samtools", "quickcheck", str(path)], cwd=path.parent)
    except ExternalToolError:
        return False
    return True


def _write_unmapped_fasta(unsorted_bam: Path, out_fasta: Path) -> int:
    """Extract unmapped records from *unsorted_bam* to *out_fasta* (FASTA).

    The transcript BAM keeps unmapped queries as ``0x4`` records before the
    ``-F 4`` filter drops them; pulling them out here (pysam, no samtools
    dependency) preserves the transcripts that failed to map for inspection.
    Returns the number written.
    """
    n = 0
    bam = pysam.AlignmentFile(str(unsorted_bam), "rb")
    try:
        with open(out_fasta, "w") as fh:
            for read in bam:
                if not read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                seq = read.query_sequence
                if not read.query_name or not seq:
                    continue
                fh.write(f">{read.query_name}\n{seq}\n")
                n += 1
    finally:
        bam.close()
    return n


def count_fasta_records(fasta: Path) -> int:
    """Number of records in *fasta* (lines starting with ``>``)."""
    with open(fasta) as fh:
        return sum(1 for line in fh if line.startswith(">"))


def _coordinate_sort_and_filter(
    unsorted: Path, out_bam: str, wd: Path, num_cpu: int
) -> None:
    """Drop unmapped reads (``-F 4``) and coordinate-sort into *out_bam*.

    Filtering before the sort is cheaper and yields the same result. On a tight
    disk this is the step that fails — samtools needs roughly twice the unsorted
    BAM in transient temp + output space — so we warn up front and, on a
    write/seek failure, re-raise with an explicit out-of-space hint.
    """
    if unsorted.exists():
        src = unsorted.stat().st_size
        free = shutil.disk_usage(wd).free
        if free < _SORT_DISK_HEADROOM * src:
            log.warning(
                "Low disk for coordinate-sort: %.1f GB free vs ~%.1f GB likely "
                "needed (temp spill + output for a %.1f GB BAM). If it fails, "
                "free space and re-run — the unsorted BAM is reused, so the "
                "mapper will not re-run.",
                free / 1e9, _SORT_DISK_HEADROOM * src / 1e9, src / 1e9,
            )
    try:
        run_piped(
            ["samtools", "view", "-bh", "-F", "4", str(unsorted)],
            ["samtools", "sort", "-@", str(num_cpu), "-o", out_bam, "-"],
            cwd=wd,
        )
    except ExternalToolError as exc:
        snippet = (exc.stderr_snippet or "").lower()
        if any(marker in snippet for marker in _DISK_FULL_MARKERS):
            raise ExternalToolError(
                "samtools sort failed writing the coordinate-sorted BAM, which "
                "on a full filesystem samtools reports as 'Illegal seek'. The "
                "sort needs roughly twice the unsorted BAM in transient temp + "
                "output space.",
                tool=exc.tool, returncode=exc.returncode, cmd=exc.cmd,
                stderr_snippet=exc.stderr_snippet,
                hint=(
                    "Free space on the run-dir filesystem (or run the assemble "
                    "step on a larger disk), then re-run: the mapper reuses the "
                    "existing unsorted BAM and skips the slow re-mapping."
                ),
            ) from exc
        raise
