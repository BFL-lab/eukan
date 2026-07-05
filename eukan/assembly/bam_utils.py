"""Shared BAM helpers and the mapped-transcript track definitions.

Aligner-neutral utilities used by :mod:`eukan.assembly.minimap2`,
:mod:`eukan.assembly.jaccard`, and :mod:`eukan.assembly.tracks`: the
``(query FASTA, output BAM)`` track table that is the single source of truth for
the mapped-transcript stems, the jaccard-sibling query resolver, a
quickcheck-based BAM resume guard, an unmapped-record extractor, and the
coordinate-sort/filter recipe with its full-disk error translation. Living here
(rather than in an aligner module) keeps them importable without pulling in the
aligner and avoids an import cycle with ``tracks``.
"""

from __future__ import annotations

import shutil
from pathlib import Path

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

# ---------------------------------------------------------------------------
# Transcript -> genome mapping tracks (consumed by minimap2.map_transcripts_*)
# ---------------------------------------------------------------------------
# Trinity transcripts are mapped to the genome SPLICED (minimap2 -x splice:hq).
# Both Trinity modes emit transcript-coordinate FASTAs (genome-guided only
# clusters reads per locus), so both are mapped here. These (query FASTA, output
# BAM) pairs and the jaccard-sibling resolver are the single source of truth for
# the mapped-transcript stems (see :mod:`eukan.assembly.tracks`).
_TRANSCRIPT_SETS: tuple[tuple[str, str], ...] = (
    ("trinity-denovo.fasta", "trinity-denovo.genome.bam"),
    ("trinity-gg.fasta", "trinity-gg.genome.bam"),
)
_GENOME_BAM_SUFFIX = ".genome.bam"


def _is_gzipped(path: Path) -> bool:
    """Check if a file is gzip-compressed by reading the magic bytes."""
    try:
        with open(path, "rb") as f:
            return f.read(2) == b"\x1f\x8b"
    except OSError:
        return False


def _resolve_query(wd: Path, query_name: str) -> Path:
    """The transcript FASTA to map: the jaccard-clipped sibling if it exists.

    The jaccard step (``assembly/jaccard.py``) rewrites each de novo / genome-
    guided FASTA into a ``<stem>.jaccard.fasta``; prefer it so fused contigs are
    split before consolidation. Falls back to the original when clipping is off.
    """
    clipped = wd / query_name.replace(".fasta", ".jaccard.fasta")
    return clipped if clipped.exists() and clipped.stat().st_size > 0 else wd / query_name


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
    import pysam

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
