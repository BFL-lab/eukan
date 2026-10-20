"""Shared BAM helpers for the mapping steps.

Aligner-neutral utilities used by :mod:`eukan.assembly.minimap2` and
:mod:`eukan.assembly.jaccard`: a quickcheck-based BAM resume guard, an
unmapped-record extractor and a FASTA record counter. Living here (rather than in
an aligner module) keeps them importable without pulling in the aligner. The
mapped-transcript track table lives in :mod:`eukan.assembly.tracks`.
"""

from __future__ import annotations

from pathlib import Path

import pysam

from eukan.exceptions import ExternalToolError
from eukan.infra.logging import get_logger
from eukan.infra.runner import run_cmd

log = get_logger(__name__)


def _bam_is_complete(path: Path) -> bool:
    """True if *path* is a non-empty BAM with a valid BGZF EOF (``quickcheck``).

    Lets a re-run reuse a BAM from a previous attempt that failed only after it
    was written (e.g. while mapping the other track), rather than re-mapping (the
    expensive step). A truncated BAM from an interrupted or killed mapper or sort
    fails quickcheck and is re-mapped.
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

