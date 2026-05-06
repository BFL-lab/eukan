"""Cross-pipeline artifact registry.

A single source of truth for files that cross pipeline boundaries —
written by one pipeline, read by another. Centralizing the filenames
means renaming an artifact is one edit instead of grepping the codebase.

Two flavours of artifact:

* **Static** — fixed filename, lives directly under ``work_dir``. Use
  the :class:`Artifact` enum and :func:`find`.
* **Dynamic** — filename derived from the genome stem (e.g.
  ``<stem>.masked.fasta``). Helper functions live next to the enum.
"""

from __future__ import annotations

from enum import StrEnum
from pathlib import Path


class Artifact(StrEnum):
    """Static cross-pipeline artifacts (filename = enum value)."""

    # --- assembly outputs (consumed by annotation auto-discovery) ---
    NR_TRANSCRIPTS_FASTA = "nr_transcripts.fasta"
    NR_TRANSCRIPTS_GFF   = "nr_transcripts.gff3"
    RNASEQ_HINTS         = "hints_rnaseq.gff"

    # --- assembly diagnostics consumed by AUGUSTUS ---
    SPLICE_SUMMARY = "splice_site_summary.json"

    # --- repeats outputs consumed by AUGUSTUS ---
    REPEATMASK_HINTS = "hints_repeatmask.gff"

    # --- annotation outputs (consumed by func-annot and prep-submission) ---
    FINAL_GFF3      = "final.gff3"
    FINAL_FUNC_GFF3 = "final.mod.gff3"


def find(work_dir: Path, artifact: Artifact) -> Path | None:
    """Resolve an artifact under *work_dir*, returning ``None`` if absent."""
    path = work_dir / artifact.value
    return path if path.exists() else None


def masked_genome(work_dir: Path, stem: str) -> Path:
    """Path to the softmasked genome produced by ``eukan mask-repeats``.

    Filename pattern is ``<stem>.masked.fasta``; the file may or may not
    exist on disk. Use :func:`pathlib.Path.exists` if you need to check.
    """
    return work_dir / f"{stem}.masked.fasta"
