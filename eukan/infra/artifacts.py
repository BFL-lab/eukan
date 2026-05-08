"""Cross-pipeline artifact registry.

A single source of truth for files that cross pipeline boundaries —
written by one pipeline, read by another. Centralizing the filenames
means renaming an artifact is one edit instead of grepping the codebase.

Two flavours of artifact:

* **Static** — fixed filename, lives in its producer step's work_dir.
  Use the :class:`Artifact` enum and :func:`find`.
* **Dynamic** — filename derived from the genome stem (e.g.
  ``<stem>.masked.fasta``). Helper functions live next to the enum.

Each artifact knows which step produces it (:data:`_PRODUCER`). Lookups
first check the caller's ``work_dir`` (for flat layouts and direct
in-step references) and then fall back to the producer's sibling step
directory under the run-dir layout — see :mod:`eukan.infra.layout`.
"""

from __future__ import annotations

from enum import StrEnum
from pathlib import Path

from eukan.infra.layout import PIPELINE_SUBDIRS, sibling_step_dir


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


# Maps each artifact to the step that produces it.
_PRODUCER: dict[Artifact, str] = {
    Artifact.NR_TRANSCRIPTS_FASTA: "assemble",
    Artifact.NR_TRANSCRIPTS_GFF:   "assemble",
    Artifact.RNASEQ_HINTS:         "assemble",
    Artifact.SPLICE_SUMMARY:       "assemble",
    Artifact.REPEATMASK_HINTS:     "mask-repeats",
    Artifact.FINAL_GFF3:           "annotate",
    Artifact.FINAL_FUNC_GFF3:      "func-annot",
}


def _candidates(work_dir: Path, filename: str, producer: str | None) -> list[Path]:
    """Return search paths for *filename*: own work_dir first, then producer's sibling."""
    paths = [work_dir / filename]
    if producer and work_dir.name != PIPELINE_SUBDIRS.get(producer):
        paths.append(sibling_step_dir(work_dir, producer) / filename)
    return paths


def find(work_dir: Path, artifact: Artifact) -> Path | None:
    """Resolve an artifact, returning the first existing match or ``None``.

    Searches the caller's ``work_dir`` first, then falls back to the
    producing step's sibling directory under the run-dir layout.
    """
    for path in _candidates(work_dir, artifact.value, _PRODUCER.get(artifact)):
        if path.exists():
            return path
    return None


def masked_genome(work_dir: Path, stem: str) -> Path:
    """Path to the softmasked genome produced by ``eukan mask-repeats``.

    Filename pattern is ``<stem>.masked.fasta``. Returns the first
    existing match (own work_dir, then sibling ``repeats/``); if neither
    exists, returns the in-step path so callers can write there.
    """
    filename = f"{stem}.masked.fasta"
    for path in _candidates(work_dir, filename, "mask-repeats"):
        if path.exists():
            return path
    return work_dir / filename
