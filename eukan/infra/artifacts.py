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
from eukan.infra.logging import get_logger
from eukan.infra.manifest import RunManifest, StepStatus

log = get_logger(__name__)


class Artifact(StrEnum):
    """Static cross-pipeline artifacts (filename = enum value)."""

    # --- assembly outputs (consumed by annotation auto-discovery) ---
    NR_TRANSCRIPTS_FASTA = "nr_transcripts.fasta"
    NR_TRANSCRIPTS_GFF   = "nr_transcripts.gff3"
    RNASEQ_HINTS         = "hints_rnaseq.gff"

    # --- assembly diagnostics consumed by AUGUSTUS ---
    SPLICE_SUMMARY = "splice_site_summary.json"

    # --- assembly diagnostic emitted by the soft-clip / intron walk ---
    SOFTCLIP_DIAGNOSTIC = "softclip_diagnostic_summary.json"

    # --- assembly poly-A / poly-T soft-clip + transcript characterization ---
    POLYA_DIAGNOSTIC = "polyA_diagnostic.json"

    # --- SL trans-splice acceptor sites (persist-only; not auto-discovered) ---
    SL_ACCEPTORS = "sl_acceptors.gff3"

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
    Artifact.SOFTCLIP_DIAGNOSTIC:  "assemble",
    Artifact.POLYA_DIAGNOSTIC:     "assemble",
    Artifact.SL_ACCEPTORS:         "assemble",
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


# Producer label (the run-dir subdir in _PRODUCER) -> manifest step-key prefix.
_PRODUCER_MANIFEST_PREFIX: dict[str, str] = {
    "assemble":     "assembly",
    "mask-repeats": "repeats",
    "annotate":     "annotation",
    "func-annot":   "functional",
}


def _producer_completed(manifest: RunManifest, producer: str) -> bool:
    """True if any step of *producer*'s pipeline is recorded ``completed``."""
    prefix = _PRODUCER_MANIFEST_PREFIX.get(producer)
    if prefix is None:
        return False
    return any(
        key.startswith(f"{prefix}/") and rec.status == StepStatus.completed
        for key, rec in manifest.steps.items()
    )


def find_or_warn(
    work_dir: Path, artifact: Artifact, manifest: RunManifest | None,
) -> Path | None:
    """Like :func:`find`, but warn when a *previously produced* artifact is gone.

    Optional cross-pipeline artifacts (RepeatMasker hints, the STAR splice
    summary, ...) are produced by one pipeline and consumed by another, which
    cannot rebuild them. :func:`find` returns ``None`` for a missing one and the
    consumer falls back silently. This variant additionally logs a WARNING when
    the producing step is recorded *completed* in the manifest yet the artifact
    is missing — genuine loss — while staying silent when the producer never ran
    (the artifact is legitimately absent). ``manifest=None`` disables the warning.
    """
    path = find(work_dir, artifact)
    if path is not None:
        return path
    producer = _PRODUCER.get(artifact)
    if manifest is not None and producer is not None and _producer_completed(manifest, producer):
        log.warning(
            "%s was produced by a completed '%s' run but is now missing; "
            "proceeding without it.",
            artifact.value, producer,
        )
    return None


def unmapped_transcripts(work_dir: Path, assembler: str) -> Path:
    """Per-assembler unmapped-transcript FASTA from the map_transcripts step.

    Filename pattern is ``<assembler>.unmapped_transcripts.fasta`` (e.g.
    ``rnaspades.unmapped_transcripts.fasta``); written for reference only.
    """
    return work_dir / f"{assembler}.unmapped_transcripts.fasta"


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
