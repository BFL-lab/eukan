"""The mapped transcript tracks and their per-track filenames.

Trinity contributes two transcript sets — de novo and genome-guided — and both
are transcript-coordinate FASTAs mapped back to the genome (genome-guided only
clusters reads per locus; it is not genome-native like StringTie was). So every
downstream evidence step (strand_correct → defuse → sl_cut → combinr) processes
the two tracks uniformly, replacing the old "one de novo set + one StringTie
GFF3" special-casing.

``segemehl._TRANSCRIPT_SETS`` is the single source of truth for which tracks
exist; everything here is derived from it so adding/removing a track is a
one-line change there.
"""

from __future__ import annotations

from pathlib import Path

from eukan.assembly.segemehl import _TRANSCRIPT_SETS

# Per-track model-file variants, latest-wins (most-processed first). Each is
# ``<stem>{suffix}`` where ``<stem>`` is e.g. ``trinity-denovo.genome``:
#   .defuse.gff3   homology de-fusion (defuse.py, only with --defuse)
#   .stranded.gff3 homology strand correction (strand_correction.py)
#   .gff3          the raw BAM->models conversion (strand_correction.py step 1)
_MODEL_VARIANTS = (".defuse.gff3", ".stranded.gff3", ".gff3")


def mapped_transcript_stems() -> tuple[str, ...]:
    """The ``<assembler>.genome`` model-file prefixes for the mapped tracks.

    Derived from ``segemehl._TRANSCRIPT_SETS`` — e.g.
    ``("trinity-denovo.genome", "trinity-gg.genome")``. Append a variant suffix
    (``.gff3``, ``.stranded.gff3``, ``.defuse.gff3``, ``.maxintron.gff3``,
    ``.cut.gff3``) to name that track's file.
    """
    # "trinity-denovo.genome.bam" -> "trinity-denovo.genome"
    return tuple(bam[: -len(".bam")] for _, bam in _TRANSCRIPT_SETS)


def resolve_model_source(wd: Path, stem: str) -> Path | None:
    """The latest model variant present for *stem*: defuse > stranded > raw GFF3.

    Returns ``None`` when the track produced no models (e.g. that Trinity mode
    found nothing, so its ``.genome.gff3`` was never written).
    """
    for suffix in _MODEL_VARIANTS:
        p = wd / f"{stem}{suffix}"
        if p.exists():
            return p
    return None
