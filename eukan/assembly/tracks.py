"""The mapped transcript tracks and their per-track filenames.

Trinity contributes two transcript sets — de novo and genome-guided — and both
are transcript-coordinate FASTAs mapped back to the genome (genome-guided only
clusters reads per locus; it is not genome-native like StringTie was). So every
downstream evidence step (strand_correct → defuse → sl_cut → combinr) processes
the two tracks uniformly, replacing the old "one de novo set + one StringTie
GFF3" special-casing.

``TRANSCRIPT_SETS`` is the single source of truth for which tracks exist; the
per-track filenames the assembly steps use are derived from it, so adding or
removing a track is a one-line change here.
"""

from __future__ import annotations

from pathlib import Path

# Trinity transcripts are mapped to the genome SPLICED (minimap2 -x splice:hq).
# Both Trinity modes emit transcript-coordinate FASTAs (genome-guided only
# clusters reads per locus), so both are mapped: one (query FASTA, output BAM)
# pair per track.
TRANSCRIPT_SETS: tuple[tuple[str, str], ...] = (
    ("trinity-denovo.fasta", "trinity-denovo.genome.bam"),
    ("trinity-gg.fasta", "trinity-gg.genome.bam"),
)
GENOME_BAM_SUFFIX = ".genome.bam"

# Per-track model files are ``<stem>{suffix}``, where ``<stem>`` is a mapped stem
# such as ``trinity-denovo.genome`` (see mapped_transcript_stems()).
RAW_GFF3 = ".gff3"  # the BAM->models conversion (strand_correction.py step 1)
STRANDED_GFF3 = ".stranded.gff3"  # homology strand correction (strand_correction.py)
DEFUSE_GFF3 = ".defuse.gff3"  # homology de-fusion (defuse.py, only with --defuse)
MAXINTRON_GFF3 = ".maxintron.gff3"  # max-intron split (max_intron.py)
CUT_GFF3 = ".cut.gff3"  # SL cut (sl_cut.py), consumed by combinr

# The model variants max_intron_split reads, latest-wins (most-processed first).
_MODEL_VARIANTS = (DEFUSE_GFF3, STRANDED_GFF3, RAW_GFF3)


def jaccard_output_name(fasta_name: str) -> str:
    """The ``.jaccard.fasta`` name the jaccard step writes for *fasta_name*."""
    return fasta_name.replace(".fasta", ".jaccard.fasta")


def resolve_query(wd: Path, query_name: str) -> Path:
    """The transcript FASTA to map: the jaccard-clipped sibling if it exists.

    The jaccard step (``assembly/jaccard.py``) rewrites each de novo / genome-
    guided FASTA into a ``<stem>.jaccard.fasta``; prefer it so fused contigs are
    split before consolidation. Falls back to the original when clipping is off.
    """
    clipped = wd / jaccard_output_name(query_name)
    return clipped if clipped.exists() and clipped.stat().st_size > 0 else wd / query_name


def mapped_transcript_stems() -> tuple[str, ...]:
    """The ``<assembler>.genome`` model-file prefixes for the mapped tracks.

    Derived from ``TRANSCRIPT_SETS`` — e.g.
    ``("trinity-denovo.genome", "trinity-gg.genome")``. Append a suffix above
    (``RAW_GFF3``, ``STRANDED_GFF3``, ``DEFUSE_GFF3``, ``MAXINTRON_GFF3``,
    ``CUT_GFF3``) to name that track's file.
    """
    # "trinity-denovo.genome.bam" -> "trinity-denovo.genome"
    return tuple(bam[: -len(".bam")] for _, bam in TRANSCRIPT_SETS)


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
