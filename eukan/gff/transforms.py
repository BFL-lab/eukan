"""Generic GFF3/GTF transform callbacks and feature-construction helpers.

Transform callbacks here are tool-agnostic — anything specific to a
particular gene predictor lives next to that tool's wrapper in
``annotation/<tool>.py``. The two shared helpers (``derived_feature``,
``_swap_id_kind``) underlie both :mod:`eukan.gff.hierarchy` and the
ORF coordinate mapping in :mod:`eukan.annotation.orf`.
"""

from __future__ import annotations

import re

import gffutils

# ---------------------------------------------------------------------------
# GTF ↔ GFF3 conversion
# ---------------------------------------------------------------------------


def gtf2gff3(f: gffutils.Feature) -> gffutils.Feature | None:
    """Convert a GTF feature to GFF3 conventions."""
    feat = f.featuretype
    if feat not in ("gene", "mRNA", "transcript", "exon", "CDS"):
        return None
    if feat == "gene":
        f.attributes["ID"] = f.attributes.pop("gene_id", [f.id])
    elif feat in ("transcript", "mRNA"):
        f.attributes["ID"] = f.attributes.pop("transcript_id", [f.id])
        f.attributes["Parent"] = f.attributes.pop("gene_id", [])
        f.featuretype = "mRNA"
    else:  # exon or CDS
        f.attributes["Parent"] = f.attributes.pop("transcript_id", [])
        f.attributes.pop("gene_id", None)
        f.attributes["ID"] = [f.id]
    return f


def gff3_it(f: gffutils.Feature) -> gffutils.Feature:
    """Normalize transcript/gene features for GFF3 output."""
    if f.featuretype == "transcript":
        f.featuretype = "mRNA"
        f.attributes["ID"] = f.attributes["Parent"][0] + "_mRNA"
    elif f.featuretype == "gene":
        parent = f.attributes.pop("Parent", None)
        if parent:
            f.attributes["ID"] = parent
    else:
        f.attributes["Parent"] = f.attributes["Parent"][0] + "_mRNA"
    if not f.attributes["ID"] and f.id:
        f.attributes["ID"] = f.id
    return f


def fix_contig_names(f: gffutils.Feature) -> gffutils.Feature:
    """Strip trailing description from contig names (GeneMark quirk)."""
    f.chrom = f.chrom.split(" ")[0]
    return f


def prot2augustus_hints(f: gffutils.Feature) -> gffutils.Feature | None:
    """Convert protein alignment CDS features to AUGUSTUS CDSpart hints."""
    if f.featuretype == "CDS":
        group = f.attributes["Parent"]
        f.attributes = {"pri": "1", "src": "P", "group": group}
        f.featuretype = "CDSpart"
        return f
    return None


# ---------------------------------------------------------------------------
# Feature construction helpers
# ---------------------------------------------------------------------------


def derived_feature(
    parent: gffutils.Feature,
    featuretype: str,
    attributes: dict,
    *,
    frame: str | int = ".",
    start: int | None = None,
    end: int | None = None,
) -> gffutils.Feature:
    """Build a Feature inheriting seqid/source/strand from a parent feature.

    Coordinates default to the parent's; pass ``start``/``end`` to override
    when projecting (e.g. ORF coords onto an exon).
    """
    return gffutils.Feature(
        seqid=parent.chrom,
        source=parent.source,
        featuretype=featuretype,
        start=parent.start if start is None else start,
        end=parent.end if end is None else end,
        strand=parent.strand,
        frame=frame,
        attributes=attributes,
    )


def _swap_id_kind(feature_id: str, from_kind: str, to_kind: str) -> str:
    """Replace a feature-kind token in an ID, falling back to a suffix.

    e.g. ``exon1`` → ``CDS1``; ``g1`` (no match) → ``g1:CDS``.
    """
    new_id = re.sub(from_kind, to_kind, feature_id)
    if new_id == feature_id:
        new_id = f"{feature_id}:{to_kind}"
    return new_id
