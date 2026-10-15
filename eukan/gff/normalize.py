"""GFF3 normalization pipeline shared by every annotation tool wrapper.

Most annotation tools emit GFF that needs the same post-processing
sequence: parse with a tool-specific transform, fill in any missing
mRNA/CDS/exon features, optionally apply a second transform pass, and
write a clean GFF3 file. This helper centralises that pipeline so each
tool wrapper only contributes the parts that are actually
tool-specific (the transform callbacks).
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import gffutils

from eukan.gff import create_gff_db, transform_db
from eukan.gff.hierarchy import add_missing_feats_to_gff3
from eukan.gff.io import featuredb2gff3_file
from eukan.gff.transforms import fix_contig_names as _fix_contig_names


def normalize_to_gff3(
    raw: str | Path | gffutils.FeatureDB,
    out: Path,
    *,
    parse_transform: Callable | None = None,
    post_transform: Callable | None = None,
    add_missing_features: bool = True,
    fix_contig_names: bool = False,
) -> Path:
    """Run the standard parse → fill-in → post-transform → write pipeline.

    Args:
        raw: Path to a GFF/GTF file or an existing FeatureDB. A FeatureDB
            passed without ``parse_transform`` is used, and updated, in place.
        out: Destination GFF3 path.
        parse_transform: Optional gffutils ``transform=`` callback applied
            on the initial parse — typically a tool-specific fixup
            (``Spaln.fix_cds_featuretype``, ``Augustus.clean``, etc.).
        post_transform: Optional second transform pass applied after the
            missing-feature fill-in. Used when a fix needs to see the
            completed hierarchy (e.g. ``Spaln.fix_ids``).
        add_missing_features: When True (the default), runs
            :func:`eukan.gff.hierarchy.add_missing_feats_to_gff3` to back-fill
            mRNA/CDS/exon features from a partial output.
        fix_contig_names: When True, applies
            :func:`eukan.gff.transforms.fix_contig_names` to strip whitespace
            descriptions from chrom names. Use for tools (e.g. GeneMark)
            that copy the full FASTA header into the GFF seqid column.

    Returns:
        ``out``, after writing the normalized GFF3.
    """
    if isinstance(raw, gffutils.FeatureDB) and parse_transform is None:
        db = raw  # a plain copy pass would rebuild an identical DB
    else:
        db = create_gff_db(raw, transform=parse_transform)

    if add_missing_features:
        db.update(
            add_missing_feats_to_gff3(db),
            merge_strategy="create_unique",
        )

    if fix_contig_names:
        first = post_transform

        def post_and_fix(f: gffutils.Feature) -> gffutils.Feature | None:
            # The seqid is read by no ID or Parent derivation, so stripping it
            # in post_transform's pass equals a second full pass.
            if first is not None:
                f = first(f)
            return _fix_contig_names(f) if f else f

        post_transform = post_and_fix

    if post_transform is not None:
        db = transform_db(db, post_transform)

    featuredb2gff3_file(db, out)
    return out
