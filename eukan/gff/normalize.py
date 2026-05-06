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
from eukan.gff import parser as gffparser
from eukan.gff.io import featuredb2gff3_file


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
        raw: Path to a GFF/GTF file or an existing FeatureDB.
        out: Destination GFF3 path.
        parse_transform: Optional gffutils ``transform=`` callback applied
            on the initial parse — typically a tool-specific fixup
            (``Spaln.fix_cds_featuretype``, ``Augustus.clean``, etc.).
        post_transform: Optional second transform pass applied after the
            missing-feature fill-in. Used when a fix needs to see the
            completed hierarchy (e.g. ``Spaln.fix_ids``).
        add_missing_features: When True (the default), runs
            :func:`eukan.gff.parser.add_missing_feats_to_gff3` to back-fill
            mRNA/CDS/exon features from a partial output.
        fix_contig_names: When True, applies
            :func:`eukan.gff.parser.fix_contig_names` to strip whitespace
            descriptions from chrom names. Use for tools (e.g. GeneMark)
            that copy the full FASTA header into the GFF seqid column.

    Returns:
        ``out``, after writing the normalized GFF3.
    """
    db = create_gff_db(raw, transform=parse_transform)

    if add_missing_features:
        db.update(
            gffparser.add_missing_feats_to_gff3(db),
            merge_strategy="create_unique",
        )

    if post_transform is not None:
        db = transform_db(db, post_transform)

    if fix_contig_names:
        db = transform_db(db, gffparser.fix_contig_names)

    featuredb2gff3_file(db, out)
    return out
