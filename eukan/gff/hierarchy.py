"""Feature hierarchy completion, CDS phase correction, and final prettification."""

from __future__ import annotations

from collections.abc import Iterator

import gffutils

from eukan.gff.transforms import _swap_id_kind, derived_feature

# ---------------------------------------------------------------------------
# Feature hierarchy completion
# ---------------------------------------------------------------------------


def add_missing_feats_to_gff3(
    gff3: gffutils.FeatureDB,
) -> Iterator[gffutils.Feature]:
    """Yield missing mRNA, CDS, or exon features to complete the hierarchy.

    Many gene predictors produce incomplete GFF3 (e.g., gene+CDS but no mRNA,
    or gene+mRNA+exon but no CDS). This generator fills in the gaps.
    """
    feats = set(gff3.featuretypes())

    # gene+CDS but no mRNA → generate mRNA features
    if "mRNA" not in feats and "gene" in feats:
        for gene in gff3.features_of_type("gene"):
            gene_id = gene.attributes["ID"][0]
            yield derived_feature(
                gene, "mRNA",
                {"ID": [f"{gene_id}_mRNA"], "Parent": [gene_id]},
            )
            gene.attributes["ID"] = [f"{gene_id}_gene"]

    # exon but no CDS → derive CDS from exons
    if "CDS" not in feats and "exon" in feats:
        for exon in gff3.features_of_type("exon"):
            phase = (exon.start - 1) % 3 if exon.strand == "-" else (exon.end + 1) % 3
            attrs = dict(exon.attributes)
            attrs["ID"] = [_swap_id_kind(exon.id, "exon", "CDS")]
            yield derived_feature(exon, "CDS", attrs, frame=phase)

    # CDS but no exon → derive exons from CDS
    if "exon" not in feats and "CDS" in feats:
        for cds in gff3.features_of_type("CDS"):
            cds_id = cds.attributes["ID"][0]
            attrs = dict(cds.attributes)
            attrs["ID"] = [_swap_id_kind(cds_id, r"CDS|cds", "exon")]
            yield derived_feature(cds, "exon", attrs)


# ---------------------------------------------------------------------------
# CDS phase correction
# ---------------------------------------------------------------------------


def fix_CDS_phases(
    gff3: gffutils.FeatureDB,
) -> Iterator[gffutils.Feature]:
    """Recalculate CDS phases for all gene models.

    Yields all features (gene, mRNA, exon, CDS) with corrected CDS frame values.
    """
    for gene in gff3.features_of_type("gene"):
        yield gene
        for mRNA in gff3.children(gene, featuretype="mRNA"):
            yield mRNA
            yield from gff3.children(mRNA, featuretype="exon")

            reverse = mRNA.strand == "-"
            all_cds = list(
                gff3.children(
                    mRNA,
                    featuretype="CDS",
                    order_by="start",
                    reverse=reverse,
                )
            )
            if not all_cds:
                continue

            next_phase = 0
            for cds in all_cds:
                phase = 0 if cds is all_cds[0] else next_phase
                cds_len = cds.end - cds.start + 1
                next_phase = (3 - ((cds_len - phase) % 3)) % 3
                cds.frame = str(phase)
                yield cds


# ---------------------------------------------------------------------------
# Final GFF3 prettification
# ---------------------------------------------------------------------------


def prettify_gff3(
    gff3: gffutils.FeatureDB, shortname: str
) -> Iterator[gffutils.Feature]:
    """Assign clean sequential locus tags to all features."""
    for gene_num, gene in enumerate(
        gff3.features_of_type("gene", order_by=("seqid", "start")), start=1
    ):
        gene_name = f"{shortname}_{gene_num:05d}"
        gene.attributes["ID"] = gene_name
        gene.attributes["locus_tag"] = gene_name
        gene.attributes["Name"] = gene_name
        gene.source = "eukannotpass"
        yield gene

        child_types = {
            f.featuretype for f in gff3.children(gene, level=2)
        }

        for mRNA_num, mRNA in enumerate(
            gff3.children(gene, featuretype="mRNA", order_by=("seqid", "start")),
            start=1,
        ):
            mRNA_name = f"{gene_name}.mRNA.{mRNA_num}"
            mRNA.attributes.update(
                ID=mRNA_name, Name=mRNA_name, Parent=gene_name, locus_tag=mRNA_name
            )
            mRNA.source = "eukannotpass"
            yield mRNA

            for child_type in child_types:
                for f_num, f in enumerate(
                    gff3.children(mRNA, featuretype=child_type), start=1
                ):
                    f_name = f"{mRNA_name}:{f.featuretype}:{f_num}"
                    f.attributes.update(
                        ID=f_name, Parent=mRNA_name, locus_tag=f_name
                    )
                    f.source = "eukannotpass"
                    yield f
