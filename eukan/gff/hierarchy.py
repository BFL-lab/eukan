"""Feature hierarchy completion and final prettification (with CDS phases)."""

from __future__ import annotations

from collections.abc import Iterator

import gffutils

from eukan.gff.transforms import derived_feature, swap_id_kind

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

    # exon but no CDS → derive CDS from exons, phased from each mRNA's 5'-most
    # exon (SNAP's exon-only models start with their start codon)
    if "CDS" not in feats and "exon" in feats:
        for mRNA in gff3.features_of_type("mRNA"):
            phase = 0
            for exon in gff3.children(
                mRNA, featuretype="exon", order_by="start", reverse=mRNA.strand == "-",
            ):
                attrs = dict(exon.attributes)
                attrs["ID"] = [swap_id_kind(exon.id, "exon", "CDS")]
                yield derived_feature(exon, "CDS", attrs, frame=str(phase))
                phase = (3 - ((exon.end - exon.start + 1 - phase) % 3)) % 3

    # CDS but no exon → derive exons from CDS
    if "exon" not in feats and "CDS" in feats:
        for cds in gff3.features_of_type("CDS"):
            cds_id = cds.attributes["ID"][0]
            attrs = dict(cds.attributes)
            attrs["ID"] = [swap_id_kind(cds_id, r"CDS|cds", "exon")]
            yield derived_feature(cds, "exon", attrs)


# ---------------------------------------------------------------------------
# Final GFF3 prettification
# ---------------------------------------------------------------------------


def _cds_phases(gff3: gffutils.FeatureDB, mRNA: gffutils.Feature) -> dict[str, str]:
    """Phase of each CDS of *mRNA*, walking from the 5'-most CDS (phase 0)."""
    phases: dict[str, str] = {}
    phase = 0
    for cds in gff3.children(
        mRNA, featuretype="CDS", order_by="start", reverse=mRNA.strand == "-",
    ):
        phases[cds.id] = str(phase)
        phase = (3 - ((cds.end - cds.start + 1 - phase) % 3)) % 3
    return phases


def prettify_gff3(
    gff3: gffutils.FeatureDB, shortname: str
) -> Iterator[gffutils.Feature]:
    """Assign clean sequential locus tags to all features.

    Yields each gene, then each of its mRNAs followed by that mRNA's exons and
    then its CDSs (other child types, such as UTRs, are left out), with every
    CDS phase recomputed from the mRNA's 5'-most CDS.
    """
    for gene_num, gene in enumerate(
        gff3.features_of_type("gene", order_by=("seqid", "start")), start=1
    ):
        gene_name = f"{shortname}_{gene_num:05d}"
        gene.attributes["ID"] = gene_name
        gene.attributes["locus_tag"] = gene_name
        gene.attributes["Name"] = gene_name
        gene.source = "eukannotpass"
        yield gene

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

            phases = _cds_phases(gff3, mRNA)
            for child_type in ("exon", "CDS"):
                for f_num, f in enumerate(
                    gff3.children(mRNA, featuretype=child_type), start=1
                ):
                    if child_type == "CDS":
                        f.frame = phases[f.id]
                    f_name = f"{mRNA_name}:{f.featuretype}:{f_num}"
                    f.attributes.update(
                        ID=f_name, Parent=mRNA_name, locus_tag=f_name
                    )
                    f.source = "eukannotpass"
                    yield f
