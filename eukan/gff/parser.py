"""GFF3/GTF format conversion and normalization.

Transform callbacks for ``gffutils.create_db(..., transform=fn)`` and
core functions for building complete gene model hierarchies.

Tool-specific transform callbacks are grouped into namespace classes
(:class:`Spaln`, :class:`Gth`, :class:`Snap`, :class:`CodingQuarry`,
:class:`Augustus`) so each tool's quirks live under a single namespace
and the ``gffparser.Tool.fix_thing`` call sites make intent obvious.
"""

from __future__ import annotations

import re
from collections.abc import Callable, Iterator
from pathlib import Path

import gffutils
from Bio import SeqIO

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
            yield gffutils.Feature(
                seqid=gene.chrom,
                source=gene.source,
                featuretype="mRNA",
                start=gene.start,
                end=gene.end,
                strand=gene.strand,
                frame=".",
                attributes={"ID": [f"{gene_id}_mRNA"], "Parent": [gene_id]},
            )
            gene.attributes["ID"] = [f"{gene_id}_gene"]

    # exon but no CDS → derive CDS from exons
    if "CDS" not in feats and "exon" in feats:
        for exon in gff3.features_of_type("exon"):
            phase = (exon.start - 1) % 3 if exon.strand == "-" else (exon.end + 1) % 3
            cds_id = re.sub("exon", "CDS", exon.id)
            if cds_id == exon.id:
                cds_id = f"{exon.id}:CDS"
            attrs = dict(exon.attributes)
            attrs["ID"] = [cds_id]
            yield gffutils.Feature(
                seqid=exon.chrom,
                source=exon.source,
                featuretype="CDS",
                start=exon.start,
                end=exon.end,
                strand=exon.strand,
                frame=phase,
                attributes=attrs,
            )

    # CDS but no exon → derive exons from CDS
    if "exon" not in feats and "CDS" in feats:
        for cds in gff3.features_of_type("CDS"):
            cds_id = cds.attributes["ID"][0]
            exon_id = re.sub(r"CDS|cds", "exon", cds_id)
            if exon_id == cds_id:
                exon_id = f"{cds_id}:exon"
            attrs = dict(cds.attributes)
            attrs["ID"] = [exon_id]
            yield gffutils.Feature(
                seqid=cds.chrom,
                source=cds.source,
                featuretype="exon",
                start=cds.start,
                end=cds.end,
                strand=cds.strand,
                frame=".",
                attributes=attrs,
            )


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


# ---------------------------------------------------------------------------
# SNAP ZFF format conversion
# ---------------------------------------------------------------------------


def gff2zff(gff3: str | Path, fasta: str | Path, output_dir: Path | None = None) -> None:
    """Convert GFF3 to SNAP's ZFF annotation format (``genome.ann``)."""
    from eukan.gff import create_gff_db
    featuredb = create_gff_db(gff3)
    contig_list = [f.id for f in SeqIO.parse(str(fasta), "fasta")]
    contigs: dict[str, list[list]] = {k: [] for k in contig_list}

    for mRNA in featuredb.features_of_type("mRNA"):
        exons = list(featuredb.children(mRNA, featuretype="exon"))
        tot_exons = len(exons)
        parent = mRNA.id

        for i, exon in enumerate(exons, start=1):
            if tot_exons == 1:
                feat = "Esngl"
            elif i == 1:
                feat = "Einit" if mRNA.strand == "+" else "Eterm"
            elif i == tot_exons:
                feat = "Eterm" if mRNA.strand == "+" else "Einit"
            else:
                feat = "Exon"

            start, end = (
                (exon.start, exon.end)
                if mRNA.strand == "+"
                else (exon.end, exon.start)
            )
            contigs[mRNA.chrom].append([feat, start, end, parent])

    out_path = (output_dir / "genome.ann") if output_dir else Path("genome.ann")
    with open(out_path, "w") as fout:
        for contig, entries in contigs.items():
            fout.write(f">{contig}\n")
            for feat, start, end, parent in entries:
                fout.write(f"{feat}  {start} {end} {parent}\n")


# ---------------------------------------------------------------------------
# Tool-specific transform callbacks, grouped by source tool.
#
# Each class is a namespace (no instances) whose staticmethods are
# ``gffutils.create_db(..., transform=fn)`` callbacks. Grouping by tool
# makes it easy to see at a call site (``gffparser.Spaln.fix_ids``)
# which tool's quirks a transform is compensating for.
# ---------------------------------------------------------------------------


class Spaln:
    """Transforms for spaln protein-alignment output."""

    @staticmethod
    def fix_cds_featuretype(f: gffutils.Feature) -> gffutils.Feature:
        """Fix spaln's lowercase ``cds`` featuretype."""
        if f.featuretype == "cds":
            f.attributes["ID"][0] = f"{f.attributes['Parent'][0]}:{f.attributes['ID'][0]}"
            f.featuretype = "CDS"
        return f

    @staticmethod
    def fix_ids(f: gffutils.Feature) -> gffutils.Feature:
        """Normalize spaln feature IDs and source."""
        if f.featuretype in ("CDS", "exon"):
            f.attributes["ID"][0] = f.id
        f.source = "prot_align"
        return f


class Gth:
    """Transforms for GenomeThreader protein-alignment output."""

    @staticmethod
    def filter(f: gffutils.Feature) -> gffutils.Feature | None:
        """Keep only gene and exon features; reparent exons to ``_mRNA``."""
        if f.featuretype in ("gene", "exon"):
            if f.featuretype == "exon":
                f.attributes["Parent"][0] += "_mRNA"
            return f
        return None

    @staticmethod
    def fix_ids(f: gffutils.Feature) -> gffutils.Feature:
        """Normalize GenomeThreader feature IDs."""
        f.source = "prot_align"
        if f.featuretype == "exon":
            f.attributes["ID"] = [f.id]
        if f.featuretype == "mRNA":
            f.id = f.attributes["ID"][0]
        return f


class Snap:
    """Transforms for SNAP gene-prediction output."""

    @staticmethod
    def make_featuretype_transform() -> Callable[[gffutils.Feature], gffutils.Feature]:
        """Return a fresh transform that normalizes SNAP output features.

        Each call creates an independent counter so that multiple
        transform passes don't share state (which would cause ID
        collisions).
        """
        counter: dict[str, int] = {}

        def fix_snap_featuretype(f: gffutils.Feature) -> gffutils.Feature:
            f.source = "snap"
            f.featuretype = "exon"
            # SNAP puts the gene group name as the first attribute key
            parent = next(iter(f.attributes))
            counter[parent] = counter.get(parent, 0) + 1
            exon_id = f"{parent}_exon_{counter[parent]}"
            f.attributes = {"ID": [exon_id], "Parent": [parent]}
            return f

        return fix_snap_featuretype

    @staticmethod
    def homogenize_source(f: gffutils.Feature) -> gffutils.Feature:
        """Set ``source`` to ``snap`` for all features."""
        f.source = "snap"
        return f


class CodingQuarry:
    """Transforms for CodingQuarry gene-prediction output."""

    @staticmethod
    def add_mRNA(f: gffutils.Feature) -> gffutils.Feature:
        """Fix CodingQuarry's incomplete output by attaching mRNA parents."""
        f.source = "codingquarry"
        if f.featuretype == "CDS":
            f.attributes["Parent"][0] += "_mRNA"
        return f

    @staticmethod
    def fix_dup_ids(f: gffutils.Feature) -> gffutils.Feature:
        """Fix duplicate IDs in CodingQuarry output."""
        if f.featuretype == "mRNA":
            f.id = f.attributes["ID"][0]
        else:
            f.attributes["ID"] = [f.id]
        return f


class Augustus:
    """Transforms for AUGUSTUS gene-prediction output."""

    @staticmethod
    def clean(f: gffutils.Feature) -> gffutils.Feature | None:
        """Keep only gene/mRNA/CDS/exon features from AUGUSTUS output."""
        if f.featuretype in ("gene", "mRNA", "CDS", "exon"):
            f.source = "augustus"
            return f
        return None
