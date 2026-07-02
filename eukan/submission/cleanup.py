"""Pre-emptive GFF3 cleanup before table2asn validation.

Applies five transforms that resolve the bulk of NCBI submission errors
seen on eukan output:

1. Strip UniProt-format metadata (`` OS=...OX=...GN=...PE=...SV=...``)
   from ``product=`` values. Without this, table2asn rewrites every
   product to "hypothetical protein", losing all functional annotation.
2. Strip ``(Fragment)`` from product names (NCBI SUSPECT_PHRASES fatal).
3. Drop mRNAs that have no CDS children (transcript isoforms whose ORFs were
   not called). Drops the parent gene if no CDS-bearing siblings remain.
4. Cap ``inference=`` to at most ``INFERENCE_CAP`` accessions per feature.
5. Move ``Dbxref=KEGG:K…`` to ``Note=KEGG:K…``. KEGG is not in NCBI's
   controlled db_xref list, so leaving it as Dbxref triggers
   SEQ_FEAT.IllegalDbXref on every KOfam-annotated feature; relocating
   to a free-text Note preserves the KO accession in the GenBank record.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import gffutils

from eukan.gff import create_gff_db

INFERENCE_CAP = 3
"""Maximum number of /inference qualifiers per feature.

table2asn caps validation at the BIOSEQ-SET level once the limit is hit
and emits TooManyInferenceAccessions; trimming earlier keeps the
strongest evidence in the record.
"""

DEFAULT_PRODUCT = "hypothetical protein"
"""Fallback when product cleanup leaves nothing usable."""


@dataclass
class CleanupReport:
    """Counts of transforms applied by :func:`clean_gff3_for_submission`."""

    products_cleaned: int = 0
    fragments_stripped: int = 0
    inferences_capped: int = 0
    mrnas_dropped: int = 0
    genes_dropped: int = 0
    kegg_dbxrefs_moved: int = 0

    def summary(self) -> str:
        """One-line human-readable summary for logs."""
        return (
            f"cleaned products: {self.products_cleaned}; "
            f"fragments stripped: {self.fragments_stripped}; "
            f"dropped CDS-less mRNAs: {self.mrnas_dropped} "
            f"(and {self.genes_dropped} orphaned genes); "
            f"capped inferences: {self.inferences_capped}; "
            f"moved KEGG Dbxrefs to Note: {self.kegg_dbxrefs_moved}"
        )


def _clean_product(product: str) -> tuple[str, bool, bool]:
    """Strip UniProt metadata and ``(Fragment)`` from a product name.

    Returns ``(cleaned, os_stripped, fragment_stripped)``. UniProt's
    canonical format is ``Name OS=org OX=tax GN=gene PE=ev SV=ver``,
    always with ``OS=`` first; stripping from `` OS=`` onward removes
    every metadata field in one cut.
    """
    cleaned = product
    os_stripped = False
    fragment_stripped = False

    idx = cleaned.find(" OS=")
    if idx != -1:
        cleaned = cleaned[:idx]
        os_stripped = True

    if "(Fragment)" in cleaned:
        cleaned = cleaned.replace(" (Fragment)", "").replace("(Fragment)", "")
        fragment_stripped = True

    cleaned = cleaned.strip()
    if not cleaned:
        cleaned = DEFAULT_PRODUCT

    return cleaned, os_stripped, fragment_stripped


def _identify_drop_set(db: gffutils.FeatureDB) -> tuple[set[str], set[str]]:
    """Find mRNA IDs without CDS children, and orphan gene IDs.

    A gene is orphaned when none of its mRNAs has a CDS child; in that
    case nothing protein-coding remains, and the locus is dropped.
    """
    mrnas_drop: set[str] = set()
    genes_keep: set[str] = set()

    for mrna in db.features_of_type("mRNA"):
        has_cds = next(iter(db.children(mrna, featuretype="CDS")), None) is not None
        if not has_cds:
            mrnas_drop.add(mrna.id)
            continue
        for parent in db.parents(mrna, featuretype="gene"):
            genes_keep.add(parent.id)

    all_genes = {g.id for g in db.features_of_type("gene")}
    genes_drop = all_genes - genes_keep
    return mrnas_drop, genes_drop


def _move_kegg_dbxref_to_note(f: gffutils.Feature) -> int:
    """Relocate ``Dbxref=KEGG:K…`` entries to ``Note=KEGG:K…``.

    Returns the number of KEGG entries moved. Non-KEGG Dbxrefs are left
    untouched; if all Dbxrefs were KEGG, the attribute is removed
    entirely so empty ``Dbxref=`` doesn't survive into the output.
    """
    if "Dbxref" not in f.attributes:
        return 0

    kegg = [v for v in f.attributes["Dbxref"] if v.startswith("KEGG:")]
    if not kegg:
        return 0

    other = [v for v in f.attributes["Dbxref"] if not v.startswith("KEGG:")]
    if other:
        f.attributes["Dbxref"] = other
    else:
        del f.attributes["Dbxref"]

    existing = list(f.attributes.get("Note", []))
    f.attributes["Note"] = existing + kegg
    return len(kegg)


def _clean_feature(f: gffutils.Feature, report: CleanupReport) -> None:
    """Mutate ``f`` in place: clean product, cap inferences, relocate KEGG. Updates ``report``."""
    if "product" in f.attributes:
        new_products: list[str] = []
        for prod in f.attributes["product"]:
            cleaned, os_strip, frag_strip = _clean_product(prod)
            new_products.append(cleaned)
            if os_strip:
                report.products_cleaned += 1
            if frag_strip:
                report.fragments_stripped += 1
        f.attributes["product"] = new_products

    if "inference" in f.attributes:
        values = f.attributes["inference"]
        if len(values) > INFERENCE_CAP:
            f.attributes["inference"] = values[:INFERENCE_CAP]
            report.inferences_capped += 1

    report.kegg_dbxrefs_moved += _move_kegg_dbxref_to_note(f)


def clean_gff3_for_submission(in_gff: Path, out_gff: Path) -> CleanupReport:
    """Apply submission cleanup transforms and write a cleaned GFF3.

    Identifies the drop set, then walks the gene → mRNA → exon/CDS
    hierarchy and writes only the surviving features with cleaned
    attributes. Returns a :class:`CleanupReport` for the caller to log.
    """
    db = create_gff_db(in_gff)
    mrnas_drop, genes_drop = _identify_drop_set(db)

    report = CleanupReport(
        mrnas_dropped=len(mrnas_drop),
        genes_dropped=len(genes_drop),
    )

    with open(out_gff, "w") as fout:
        fout.write("##gff-version 3\n")
        for gene in db.features_of_type("gene", order_by=("seqid", "start")):
            if gene.id in genes_drop:
                continue
            _clean_feature(gene, report)
            fout.write(f"{gene}\n")
            for mrna in db.children(gene, featuretype="mRNA", order_by="start"):
                if mrna.id in mrnas_drop:
                    continue
                _clean_feature(mrna, report)
                fout.write(f"{mrna}\n")
                for child_type in ("exon", "CDS"):
                    for child in db.children(
                        mrna, featuretype=child_type, order_by="start",
                    ):
                        _clean_feature(child, report)
                        fout.write(f"{child}\n")

    return report
