"""GFF3 I/O: serialization, counting, and sequence extraction.

Canonical implementations of featuredb2gff3_file, count_gff3_features,
and extract_sequences.
"""

from __future__ import annotations

from collections.abc import Callable, Iterator
from pathlib import Path

import gffutils
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from eukan.gff import create_gff_db
from eukan.gff.rows import iter_gff3_rows
from eukan.infra.genome import ContigIndex


def count_gff3_features(gff3_path: str | Path, feature_type: str = "gene") -> int:
    """Count features of a given type in a GFF3 file by scanning column 3.

    Fast line-based parsing — does not load the file into a database.
    """
    return sum(1 for cols in iter_gff3_rows(gff3_path) if cols[2] == feature_type)


def featuredb2gff3_file(
    featuredb: gffutils.FeatureDB,
    out: str | Path,
    *,
    header: str = "",
    visit: Callable[[gffutils.Feature], bool] | None = None,
) -> None:
    """Write a FeatureDB to a GFF3 file with proper gene>mRNA>exon/CDS hierarchy.

    *header* is written first. *visit* is called on each feature before it is
    written and may modify it; when it returns False the feature and
    everything under it are left out.
    """
    keep = visit or (lambda feature: True)
    with open(out, "w") as fout:
        fout.write(header)
        for gene in featuredb.features_of_type("gene", order_by=("seqid", "start")):
            if not keep(gene):
                continue
            fout.write(f"{gene}\n")
            for mRNA in featuredb.children(gene, featuretype="mRNA", order_by="start"):
                if not keep(mRNA):
                    continue
                fout.write(f"{mRNA}\n")
                for child_type in ("exon", "CDS"):
                    for f in featuredb.children(
                        mRNA, featuretype=child_type, order_by="start"
                    ):
                        if keep(f):
                            fout.write(f"{f}\n")


def iter_assembled_sequences(
    gff3db: gffutils.FeatureDB,
    fasta: str | Path,
    *,
    child_featuretype: str = "CDS",
) -> Iterator[tuple[gffutils.Feature, str]]:
    """Yield ``(mRNA, assembled_DNA)`` pairs for every mRNA in *gff3db*.

    Concatenates child features (default ``CDS``) ordered by genome start,
    then reverse-complements on the minus strand. The mRNA list is sorted
    by ``(chrom, start)`` so the ContigIndex single-record cache stays warm.
    Skips mRNAs with no matching children. Case is preserved — callers
    that need uppercase (e.g. for ORF regex matching) should ``.upper()``.
    """
    mrnas = sorted(gff3db.features_of_type("mRNA"), key=lambda m: (m.chrom, m.start))
    with ContigIndex(fasta) as contigs:
        for mrna in mrnas:
            seq_parts = [
                str(contigs[child.chrom][child.start - 1 : child.end].seq)
                for child in gff3db.children(
                    mrna, featuretype=child_featuretype, order_by="start",
                )
            ]
            if not seq_parts:
                continue
            seq = "".join(seq_parts)
            if mrna.strand == "-":
                seq = str(Seq(seq).reverse_complement())
            yield mrna, seq


def extract_sequences(
    gff3: str | Path,
    genome: str | Path,
    extract_to: str = "protein",
    genetic_code: int = 1,
) -> Iterator[SeqRecord]:
    """Extract protein or cDNA sequences from a GFF3 + genome.

    Yields SeqRecord objects (caller decides how to write them).
    """
    gff3db = create_gff_db(gff3)
    codon_table = CodonTable.unambiguous_dna_by_id[genetic_code]

    for mrna, dna in iter_assembled_sequences(gff3db, genome, child_featuretype="CDS"):
        seq = Seq(dna)
        if extract_to == "protein":
            seq = seq.translate(table=codon_table)
        yield SeqRecord(seq, id=mrna.id, description="")
