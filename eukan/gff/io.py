"""GFF3 I/O: serialization and sequence extraction.

Canonical implementations of featuredb2gff3_file and gff3_to_fasta,
deduplicated from gffparser.py, func-annot, and gff3toseq.
"""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

import gffutils
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from eukan.gff import create_gff_db


def featuredb2gff3_file(featuredb: gffutils.FeatureDB, out: str | Path) -> None:
    """Write a FeatureDB to a GFF3 file with proper gene>mRNA>exon/CDS hierarchy."""
    with open(out, "w") as fout:
        for gene in featuredb.features_of_type("gene", order_by=("seqid", "start")):
            fout.write(f"{gene}\n")
            for mRNA in featuredb.children(gene, featuretype="mRNA", order_by="start"):
                fout.write(f"{mRNA}\n")
                for child_type in ("exon", "CDS"):
                    for f in featuredb.children(
                        mRNA, featuretype=child_type, order_by="start"
                    ):
                        fout.write(f"{f}\n")


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
    contigs = SeqIO.to_dict(SeqIO.parse(str(genome), "fasta"))
    codon_table = CodonTable.unambiguous_dna_by_id[genetic_code]

    for mrna in gff3db.features_of_type("mRNA"):
        cds_seqs = [
            contigs[child.chrom][child.start - 1 : child.end].seq
            for child in gff3db.children(mrna, featuretype="CDS", order_by="start")
        ]
        if not cds_seqs:
            continue

        seq_concat = Seq("".join(str(s) for s in cds_seqs))
        if mrna.strand == "-":
            seq_concat = seq_concat.reverse_complement()

        if extract_to == "protein":
            translated = seq_concat.translate(table=codon_table)
            yield SeqRecord(translated, id=mrna.id, description="")
        else:
            yield SeqRecord(seq_concat, id=mrna.id, description="")
