"""Tests for eukan.annotation.alignment line scans (protein hints, gene stats)."""

import statistics

import pytest

from eukan.annotation.alignment import _gene_lengths_and_introns, _write_protein_hints
from eukan.gff import create_gff_db
from eukan.gff.transforms import prot2augustus_hints

# gth- and spaln-style prot.gff3 rows as eukan writes them, plus a
# multi-parent CDS and URL-escaped Parent values.
_PROT = (
    "##gff-version 3\n"
    "c1\tprot_align\tgene\t448\t2950\t0.465\t-\t.\tID=gene1;Target=KQC44209.1 1 964 +;\n"
    "c1\tprot_align\tmRNA\t448\t2950\t.\t-\t.\tID=gene1_mRNA;Parent=gene1;\n"
    "c1\tprot_align\texon\t448\t2950\t0.465\t-\t.\tParent=gene1_mRNA;ID=exon_1;\n"
    "c1\tprot_align\tCDS\t448\t2950\t.\t-\t0\tParent=gene1_mRNA;ID=CDS_1;\n"
    "# a comment\n"
    "c2\tprot_align\tgene\t21000\t23371\t.\t+\t.\tID=gene00020;Name=prot20;\n"
    "c2\tprot_align\tmRNA\t21000\t23371\t.\t+\t.\tID=mRNA00020;Parent=gene00020;Name=prot20;\n"
    "c2\tprot_align\tCDS\t21000\t21104\t.\t+\t0\tID=mRNA00020:cds00020_1;Parent=mRNA00020;"
    "Name=prot20;Target=prot20 1 67 +;\n"
    "c2\tprot_align\tCDS\t21909\t22100\t12.5\t+\t2\tID=mRNA00020:cds00020_2;Parent=mRNA00020;"
    "Name=prot20;Target=prot20 1 67 +;\n"
    "c3\tprot_align\tCDS\t10\t90\t.\t+\t1\tID=shared;Parent=mA,mB;\n"
    "c3\tprot_align\tCDS\t100\t190\t.\t-\t0\tID=esc;Parent=m%3Bx%3Dy%26z%25w%2Cv;\n"
)


@pytest.mark.parametrize("trailing", [True, False])
def test_protein_hints_match_the_gffutils_transform(tmp_path, trailing):
    """The line scan writes the bytes the prot2augustus_hints DB pass wrote."""
    text = _PROT if trailing else _PROT.replace(";\n", "\n")
    prot = tmp_path / "prot.gff3"
    prot.write_text(text)
    ref = "".join(
        f"{f}\n" for f in create_gff_db(prot, transform=prot2augustus_hints).all_features()
    )
    out = tmp_path / "hints_protein.gff"

    _write_protein_hints(prot, out)

    assert out.read_text() == ref
    assert len(ref.splitlines()) == 5


@pytest.mark.parametrize("trailing", [True, False])
def test_gene_lengths_and_introns_match_the_db_walk(tmp_path, trailing):
    # Without trailing semicolons the old scan kept the newline on the last
    # attribute, so g1.t1's CDS (Parent last) missed g1 and g1 lost an intron.
    gff = tmp_path / "genemark.gff3"
    text = (
        "c1\tgenemark\tgene\t100\t900\t.\t+\t.\tID=g1;\n"
        "c1\tgenemark\tmRNA\t100\t900\t.\t+\t.\tID=g1.t1;Parent=g1;\n"
        "c1\tgenemark\tCDS\t100\t200\t.\t+\t0\tID=g1.c1;Parent=g1.t1;\n"
        "c1\tgenemark\tCDS\t300\t400\t.\t+\t0\tID=g1.c2;Parent=g1.t1;\n"
        "c1\tgenemark\tCDS\t800\t900\t.\t+\t0\tID=g1.c3;Parent=g1.t1;\n"
        "c1\tgenemark\tmRNA\t100\t400\t.\t+\t.\tID=g1.t2;Parent=g1;\n"
        "c1\tgenemark\tCDS\t100\t200\t.\t+\t0\tID=g1.c4;Parent=g1.t2;\n"
        "c1\tgenemark\tgene\t2000\t2300\t.\t-\t.\tID=g2;\n"
        "c1\tgenemark\tmRNA\t2000\t2300\t.\t-\t.\tID=g2.t1;Parent=g2;\n"
        "c1\tgenemark\tCDS\t2000\t2300\t.\t-\t0\tID=g2.c1;Parent=g2.t1;\n"
        "c1\tgenemark\tgene\t5000\t5600\t.\t+\t.\tID=g3;\n"
        "c1\tgenemark\tmRNA\t5000\t5600\t.\t+\t.\tID=g3.t1;Parent=g3;\n"
    )
    gff.write_text(text if trailing else text.replace(";\n", "\n"))
    db = create_gff_db(gff)
    genes = list(db.features_of_type("gene"))
    db_lengths = [abs(g.end - g.start) for g in genes]
    db_introns = sum(max(0, len(list(db.children(g, featuretype="CDS"))) - 1) for g in genes)

    lengths, introns = _gene_lengths_and_introns(gff)

    assert (lengths, introns) == (db_lengths, db_introns) == ([800, 300, 600], 3)
    assert statistics.stdev(lengths) == statistics.stdev(db_lengths)
