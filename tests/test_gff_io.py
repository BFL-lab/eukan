"""Tests for eukan.gff.io — GFF3 serialization and sequence extraction."""

import warnings

from eukan.gff.io import extract_sequences, featuredb2gff3_file


class TestFeaturedb2Gff3File:
    def test_writes_hierarchical_gff3(self, tmp_path, db_from_string):
        """Output should contain gene > mRNA > exon + CDS hierarchy."""
        gff = (
            "chr1\ttest\tgene\t100\t400\t.\t+\t.\tID=gene1\n"
            "chr1\ttest\tmRNA\t100\t400\t.\t+\t.\tID=mrna1;Parent=gene1\n"
            "chr1\ttest\texon\t100\t400\t.\t+\t.\tID=exon1;Parent=mrna1\n"
            "chr1\ttest\tCDS\t100\t400\t.\t+\t0\tID=cds1;Parent=mrna1\n"
        )
        db = db_from_string(gff)
        out_path = tmp_path / "test.gff3"
        featuredb2gff3_file(db, out_path)

        lines = out_path.read_text().strip().split("\n")
        assert len(lines) == 4
        assert "gene" in lines[0]
        assert "mRNA" in lines[1]
        assert "exon" in lines[2]
        assert "CDS" in lines[3]


class TestExtractSequences:
    def test_protein_extraction(self, tmp_path):
        """Should extract translated protein sequences from GFF3 + genome."""
        # Create a simple genome with one contig
        genome_path = tmp_path / "genome.fa"
        # ATG AAA AAA TAA = M K K *
        # Start at pos 10, end at pos 21 (1-based, inclusive)
        seq = "N" * 9 + "ATGAAAAAA" + "TAA" + "N" * 78
        genome_path.write_text(f">chr1\n{seq}\n")

        gff_path = tmp_path / "test.gff3"
        gff_path.write_text(
            "chr1\ttest\tgene\t10\t21\t.\t+\t.\tID=gene1\n"
            "chr1\ttest\tmRNA\t10\t21\t.\t+\t.\tID=mrna1;Parent=gene1\n"
            "chr1\ttest\tCDS\t10\t21\t.\t+\t0\tID=cds1;Parent=mrna1\n"
        )

        records = list(extract_sequences(gff_path, genome_path, "protein", 1))
        assert len(records) == 1
        # ATG AAA AAA TAA -> MKK*
        assert str(records[0].seq).startswith("MKK")

    def test_protein_is_read_from_the_cds_phase(self, tmp_path):
        """A 5'-partial CDS in phase 1 skips one base before its first codon,
        and the trailing partial codon is left untranslated."""
        genome_path = tmp_path / "genome.fa"
        # C | ATG AAA AAA | TA, at 10-21
        seq = "N" * 9 + "CATGAAAAAATA" + "N" * 79
        genome_path.write_text(f">chr1\n{seq}\n")
        gff_path = tmp_path / "test.gff3"
        gff_path.write_text(
            "chr1\ttest\tgene\t10\t21\t.\t+\t.\tID=gene1\n"
            "chr1\ttest\tmRNA\t10\t21\t.\t+\t.\tID=mrna1;Parent=gene1\n"
            "chr1\ttest\tCDS\t10\t21\t.\t+\t1\tID=cds1;Parent=mrna1\n"
        )

        with warnings.catch_warnings():
            warnings.simplefilter("error")
            records = list(extract_sequences(gff_path, genome_path, "protein", 1))

        assert str(records[0].seq) == "MKK"

    def test_minus_strand_phase_is_read_at_the_5_prime_end(self, tmp_path):
        genome_path = tmp_path / "genome.fa"
        # Reverse complement of GC | ATG AAA TAA, at 10-20
        seq = "N" * 9 + "TTATTTCATGC" + "N" * 80
        genome_path.write_text(f">chr1\n{seq}\n")
        gff_path = tmp_path / "test.gff3"
        gff_path.write_text(
            "chr1\ttest\tgene\t10\t20\t.\t-\t.\tID=gene1\n"
            "chr1\ttest\tmRNA\t10\t20\t.\t-\t.\tID=mrna1;Parent=gene1\n"
            "chr1\ttest\tCDS\t10\t20\t.\t-\t2\tID=cds1;Parent=mrna1\n"
        )

        records = list(extract_sequences(gff_path, genome_path, "protein", 1))

        assert str(records[0].seq) == "MK*"

    def test_cdna_extraction(self, tmp_path):
        """Should extract cDNA (nucleotide) sequences."""
        genome_path = tmp_path / "genome.fa"
        seq = "N" * 9 + "ATGAAATAA" + "N" * 82
        genome_path.write_text(f">chr1\n{seq}\n")

        gff_path = tmp_path / "test.gff3"
        gff_path.write_text(
            "chr1\ttest\tgene\t10\t18\t.\t+\t.\tID=gene1\n"
            "chr1\ttest\tmRNA\t10\t18\t.\t+\t.\tID=mrna1;Parent=gene1\n"
            "chr1\ttest\tCDS\t10\t18\t.\t+\t0\tID=cds1;Parent=mrna1\n"
        )

        records = list(extract_sequences(gff_path, genome_path, "cdna", 1))
        assert len(records) == 1
        assert str(records[0].seq) == "ATGAAATAA"
