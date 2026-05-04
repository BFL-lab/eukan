"""Tests for eukan.gff.parser — GFF3 transforms and feature hierarchy."""

from eukan.gff.parser import (
    add_missing_feats_to_gff3,
    fix_CDS_phases,
    prettify_gff3,
)

# ---------------------------------------------------------------------------
# add_missing_feats_to_gff3
# ---------------------------------------------------------------------------


class TestAddMissingFeats:
    def test_adds_mRNA_when_missing(self, db_from_string):
        """Gene + CDS but no mRNA → should generate mRNA features."""
        gff = (
            "chr1\ttest\tgene\t100\t500\t.\t+\t.\tID=gene1\n"
            "chr1\ttest\tCDS\t100\t300\t.\t+\t0\tID=cds1;Parent=gene1\n"
            "chr1\ttest\tCDS\t400\t500\t.\t+\t0\tID=cds2;Parent=gene1\n"
        )
        db = db_from_string(gff)
        new_feats = list(add_missing_feats_to_gff3(db))

        mrna_feats = [f for f in new_feats if f.featuretype == "mRNA"]
        assert len(mrna_feats) == 1
        assert mrna_feats[0].start == 100
        assert mrna_feats[0].end == 500

    def test_adds_exon_from_CDS(self, db_from_string):
        """Gene + mRNA + CDS but no exon → should derive exons from CDS."""
        gff = (
            "chr1\ttest\tgene\t100\t500\t.\t+\t.\tID=gene1\n"
            "chr1\ttest\tmRNA\t100\t500\t.\t+\t.\tID=mrna1;Parent=gene1\n"
            "chr1\ttest\tCDS\t100\t300\t.\t+\t0\tID=cds1;Parent=mrna1\n"
            "chr1\ttest\tCDS\t400\t500\t.\t+\t1\tID=cds2;Parent=mrna1\n"
        )
        db = db_from_string(gff)
        new_feats = list(add_missing_feats_to_gff3(db))

        exon_feats = [f for f in new_feats if f.featuretype == "exon"]
        assert len(exon_feats) == 2
        assert exon_feats[0].start == 100
        assert exon_feats[1].start == 400

    def test_adds_CDS_from_exon(self, db_from_string):
        """Gene + mRNA + exon but no CDS → should derive CDS from exons."""
        gff = (
            "chr1\ttest\tgene\t100\t500\t.\t+\t.\tID=gene1\n"
            "chr1\ttest\tmRNA\t100\t500\t.\t+\t.\tID=mrna1;Parent=gene1\n"
            "chr1\ttest\texon\t100\t300\t.\t+\t.\tID=exon1;Parent=mrna1\n"
            "chr1\ttest\texon\t400\t500\t.\t+\t.\tID=exon2;Parent=mrna1\n"
        )
        db = db_from_string(gff)
        new_feats = list(add_missing_feats_to_gff3(db))

        cds_feats = [f for f in new_feats if f.featuretype == "CDS"]
        assert len(cds_feats) == 2

    def test_no_changes_when_complete(self, db_from_string):
        """Complete hierarchy should not generate any new features."""
        gff = (
            "chr1\ttest\tgene\t100\t500\t.\t+\t.\tID=gene1\n"
            "chr1\ttest\tmRNA\t100\t500\t.\t+\t.\tID=mrna1;Parent=gene1\n"
            "chr1\ttest\texon\t100\t300\t.\t+\t.\tID=exon1;Parent=mrna1\n"
            "chr1\ttest\tCDS\t100\t300\t.\t+\t0\tID=cds1;Parent=mrna1\n"
        )
        db = db_from_string(gff)
        new_feats = list(add_missing_feats_to_gff3(db))
        assert len(new_feats) == 0


# ---------------------------------------------------------------------------
# fix_CDS_phases
# ---------------------------------------------------------------------------


class TestFixCDSPhases:
    def test_single_cds_phase_zero(self, db_from_string):
        """Single-CDS gene should have phase 0."""
        gff = (
            "chr1\ttest\tgene\t100\t400\t.\t+\t.\tID=gene1\n"
            "chr1\ttest\tmRNA\t100\t400\t.\t+\t.\tID=mrna1;Parent=gene1\n"
            "chr1\ttest\texon\t100\t400\t.\t+\t.\tID=exon1;Parent=mrna1\n"
            "chr1\ttest\tCDS\t100\t400\t.\t+\t.\tID=cds1;Parent=mrna1\n"
        )
        db = db_from_string(gff)
        features = list(fix_CDS_phases(db))

        cds_feats = [f for f in features if f.featuretype == "CDS"]
        assert len(cds_feats) == 1
        assert cds_feats[0].frame == "0"

    def test_multi_cds_phase_propagation(self, db_from_string):
        """Multi-CDS gene should have correct phase propagation."""
        gff = (
            "chr1\ttest\tgene\t100\t600\t.\t+\t.\tID=gene1\n"
            "chr1\ttest\tmRNA\t100\t600\t.\t+\t.\tID=mrna1;Parent=gene1\n"
            "chr1\ttest\texon\t100\t200\t.\t+\t.\tID=exon1;Parent=mrna1\n"
            "chr1\ttest\texon\t400\t600\t.\t+\t.\tID=exon2;Parent=mrna1\n"
            "chr1\ttest\tCDS\t100\t200\t.\t+\t.\tID=cds1;Parent=mrna1\n"
            "chr1\ttest\tCDS\t400\t600\t.\t+\t.\tID=cds2;Parent=mrna1\n"
        )
        db = db_from_string(gff)
        features = list(fix_CDS_phases(db))

        cds_feats = [f for f in features if f.featuretype == "CDS"]
        assert len(cds_feats) == 2
        assert cds_feats[0].frame == "0"

        # First CDS is 101 nt (100-200 inclusive), phase 0
        # next_phase = (3 - ((101 - 0) % 3)) % 3 = (3 - 2) % 3 = 1
        assert cds_feats[1].frame == "1"

    def test_minus_strand_phase(self, db_from_string):
        """Minus strand: phases should be computed from the last CDS first."""
        gff = (
            "chr1\ttest\tgene\t100\t600\t.\t-\t.\tID=gene1\n"
            "chr1\ttest\tmRNA\t100\t600\t.\t-\t.\tID=mrna1;Parent=gene1\n"
            "chr1\ttest\texon\t100\t200\t.\t-\t.\tID=exon1;Parent=mrna1\n"
            "chr1\ttest\texon\t400\t600\t.\t-\t.\tID=exon2;Parent=mrna1\n"
            "chr1\ttest\tCDS\t100\t200\t.\t-\t.\tID=cds1;Parent=mrna1\n"
            "chr1\ttest\tCDS\t400\t600\t.\t-\t.\tID=cds2;Parent=mrna1\n"
        )
        db = db_from_string(gff)
        features = list(fix_CDS_phases(db))

        cds_feats = [f for f in features if f.featuretype == "CDS"]
        assert len(cds_feats) == 2
        # On minus strand, the last CDS (400-600, 201 nt) is processed first
        assert cds_feats[0].frame == "0"


# ---------------------------------------------------------------------------
# prettify_gff3
# ---------------------------------------------------------------------------


class TestPrettifyGFF3:
    def test_sequential_locus_tags(self, db_from_string):
        """Genes should get sequential locus tags."""
        gff = (
            "chr1\ttest\tgene\t100\t400\t.\t+\t.\tID=gene1\n"
            "chr1\ttest\tmRNA\t100\t400\t.\t+\t.\tID=mrna1;Parent=gene1\n"
            "chr1\ttest\texon\t100\t400\t.\t+\t.\tID=exon1;Parent=mrna1\n"
            "chr1\ttest\tCDS\t100\t400\t.\t+\t0\tID=cds1;Parent=mrna1\n"
            "chr1\ttest\tgene\t600\t900\t.\t+\t.\tID=gene2\n"
            "chr1\ttest\tmRNA\t600\t900\t.\t+\t.\tID=mrna2;Parent=gene2\n"
            "chr1\ttest\texon\t600\t900\t.\t+\t.\tID=exon2;Parent=mrna2\n"
            "chr1\ttest\tCDS\t600\t900\t.\t+\t0\tID=cds2;Parent=mrna2\n"
        )
        db = db_from_string(gff)
        features = list(prettify_gff3(db, "TEST"))

        genes = [f for f in features if f.featuretype == "gene"]
        assert len(genes) == 2
        # gffutils stores attribute values as lists.
        assert genes[0].attributes["ID"] == ["TEST_00001"]
        assert genes[1].attributes["ID"] == ["TEST_00002"]
        assert all(f.source == "eukannotpass" for f in features)
