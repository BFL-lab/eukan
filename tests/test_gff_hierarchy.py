"""Tests for eukan.gff.hierarchy — feature hierarchy operations."""

import os
import subprocess
import sys

from eukan.gff.hierarchy import (
    add_missing_feats_to_gff3,
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

    def test_CDS_from_exons_is_phased_from_the_5_prime_exon(self, db_from_string):
        """Exon-only models (SNAP's output) are phased from each mRNA's
        5'-most exon, where the start codon is, not from genomic coordinates."""
        gff = (
            "chr1\ttest\tgene\t100\t600\t.\t+\t.\tID=gene1\n"
            "chr1\ttest\tmRNA\t100\t600\t.\t+\t.\tID=mrna1;Parent=gene1\n"
            "chr1\ttest\texon\t100\t201\t.\t+\t.\tID=exon1;Parent=mrna1\n"
            "chr1\ttest\texon\t400\t600\t.\t+\t.\tID=exon2;Parent=mrna1\n"
            "chr1\ttest\tgene\t1000\t1601\t.\t-\t.\tID=gene2\n"
            "chr1\ttest\tmRNA\t1000\t1601\t.\t-\t.\tID=mrna2;Parent=gene2\n"
            "chr1\ttest\texon\t1000\t1200\t.\t-\t.\tID=exon3;Parent=mrna2\n"
            "chr1\ttest\texon\t1400\t1601\t.\t-\t.\tID=exon4;Parent=mrna2\n"
        )
        phases = {
            (f.start, f.end): str(f.frame)
            for f in add_missing_feats_to_gff3(db_from_string(gff))
            if f.featuretype == "CDS"
        }
        # Plus: 102 nt at phase 0 leaves 0 for the next exon. Minus: 1400-1601
        # (202 nt) is the 5'-most exon, leaving (3 - 202 % 3) % 3 = 2.
        assert phases == {
            (100, 201): "0", (400, 600): "0", (1400, 1601): "0", (1000, 1200): "2",
        }

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

    def _cds_frames(self, db):
        return {
            (f.start, f.end): f.frame
            for f in prettify_gff3(db, "TEST") if f.featuretype == "CDS"
        }

    def test_single_cds_phase_zero(self, db_from_string):
        gff = (
            "chr1\ttest\tgene\t100\t400\t.\t+\t.\tID=gene1\n"
            "chr1\ttest\tmRNA\t100\t400\t.\t+\t.\tID=mrna1;Parent=gene1\n"
            "chr1\ttest\texon\t100\t400\t.\t+\t.\tID=exon1;Parent=mrna1\n"
            "chr1\ttest\tCDS\t100\t400\t.\t+\t2\tID=cds1;Parent=mrna1\n"
        )
        assert self._cds_frames(db_from_string(gff)) == {(100, 400): "0"}

    def test_multi_cds_phase_propagation(self, db_from_string):
        gff = (
            "chr1\ttest\tgene\t100\t600\t.\t+\t.\tID=gene1\n"
            "chr1\ttest\tmRNA\t100\t600\t.\t+\t.\tID=mrna1;Parent=gene1\n"
            "chr1\ttest\texon\t100\t200\t.\t+\t.\tID=exon1;Parent=mrna1\n"
            "chr1\ttest\texon\t400\t600\t.\t+\t.\tID=exon2;Parent=mrna1\n"
            "chr1\ttest\tCDS\t100\t200\t.\t+\t.\tID=cds1;Parent=mrna1\n"
            "chr1\ttest\tCDS\t400\t600\t.\t+\t.\tID=cds2;Parent=mrna1\n"
        )
        # 101 nt at phase 0 leaves (3 - 101 % 3) % 3 = 1 for the next CDS.
        assert self._cds_frames(db_from_string(gff)) == {(100, 200): "0", (400, 600): "1"}

    def test_minus_strand_phase_starts_at_the_last_cds(self, db_from_string):
        gff = (
            "chr1\ttest\tgene\t100\t601\t.\t-\t.\tID=gene1\n"
            "chr1\ttest\tmRNA\t100\t601\t.\t-\t.\tID=mrna1;Parent=gene1\n"
            "chr1\ttest\texon\t100\t200\t.\t-\t.\tID=exon1;Parent=mrna1\n"
            "chr1\ttest\texon\t400\t601\t.\t-\t.\tID=exon2;Parent=mrna1\n"
            "chr1\ttest\tCDS\t100\t200\t.\t-\t.\tID=cds1;Parent=mrna1\n"
            "chr1\ttest\tCDS\t400\t601\t.\t-\t.\tID=cds2;Parent=mrna1\n"
        )
        # The 5'-most CDS on the minus strand (400-601, 202 nt) is phase 0.
        assert self._cds_frames(db_from_string(gff)) == {(400, 601): "0", (100, 200): "2"}

    def test_only_exon_and_cds_children_are_written(self, db_from_string):
        """UTR features are left out; exons extending past the CDS carry them."""
        gff = (
            "chr1\ttest\tgene\t100\t600\t.\t+\t.\tID=gene1\n"
            "chr1\ttest\tmRNA\t100\t600\t.\t+\t.\tID=mrna1;Parent=gene1\n"
            "chr1\ttest\texon\t100\t600\t.\t+\t.\tID=exon1;Parent=mrna1\n"
            "chr1\ttest\tfive_prime_UTR\t100\t199\t.\t+\t.\tID=u5;Parent=mrna1\n"
            "chr1\ttest\tCDS\t200\t500\t.\t+\t0\tID=cds1;Parent=mrna1\n"
            "chr1\ttest\tthree_prime_UTR\t501\t600\t.\t+\t.\tID=u3;Parent=mrna1\n"
        )
        types = sorted(f.featuretype for f in prettify_gff3(db_from_string(gff), "TEST"))
        assert types == ["CDS", "exon", "gene", "mRNA"]

    def test_block_order_does_not_depend_on_the_hash_seed(self):
        """Each mRNA's exons come before its CDSs under any PYTHONHASHSEED, so
        final.gff3 is byte-stable across runs."""
        gff = (
            "chr1\ttest\tgene\t100\t600\t.\t+\t.\tID=gene1\n"
            "chr1\ttest\tmRNA\t100\t600\t.\t+\t.\tID=mrna1;Parent=gene1\n"
            "chr1\ttest\tCDS\t100\t200\t.\t+\t0\tID=cds1;Parent=mrna1\n"
            "chr1\ttest\texon\t100\t200\t.\t+\t.\tID=exon1;Parent=mrna1\n"
            "chr1\ttest\tCDS\t400\t600\t.\t+\t1\tID=cds2;Parent=mrna1\n"
            "chr1\ttest\texon\t400\t600\t.\t+\t.\tID=exon2;Parent=mrna1\n"
        )
        code = (
            "import sys, gffutils\n"
            "from eukan.gff.hierarchy import prettify_gff3\n"
            "db = gffutils.create_db(sys.stdin.read(), ':memory:', from_string=True)\n"
            "for f in prettify_gff3(db, 'TEST'):\n"
            "    print(f)\n"
        )
        # Seeds 0 and 3 iterate {"exon", "CDS"} in opposite orders.
        outs = [
            subprocess.run(
                [sys.executable, "-c", code], input=gff, capture_output=True, text=True,
                check=True, env={**os.environ, "PYTHONHASHSEED": seed},
            ).stdout
            for seed in ("0", "3")
        ]
        assert outs[0] == outs[1]
        types = [line.split("\t")[2] for line in outs[0].splitlines()]
        assert types == ["gene", "mRNA", "exon", "exon", "CDS", "CDS"]
