"""Tests for eukan.annotation.consensus — transcript-ORF reintroduction.

``build_consensus_models`` itself is exercised end-to-end (combinr through the
prettify tail) in test_combinr_consensus.py. This file covers the consensus-only
``_patch_in_unmatched_orfs`` step, which folds transcript ORFs that don't overlap
any consensus gene back into the model set.
"""

from __future__ import annotations

from eukan.annotation import consensus
from eukan.annotation.consensus import _patch_in_unmatched_orfs, _write_prettified_gff3
from eukan.gff import create_gff_db


def _gene_gff(
    chrom: str, start: int, end: int, gid: str, *, cds: bool = True, strand: str = "+",
) -> str:
    rows = [
        f"{chrom}\tx\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gid}",
        f"{chrom}\tx\tmRNA\t{start}\t{end}\t.\t{strand}\t.\tID={gid}.t1;Parent={gid}",
    ]
    if cds:
        rows.append(
            f"{chrom}\tx\tCDS\t{start}\t{end}\t.\t{strand}\t0\tID={gid}.cds;Parent={gid}.t1"
        )
    return "##gff-version 3\n" + "\n".join(rows) + "\n"


def _gene_ids(db) -> set[str]:
    return {f.id for f in db.features_of_type("gene")}


# One consensus gene as combinr writes it.
_CONSENSUS_GFF = (
    "##gff-version 3\n"
    "chr1\tcombinr\tgene\t1\t300\t.\t+\t.\tID=g1\n"
    "chr1\tcombinr\tmRNA\t1\t300\t.\t+\t.\tID=g1.t1;Parent=g1;support=consensus\n"
    "chr1\tcombinr\texon\t1\t300\t.\t+\t.\tID=g1.exon1;Parent=g1.t1\n"
    "chr1\tcombinr\tCDS\t1\t300\t.\t+\t0\tID=g1.cds1;Parent=g1.t1\n"
)


def _orf_gff(start: int, end: int, gid: str) -> str:
    """One transcript ORF gene on chr1 as find_orfs writes it: gene and mRNA rows
    derived by gffutils, exon and CDS rows from the assembly."""
    rows = [
        f"chr1\tgffutils_derived\tgene\t{start}\t{end}\t.\t+\t.\tID={gid}",
        f"chr1\tgffutils_derived\tmRNA\t{start}\t{end}\t.\t+\t.\tParent={gid};ID={gid}_mRNA",
        f"chr1\tcombinr-assembly\texon\t{start}\t{end}\t.\t+\t.\tID={gid}:exon:1;Parent={gid}_mRNA",
        f"chr1\tcombinr-assembly\tCDS\t{start}\t{end}\t.\t+\t0\tID={gid}:CDS:1;Parent={gid}_mRNA",
    ]
    return "##gff-version 3\n" + "\n".join(rows) + "\n"


class TestPatchInUnmatchedOrfs:
    def _consdb(self, tmp_path):
        path = tmp_path / "consensus_models.gff3"
        path.write_text(_gene_gff("chr1", 1, 300, "g1"))
        return create_gff_db(path)

    def test_no_orf_file_returns_unchanged(self, tmp_path):
        consdb = self._consdb(tmp_path)
        result = _patch_in_unmatched_orfs(consdb, tmp_path / "missing.gff3")
        assert result is consdb

    def test_nonoverlapping_orf_is_added(self, tmp_path):
        consdb = self._consdb(tmp_path)
        orf = tmp_path / "transcript_orfs.gff3"
        orf.write_text(_gene_gff("chr1", 1000, 1300, "orf1"))

        result = _patch_in_unmatched_orfs(consdb, orf)

        assert _gene_ids(result) == {"g1", "orf1"}

    def test_orfs_are_added_in_place(self, tmp_path):
        """The ORFs go into the consensus DB itself, not a rebuilt copy."""
        consdb = self._consdb(tmp_path)
        orf = tmp_path / "transcript_orfs.gff3"
        orf.write_text(_gene_gff("chr1", 1000, 1300, "orf1"))

        result = _patch_in_unmatched_orfs(consdb, orf)

        assert result is consdb
        assert _gene_ids(consdb) == {"g1", "orf1"}
        assert [f.id for f in consdb.children("orf1.t1")] == ["orf1.cds"]

    def test_overlapping_orf_not_added(self, tmp_path):
        consdb = self._consdb(tmp_path)
        orf = tmp_path / "transcript_orfs.gff3"
        orf.write_text(_gene_gff("chr1", 50, 250, "orf1"))  # overlaps g1

        result = _patch_in_unmatched_orfs(consdb, orf)

        assert result is consdb

    def test_antisense_orf_not_added(self, tmp_path):
        """An ORF on the opposite strand of a consensus gene overlaps it too."""
        consdb = self._consdb(tmp_path)
        orf = tmp_path / "transcript_orfs.gff3"
        orf.write_text(_gene_gff("chr1", 50, 250, "orf1", strand="-"))  # antisense to g1

        result = _patch_in_unmatched_orfs(consdb, orf)

        assert result is consdb
        assert _gene_ids(consdb) == {"g1"}

    def test_log_counts_genes_not_features(self, tmp_path, caplog):
        """One ORF gene (gene + mRNA + CDS) is one reintroduced ORF, not three."""
        consdb = self._consdb(tmp_path)
        orf = tmp_path / "transcript_orfs.gff3"
        orf.write_text(_gene_gff("chr1", 1000, 1300, "orf1"))

        with caplog.at_level("INFO", logger="eukan.annotation.consensus"):
            _patch_in_unmatched_orfs(consdb, orf)

        assert "Reintroduced 1 transcript ORFs" in caplog.text

    def test_orf_without_cds_not_added(self, tmp_path):
        consdb = self._consdb(tmp_path)
        orf = tmp_path / "transcript_orfs.gff3"
        orf.write_text(_gene_gff("chr1", 1000, 1300, "orf1", cds=False))

        result = _patch_in_unmatched_orfs(consdb, orf)

        assert result is consdb


class TestPatchLabels:
    """The patched genes stay recognisable: every row the patch adds says
    eukan-orf-patch in column 2 and every mRNA it adds carries
    support=transcript_orf, beside combinr's support=consensus."""

    def _patched(self, tmp_path):
        cons = tmp_path / "consensus_models.gff3"
        cons.write_text(_CONSENSUS_GFF)
        orf = tmp_path / "transcript_orfs.gff3"
        orf.write_text(_orf_gff(1000, 1299, "orf1"))
        return _patch_in_unmatched_orfs(create_gff_db(cons), orf)

    def test_added_rows_are_labelled_and_consensus_rows_are_not(self, tmp_path):
        db = self._patched(tmp_path)

        assert {f.id: f.source for f in db.all_features()} == {
            "g1": "combinr", "g1.t1": "combinr", "g1.exon1": "combinr", "g1.cds1": "combinr",
            "orf1": "eukan-orf-patch", "orf1_mRNA": "eukan-orf-patch",
            "orf1:exon:1": "eukan-orf-patch", "orf1:CDS:1": "eukan-orf-patch",
        }

    def test_support_attribute_on_added_mrnas_only(self, tmp_path):
        db = self._patched(tmp_path)

        assert {f.id: f.attributes.get("support") for f in db.all_features()} == {
            "g1": None, "g1.t1": ["consensus"], "g1.exon1": None, "g1.cds1": None,
            "orf1": None, "orf1_mRNA": ["transcript_orf"], "orf1:exon:1": None,
            "orf1:CDS:1": None,
        }

    def test_labels_reach_final_gff3(self, tmp_path):
        """Every final.gff3 row is eukan or, for a patched gene, eukan-orf-patch."""
        _write_prettified_gff3(self._patched(tmp_path), "TEST", tmp_path / "final.gff3")

        rows = [ln.split("\t") for ln in (tmp_path / "final.gff3").read_text().splitlines()]
        assert [(c[2], c[1]) for c in rows] == [
            ("gene", "eukan"), ("mRNA", "eukan"), ("exon", "eukan"), ("CDS", "eukan"),
            ("gene", "eukan-orf-patch"), ("mRNA", "eukan-orf-patch"),
            ("exon", "eukan-orf-patch"), ("CDS", "eukan-orf-patch"),
        ]
        consensus_mrna, patched_mrna = (c[8] for c in rows if c[2] == "mRNA")
        assert "support=consensus" in consensus_mrna
        assert "support=transcript_orf" in patched_mrna


class TestWritePrettifiedGff3:
    def test_phases_are_fixed_without_rebuilding_the_db(self, tmp_path, monkeypatch):
        path = tmp_path / "consensus_models.gff3"
        path.write_text(
            "chr1\tx\tgene\t1\t600\t.\t+\t.\tID=g1\n"
            "chr1\tx\tmRNA\t1\t600\t.\t+\t.\tID=g1.t1;Parent=g1\n"
            "chr1\tx\tCDS\t1\t100\t.\t+\t0\tID=g1.cds1;Parent=g1.t1\n"
            "chr1\tx\tCDS\t301\t600\t.\t+\t0\tID=g1.cds2;Parent=g1.t1\n"
        )
        consdb = create_gff_db(path)
        builds = []
        real_create = consensus.create_gff_db
        monkeypatch.setattr(
            consensus, "create_gff_db",
            lambda *a, **k: builds.append(a) or real_create(*a, **k),
        )

        _write_prettified_gff3(consdb, "TEST", tmp_path / "final.gff3")

        assert builds == []
        cds = [ln.split("\t") for ln in (tmp_path / "final.gff3").read_text().splitlines()
               if ln.split("\t")[2] == "CDS"]
        # 100 nt at phase 0 leaves phase 2 for the second CDS.
        assert [(c[3], c[7]) for c in cds] == [("1", "0"), ("301", "2")]
