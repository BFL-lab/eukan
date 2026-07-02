"""Tests for eukan.annotation.consensus — transcript-ORF reintroduction.

``build_consensus_models`` itself is exercised end-to-end (combinr through the
prettify tail) in test_combinr_consensus.py. This file covers the consensus-only
``_patch_in_unmatched_orfs`` step, which folds transcript ORFs that don't overlap
any consensus gene back into the model set.
"""

from __future__ import annotations

from eukan.annotation.consensus import _patch_in_unmatched_orfs
from eukan.gff import create_gff_db


def _gene_gff(chrom: str, start: int, end: int, gid: str, *, cds: bool = True) -> str:
    rows = [
        f"{chrom}\tx\tgene\t{start}\t{end}\t.\t+\t.\tID={gid}",
        f"{chrom}\tx\tmRNA\t{start}\t{end}\t.\t+\t.\tID={gid}.t1;Parent={gid}",
    ]
    if cds:
        rows.append(f"{chrom}\tx\tCDS\t{start}\t{end}\t.\t+\t0\tID={gid}.cds;Parent={gid}.t1")
    return "##gff-version 3\n" + "\n".join(rows) + "\n"


def _gene_ids(db) -> set[str]:
    return {f.id for f in db.features_of_type("gene")}


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

    def test_overlapping_orf_not_added(self, tmp_path):
        consdb = self._consdb(tmp_path)
        orf = tmp_path / "transcript_orfs.gff3"
        orf.write_text(_gene_gff("chr1", 50, 250, "orf1"))  # overlaps g1

        result = _patch_in_unmatched_orfs(consdb, orf)

        assert result is consdb

    def test_orf_without_cds_not_added(self, tmp_path):
        consdb = self._consdb(tmp_path)
        orf = tmp_path / "transcript_orfs.gff3"
        orf.write_text(_gene_gff("chr1", 1000, 1300, "orf1", cds=False))

        result = _patch_in_unmatched_orfs(consdb, orf)

        assert result is consdb
