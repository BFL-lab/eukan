"""Tests for eukan.stats.compare — annotation comparison."""

from __future__ import annotations

from pathlib import Path

from eukan.stats import compare_annotations
from eukan.stats.compare import _stream_features

# A minimal GFF3 with two genes: one single-CDS, one with two CDS (one intron).
_REF_GFF = """##gff-version 3
chr1\ttest\tgene\t100\t500\t.\t+\t.\tID=g1
chr1\ttest\tmRNA\t100\t500\t.\t+\t.\tID=m1;Parent=g1
chr1\ttest\tCDS\t100\t500\t.\t+\t0\tID=c1;Parent=m1
chr1\ttest\tgene\t1000\t2500\t.\t-\t.\tID=g2
chr1\ttest\tmRNA\t1000\t2500\t.\t-\t.\tID=m2;Parent=g2
chr1\ttest\tCDS\t1000\t1400\t.\t-\t0\tID=c2a;Parent=m2
chr1\ttest\tCDS\t2000\t2500\t.\t-\t1\tID=c2b;Parent=m2
"""


def _write(path: Path, content: str) -> Path:
    path.write_text(content)
    return path


class TestStreamFeatures:
    def test_extracts_gene_mrna_cds(self, tmp_path):
        gff = _write(tmp_path / "ref.gff3", _REF_GFF)
        genes, mrnas, cdss = _stream_features(gff)
        assert [g.feat_id for g in genes] == ["g1", "g2"]
        assert [m.parent_id for m in mrnas] == ["g1", "g2"]
        assert [c.feat_id for c in cdss] == ["c1", "c2a", "c2b"]
        assert genes[1].strand == "-"
        assert cdss[2].start == 2000

    def test_skips_other_feature_types(self, tmp_path):
        extra = _REF_GFF + "chr1\ttest\texon\t100\t500\t.\t+\t.\tID=e1;Parent=m1\n"
        gff = _write(tmp_path / "ref.gff3", extra)
        genes, mrnas, cdss = _stream_features(gff)
        assert len(genes) == 2
        assert len(mrnas) == 2
        assert len(cdss) == 3  # exon ignored

    def test_stops_at_fasta_directive(self, tmp_path):
        fasta_section = "##FASTA\n>chr1\nACGT\n" * 1000
        gff = _write(tmp_path / "ref.gff3", _REF_GFF + fasta_section)
        genes, _, _ = _stream_features(gff)
        # Should not have parsed the FASTA section as feature lines
        assert len(genes) == 2

    def test_first_parent_for_multi_parent(self, tmp_path):
        gff_text = (
            "##gff-version 3\n"
            "chr1\ttest\tgene\t100\t500\t.\t+\t.\tID=g1\n"
            "chr1\ttest\tgene\t100\t500\t.\t+\t.\tID=g2\n"
            "chr1\ttest\tmRNA\t100\t500\t.\t+\t.\tID=m1;Parent=g1,g2\n"
        )
        gff = _write(tmp_path / "ref.gff3", gff_text)
        _, mrnas, _ = _stream_features(gff)
        assert mrnas[0].parent_id == "g1"


class TestCompareAnnotations:
    def test_identical_inputs_all_match(self, tmp_path):
        gff = _write(tmp_path / "ref.gff3", _REF_GFF)
        result = compare_annotations(gff, gff)
        assert result.gene_stats.exact == 2
        assert result.gene_stats.inexact == 0
        assert result.gene_stats.missing == 0
        assert result.mrna_stats.match == 2
        assert result.cds_stats.match == 3
        assert result.intron_stats.match == 1  # one intron in g2

    def test_missing_gene_in_pred(self, tmp_path):
        ref = _write(tmp_path / "ref.gff3", _REF_GFF)
        # Predicted has only g1, not g2
        pred_text = "\n".join(_REF_GFF.splitlines()[:4]) + "\n"
        pred = _write(tmp_path / "pred.gff3", pred_text)
        result = compare_annotations(ref, pred)
        assert result.gene_stats.exact == 1
        assert result.gene_stats.missing == 1

    def test_inexact_boundary(self, tmp_path):
        ref = _write(tmp_path / "ref.gff3", _REF_GFF)
        # Shift g1's gene boundary by +10
        pred_text = _REF_GFF.replace(
            "gene\t100\t500\t.\t+\t.\tID=g1",
            "gene\t110\t510\t.\t+\t.\tID=g1",
        )
        pred = _write(tmp_path / "pred.gff3", pred_text)
        result = compare_annotations(ref, pred)
        assert result.gene_stats.exact == 1  # g2 still exact
        assert result.gene_stats.inexact == 1  # g1 shifted
