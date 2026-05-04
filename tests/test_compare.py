"""Tests for eukan.stats.compare — annotation comparison."""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from eukan.stats import compare_annotations, compare_multiple
from eukan.stats.compare import _stream_features
from eukan.stats.models import MultiComparisonResult

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


# Reference with 4 genes — g2 has 2 CDS (one intron); the rest are single-CDS.
# Used by the multi-prediction test class below.
_MULTI_REF_GFF = """##gff-version 3
chr1\t.\tgene\t100\t500\t.\t+\t.\tID=g1
chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=m1;Parent=g1
chr1\t.\tCDS\t100\t500\t.\t+\t0\tID=c1;Parent=m1
chr1\t.\tgene\t1000\t2500\t.\t-\t.\tID=g2
chr1\t.\tmRNA\t1000\t2500\t.\t-\t.\tID=m2;Parent=g2
chr1\t.\tCDS\t1000\t1400\t.\t-\t0\tID=c2a;Parent=m2
chr1\t.\tCDS\t2000\t2500\t.\t-\t1\tID=c2b;Parent=m2
chr1\t.\tgene\t3000\t3500\t.\t+\t.\tID=g3
chr1\t.\tmRNA\t3000\t3500\t.\t+\t.\tID=m3;Parent=g3
chr1\t.\tCDS\t3000\t3500\t.\t+\t0\tID=c3;Parent=m3
chr1\t.\tgene\t5000\t6000\t.\t-\t.\tID=g4
chr1\t.\tmRNA\t5000\t6000\t.\t-\t.\tID=m4;Parent=g4
chr1\t.\tCDS\t5000\t6000\t.\t-\t0\tID=c4;Parent=m4
"""

# p1: identical to ref -> 4 exact genes
_MULTI_P1_GFF = _MULTI_REF_GFF

# p2: g4 missing; g1, g2, g3 exact
_MULTI_P2_GFF = """##gff-version 3
chr1\t.\tgene\t100\t500\t.\t+\t.\tID=g1
chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=m1;Parent=g1
chr1\t.\tCDS\t100\t500\t.\t+\t0\tID=c1;Parent=m1
chr1\t.\tgene\t1000\t2500\t.\t-\t.\tID=g2
chr1\t.\tmRNA\t1000\t2500\t.\t-\t.\tID=m2;Parent=g2
chr1\t.\tCDS\t1000\t1400\t.\t-\t0\tID=c2a;Parent=m2
chr1\t.\tCDS\t2000\t2500\t.\t-\t1\tID=c2b;Parent=m2
chr1\t.\tgene\t3000\t3500\t.\t+\t.\tID=g3
chr1\t.\tmRNA\t3000\t3500\t.\t+\t.\tID=m3;Parent=g3
chr1\t.\tCDS\t3000\t3500\t.\t+\t0\tID=c3;Parent=m3
"""

# p3: g3 shifted by +100 (inexact); g4 missing
_MULTI_P3_GFF = """##gff-version 3
chr1\t.\tgene\t100\t500\t.\t+\t.\tID=g1
chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=m1;Parent=g1
chr1\t.\tCDS\t100\t500\t.\t+\t0\tID=c1;Parent=m1
chr1\t.\tgene\t1000\t2500\t.\t-\t.\tID=g2
chr1\t.\tmRNA\t1000\t2500\t.\t-\t.\tID=m2;Parent=g2
chr1\t.\tCDS\t1000\t1400\t.\t-\t0\tID=c2a;Parent=m2
chr1\t.\tCDS\t2000\t2500\t.\t-\t1\tID=c2b;Parent=m2
chr1\t.\tgene\t3100\t3600\t.\t+\t.\tID=g3
chr1\t.\tmRNA\t3100\t3600\t.\t+\t.\tID=m3;Parent=g3
chr1\t.\tCDS\t3100\t3600\t.\t+\t0\tID=c3;Parent=m3
"""


class TestCompareMultiple:
    """Multi-prediction comparison driver."""

    @pytest.fixture
    def three_preds(self, tmp_path) -> tuple[Path, list[tuple[str, Path]]]:
        ref = _write(tmp_path / "ref.gff3", _MULTI_REF_GFF)
        return ref, [
            ("p1", _write(tmp_path / "p1.gff3", _MULTI_P1_GFF)),
            ("p2", _write(tmp_path / "p2.gff3", _MULTI_P2_GFF)),
            ("p3", _write(tmp_path / "p3.gff3", _MULTI_P3_GFF)),
        ]

    def test_returns_multi_result(self, three_preds):
        ref, preds = three_preds
        result = compare_multiple(ref, preds)
        assert isinstance(result, MultiComparisonResult)
        assert len(result.per_prediction) == 3
        assert [p.label for p in result.per_prediction] == ["p1", "p2", "p3"]

    def test_per_prediction_classifications(self, three_preds):
        ref, preds = three_preds
        result = compare_multiple(ref, preds)
        # p1 perfect: 4 exact
        assert result.per_prediction[0].gene_stats.exact == 4
        # p2: g1, g2, g3 exact, g4 missing
        assert result.per_prediction[1].gene_stats.exact == 3
        assert result.per_prediction[1].gene_stats.missing == 1
        # p3: g1, g2 exact, g3 inexact, g4 missing
        assert result.per_prediction[2].gene_stats.exact == 2
        assert result.per_prediction[2].gene_stats.inexact == 1
        assert result.per_prediction[2].gene_stats.missing == 1

    def test_pair_test_count_default_metric(self, three_preds):
        ref, preds = three_preds
        result = compare_multiple(ref, preds)
        # 4 levels * 3 pairs * (1 KS metric + 1 chi2) = 24
        assert len(result.pair_tests) == 4 * 3 * 2

    def test_pair_test_count_three_metrics(self, three_preds):
        ref, preds = three_preds
        result = compare_multiple(ref, preds, ecdf_metrics=("sn", "sp", "f1"))
        # 4 levels * 3 pairs * (3 KS metrics + 1 chi2) = 48
        assert len(result.pair_tests) == 4 * 3 * 4

    def test_pair_test_p_adj_set(self, three_preds):
        ref, preds = three_preds
        result = compare_multiple(ref, preds)
        # p_adj must be set (>= p_value, since BH only inflates)
        for t in result.pair_tests:
            assert t.p_adj >= t.p_value - 1e-12
            assert 0.0 <= t.p_adj <= 1.0

    def test_pair_test_levels_and_kinds(self, three_preds):
        ref, preds = three_preds
        result = compare_multiple(ref, preds, ecdf_metrics=("f1",))
        levels = {t.level for t in result.pair_tests}
        assert levels == {"gene", "mrna", "cds", "intron"}
        kinds = {t.test for t in result.pair_tests}
        assert kinds == {"ks_f1", "chi2_classification"}

    def test_kappa_matrix_off_diagonal_pairs(self, three_preds):
        ref, preds = three_preds
        result = compare_multiple(ref, preds)
        # 3 preds -> 3 off-diagonal pairs
        assert len(result.kappa_matrix) == 3
        assert set(result.kappa_matrix.keys()) == {
            ("p1", "p2"), ("p1", "p3"), ("p2", "p3"),
        }
        for k in result.kappa_matrix.values():
            assert -1.0 <= k <= 1.0 or math.isnan(k)

    def test_powerset_sums_to_ref_total(self, three_preds):
        ref, preds = three_preds
        result = compare_multiple(ref, preds)
        # Every ref gene falls in exactly one powerset bucket.
        total = sum(result.powerset_matched.values())
        assert total == result.per_prediction[0].gene_stats.ref_total

    def test_powerset_buckets(self, three_preds):
        ref, preds = three_preds
        result = compare_multiple(ref, preds)
        # g1, g2 matched by all three preds; g3 matched by all three (p3 inexact);
        # g4 matched only by p1.
        assert result.powerset_matched[("p1", "p2", "p3")] == 3
        assert result.powerset_matched[("p1",)] == 1

    def test_n_one_yields_no_inter_pred_stats(self, tmp_path):
        ref = _write(tmp_path / "ref.gff3", _MULTI_REF_GFF)
        p1 = _write(tmp_path / "p1.gff3", _MULTI_P1_GFF)
        result = compare_multiple(ref, [("only", p1)])
        assert len(result.per_prediction) == 1
        assert result.pair_tests == []
        assert result.kappa_matrix == {}
        assert result.powerset_matched == {}

    def test_empty_predictions_raises(self, tmp_path):
        ref = _write(tmp_path / "ref.gff3", _MULTI_REF_GFF)
        with pytest.raises(ValueError, match="non-empty"):
            compare_multiple(ref, [])

    def test_duplicate_label_raises(self, tmp_path):
        ref = _write(tmp_path / "ref.gff3", _MULTI_REF_GFF)
        p1 = _write(tmp_path / "p1.gff3", _MULTI_P1_GFF)
        with pytest.raises(ValueError, match="unique"):
            compare_multiple(ref, [("a", p1), ("a", p1)])

    def test_invalid_metric_raises(self, tmp_path):
        ref = _write(tmp_path / "ref.gff3", _MULTI_REF_GFF)
        p1 = _write(tmp_path / "p1.gff3", _MULTI_P1_GFF)
        p2 = _write(tmp_path / "p2.gff3", _MULTI_P2_GFF)
        with pytest.raises(ValueError, match="ECDF metric"):
            compare_multiple(ref, [("a", p1), ("b", p2)], ecdf_metrics=("foo",))
