"""Tests for eukan.gff.concordance — multi-source gene-model operations."""

from collections import Counter
from pathlib import Path

import gffutils
import pytest

from eukan.gff.concordance import (
    _cds_boundaries_match,
    _features_overlap,
    _find_overlapping_genes,
    combine_nonredundant_models,
    extract_supported_models,
    find_concordant_models,
    find_nonoverlapping_genes,
    merge_fully_overlapping_transcript_genes,
)


def _feature(chrom, start, end, strand="+", featuretype="gene", id="g1", **attrs):
    return gffutils.Feature(
        seqid=chrom, source="test", featuretype=featuretype,
        start=start, end=end, strand=strand, frame=".",
        attributes={"ID": id, **attrs},
    )


# ---------------------------------------------------------------------------
# Overlap detection
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("a", "b", "expected"),
    [
        (("chr1", 100, 500, "+"), ("chr1", 400, 800, "+"), True),
        (("chr1", 100, 300, "+"), ("chr1", 400, 800, "+"), False),
        (("chr1", 100, 500, "+"), ("chr1", 100, 500, "-"), False),
        (("chr1", 100, 500, "+"), ("chr2", 100, 500, "+"), False),
        (("chr1", 100, 200, "+"), ("chr1", 201, 300, "+"), False),
        (("chr1", 100, 200, "+"), ("chr1", 200, 300, "+"), True),
    ],
    ids=[
        "overlapping",
        "non_overlapping",
        "different_strand",
        "different_chrom",
        "adjacent",
        "touching",
    ],
)
def test_features_overlap(a, b, expected):
    assert _features_overlap(_feature(*a), _feature(*b)) == expected


class TestFindOverlappingGenes:
    def test_two_overlapping(self, db_from_string):
        db = db_from_string(
            "chr1\ttest\tgene\t100\t500\t.\t+\t.\tID=g1\n"
            "chr1\ttest\tgene\t300\t700\t.\t+\t.\tID=g2\n"
        )
        clusters = _find_overlapping_genes(db)
        assert len(clusters) == 1
        assert len(clusters[0]) == 2

    def test_no_overlap(self, db_from_string):
        db = db_from_string(
            "chr1\ttest\tgene\t100\t200\t.\t+\t.\tID=g1\n"
            "chr1\ttest\tgene\t500\t700\t.\t+\t.\tID=g2\n"
        )
        clusters = _find_overlapping_genes(db)
        assert len(clusters) == 0

    def test_different_strands_not_merged(self, db_from_string):
        db = db_from_string(
            "chr1\ttest\tgene\t100\t500\t.\t+\t.\tID=g1\n"
            "chr1\ttest\tgene\t100\t500\t.\t-\t.\tID=g2\n"
        )
        clusters = _find_overlapping_genes(db)
        assert len(clusters) == 0


# ---------------------------------------------------------------------------
# CDS boundary matching
# ---------------------------------------------------------------------------


class TestCdsBoundaryMatch:
    def test_single_cds_always_matches(self):
        a = [_feature("chr1", 100, 300, featuretype="CDS")]
        b = [_feature("chr1", 105, 295, featuretype="CDS")]
        assert _cds_boundaries_match(a, b)

    def test_matching_multi_cds(self):
        a = [
            _feature("chr1", 100, 200, featuretype="CDS"),
            _feature("chr1", 400, 600, featuretype="CDS"),
        ]
        b = [
            _feature("chr1", 95, 201, featuretype="CDS"),   # end within 3bp
            _feature("chr1", 399, 605, featuretype="CDS"),   # start within 3bp
        ]
        assert _cds_boundaries_match(a, b)

    def test_mismatched_boundaries(self):
        a = [
            _feature("chr1", 100, 200, featuretype="CDS"),
            _feature("chr1", 400, 600, featuretype="CDS"),
        ]
        b = [
            _feature("chr1", 100, 210, featuretype="CDS"),  # end off by 10
            _feature("chr1", 400, 600, featuretype="CDS"),
        ]
        assert not _cds_boundaries_match(a, b)

    def test_internal_cds_both_boundaries(self):
        """Internal CDS (not first or last) must match both boundaries."""
        a = [
            _feature("chr1", 100, 200, featuretype="CDS"),
            _feature("chr1", 300, 400, featuretype="CDS"),
            _feature("chr1", 500, 600, featuretype="CDS"),
        ]
        b = [
            _feature("chr1", 100, 201, featuretype="CDS"),
            _feature("chr1", 301, 401, featuretype="CDS"),  # both within 3bp
            _feature("chr1", 499, 600, featuretype="CDS"),
        ]
        assert _cds_boundaries_match(a, b)


# ---------------------------------------------------------------------------
# Merge overlapping genes
# ---------------------------------------------------------------------------


class TestMergeOverlappingGenes:
    def test_reparents_mrna(self, db_from_string):
        db = db_from_string(
            "chr1\ttest\tgene\t100\t500\t.\t+\t.\tID=g1\n"
            "chr1\ttest\tmRNA\t100\t500\t.\t+\t.\tID=m1;Parent=g1\n"
            "chr1\ttest\tgene\t200\t500\t.\t+\t.\tID=g2\n"
            "chr1\ttest\tmRNA\t200\t500\t.\t+\t.\tID=m2;Parent=g2\n"
        )
        merged = merge_fully_overlapping_transcript_genes(db)
        genes = list(merged.features_of_type("gene"))
        mrnas = list(merged.features_of_type("mRNA"))

        assert len(genes) == 1
        assert genes[0].id == "g1"
        assert len(mrnas) == 2
        # m2 should now have g1 as parent
        m2 = merged["m2"]
        assert m2.attributes["Parent"] == ["g1"]


# ---------------------------------------------------------------------------
# Non-overlapping genes
# ---------------------------------------------------------------------------


class TestFindNonoverlapping:
    def test_finds_non_overlapping_with_orf(self, db_from_string):
        db1 = db_from_string(
            "chr1\ttest\tgene\t100\t300\t.\t+\t.\tID=g1\n"
            "chr1\ttest\tmRNA\t100\t300\t.\t+\t.\tID=m1;Parent=g1\n"
            "chr1\ttest\tCDS\t100\t300\t.\t+\t0\tID=c1;Parent=m1\n"
            "chr1\ttest\tgene\t800\t1000\t.\t+\t.\tID=g2\n"
            "chr1\ttest\tmRNA\t800\t1000\t.\t+\t.\tID=m2;Parent=g2\n"
            "chr1\ttest\tCDS\t800\t1000\t.\t+\t0\tID=c2;Parent=m2\n"
        )
        db2 = db_from_string(
            "chr1\ttest\tgene\t100\t300\t.\t+\t.\tID=t1\n"
        )
        result = find_nonoverlapping_genes(db1, db2)
        gene_ids = [f.id for f in result if f.featuretype == "gene"]
        # g1 overlaps t1, g2 does not
        assert "g2" in gene_ids
        assert "g1" not in gene_ids

    def test_skips_genes_without_cds(self, db_from_string):
        db1 = db_from_string(
            "chr1\ttest\tgene\t800\t1000\t.\t+\t.\tID=g1\n"
            "chr1\ttest\tmRNA\t800\t1000\t.\t+\t.\tID=m1;Parent=g1\n"
            "chr1\ttest\texon\t800\t1000\t.\t+\t.\tID=e1;Parent=m1\n"
        )
        db2 = db_from_string(
            "chr1\ttest\tgene\t100\t300\t.\t+\t.\tID=t1\n"
        )
        result = find_nonoverlapping_genes(db1, db2)
        assert len([f for f in result if f.featuretype == "gene"]) == 0

    def test_long_spanning_target_engulfs_source(self, db_from_string):
        """A long-spanning target before short non-overlapping ones must
        still be detected as overlapping. Exercises the prefix-max-end
        bound used by the index."""
        db_source = db_from_string(
            "chr1\ttest\tgene\t5000\t5100\t.\t+\t.\tID=src\n"
            "chr1\ttest\tmRNA\t5000\t5100\t.\t+\t.\tID=m;Parent=src\n"
            "chr1\ttest\tCDS\t5000\t5100\t.\t+\t0\tID=c;Parent=m\n"
        )
        # Long-spanning gene starts well before src and engulfs it.
        # Several short genes with start > src.end serve as decoys.
        db_target = db_from_string(
            "chr1\ttest\tgene\t100\t9000\t.\t+\t.\tID=engulf\n"
            "chr1\ttest\tgene\t6000\t6100\t.\t+\t.\tID=after1\n"
            "chr1\ttest\tgene\t7000\t7100\t.\t+\t.\tID=after2\n"
        )
        result = find_nonoverlapping_genes(db_source, db_target)
        assert [f.id for f in result if f.featuretype == "gene"] == []


# ---------------------------------------------------------------------------
# Concordant models
# ---------------------------------------------------------------------------


class TestFindConcordantModels:
    def test_matching_single_cds_models(self, tmp_path):
        gff1 = tmp_path / "a.gff3"
        gff2 = tmp_path / "b.gff3"

        gff1.write_text(
            "chr1\tpred\tgene\t100\t400\t.\t+\t.\tID=ga\n"
            "chr1\tpred\tmRNA\t100\t400\t.\t+\t.\tID=ma;Parent=ga\n"
            "chr1\tpred\tCDS\t100\t400\t.\t+\t0\tID=ca;Parent=ma\n"
        )
        gff2.write_text(
            "chr1\taln\tgene\t105\t395\t.\t+\t.\tID=gb\n"
            "chr1\taln\tmRNA\t105\t395\t.\t+\t.\tID=mb;Parent=gb\n"
            "chr1\taln\tCDS\t105\t395\t.\t+\t0\tID=cb;Parent=mb\n"
        )

        result = find_concordant_models(gff1, gff2)
        genes = list(result.features_of_type("gene"))
        assert len(genes) == 1
        assert genes[0].id == "ga"

    def test_different_cds_count_not_concordant(self, tmp_path):
        gff1 = tmp_path / "a.gff3"
        gff2 = tmp_path / "b.gff3"

        gff1.write_text(
            "chr1\tpred\tgene\t100\t800\t.\t+\t.\tID=ga\n"
            "chr1\tpred\tmRNA\t100\t800\t.\t+\t.\tID=ma;Parent=ga\n"
            "chr1\tpred\tCDS\t100\t300\t.\t+\t0\tID=ca1;Parent=ma\n"
            "chr1\tpred\tCDS\t500\t800\t.\t+\t1\tID=ca2;Parent=ma\n"
        )
        gff2.write_text(
            "chr1\taln\tgene\t100\t800\t.\t+\t.\tID=gb\n"
            "chr1\taln\tmRNA\t100\t800\t.\t+\t.\tID=mb;Parent=gb\n"
            "chr1\taln\tCDS\t100\t800\t.\t+\t0\tID=cb;Parent=mb\n"
        )

        result = find_concordant_models(gff1, gff2)
        assert len(list(result.features_of_type("gene"))) == 0

    def test_no_overlap_not_concordant(self, tmp_path):
        gff1 = tmp_path / "a.gff3"
        gff2 = tmp_path / "b.gff3"

        gff1.write_text(
            "chr1\tpred\tgene\t100\t300\t.\t+\t.\tID=ga\n"
            "chr1\tpred\tmRNA\t100\t300\t.\t+\t.\tID=ma;Parent=ga\n"
            "chr1\tpred\tCDS\t100\t300\t.\t+\t0\tID=ca;Parent=ma\n"
        )
        gff2.write_text(
            "chr1\taln\tgene\t5000\t6000\t.\t+\t.\tID=gb\n"
            "chr1\taln\tmRNA\t5000\t6000\t.\t+\t.\tID=mb;Parent=gb\n"
            "chr1\taln\tCDS\t5000\t6000\t.\t+\t0\tID=cb;Parent=mb\n"
        )

        result = find_concordant_models(gff1, gff2)
        assert len(list(result.features_of_type("gene"))) == 0


# ---------------------------------------------------------------------------
# Non-redundant combination
# ---------------------------------------------------------------------------


class TestCombineNonredundant:
    def test_deduplicates_by_id(self, db_from_string):
        db1 = db_from_string(
            "chr1\ttest\tgene\t100\t300\t.\t+\t.\tID=g1\n"
            "chr1\ttest\tmRNA\t100\t300\t.\t+\t.\tID=m1;Parent=g1\n"
        )
        db2 = db_from_string(
            "chr1\ttest\tgene\t100\t300\t.\t+\t.\tID=g1\n"
            "chr1\ttest\tmRNA\t100\t300\t.\t+\t.\tID=m1;Parent=g1\n"
            "chr1\ttest\tgene\t500\t700\t.\t+\t.\tID=g2\n"
            "chr1\ttest\tmRNA\t500\t700\t.\t+\t.\tID=m2;Parent=g2\n"
        )
        result = combine_nonredundant_models(db1, db2)
        genes = list(result.features_of_type("gene"))
        assert len(genes) == 2
        ids = {g.id for g in genes}
        assert ids == {"g1", "g2"}


# ---------------------------------------------------------------------------
# Training-set extraction
# ---------------------------------------------------------------------------


class TestExtractSupportedModels:

    def _write_simple_pair(self, dirp, name, gene_id, start=100, end=400):
        path = dirp / name
        path.write_text(
            f"chr1\tpred\tgene\t{start}\t{end}\t.\t+\t.\tID={gene_id}\n"
            f"chr1\tpred\tmRNA\t{start}\t{end}\t.\t+\t.\tID=m_{gene_id};Parent={gene_id}\n"
            f"chr1\tpred\tCDS\t{start}\t{end}\t.\t+\t0\tID=c_{gene_id};Parent=m_{gene_id}\n"
        )
        return path

    def test_three_paths_parses_each_source_once(self, tmp_path, monkeypatch):
        """The three concordance pairs share one parse per source file
        (3 file parses, not one parse per side of every pair: 6)."""
        from eukan.gff import concordance

        parses: Counter[str] = Counter()
        real_create = concordance.create_gff_db

        def counting(source, **kwargs):
            if isinstance(source, (str, Path)):
                parses[Path(source).name] += 1
            return real_create(source, **kwargs)

        monkeypatch.setattr(concordance, "create_gff_db", counting)
        a = self._write_simple_pair(tmp_path, "a.gff3", "ga")
        b = self._write_simple_pair(tmp_path, "b.gff3", "ga")
        c = self._write_simple_pair(tmp_path, "c.gff3", "ga")

        result = extract_supported_models(a, b, c, output_dir=tmp_path)

        assert [parses[n] for n in ("a.gff3", "b.gff3", "c.gff3")] == [1, 1, 1]
        assert [g.id for g in result.features_of_type("gene")] == ["ga"]
        assert (tmp_path / "training_set.gff3").exists()

    def test_three_paths_near_identical_locus_kept_once(self, tmp_path):
        """The pairs return different sources' genes for one locus within
        the boundary tolerance; the training set keeps it once."""
        # (a,b) and (a,c) return a's ga; (c,b) returns c's gc, because
        # find_concordant_models compares structure, not IDs.
        a = self._write_simple_pair(tmp_path, "a.gff3", "ga", start=100, end=400)
        b = self._write_simple_pair(tmp_path, "b.gff3", "ga", start=105, end=395)
        # c shares ga's locus with a, but under a different ID.
        c_path = tmp_path / "c.gff3"
        c_path.write_text(
            "chr1\tpred\tgene\t102\t398\t.\t+\t.\tID=gc\n"
            "chr1\tpred\tmRNA\t102\t398\t.\t+\t.\tID=mc;Parent=gc\n"
            "chr1\tpred\tCDS\t102\t398\t.\t+\t0\tID=cc;Parent=mc\n"
        )

        result = extract_supported_models(a, b, c_path, output_dir=tmp_path)

        # gc is the same locus as ga: the first source's model is kept, once.
        genes = list(result.features_of_type("gene"))
        assert [g.id for g in genes] == ["ga"]

    @staticmethod
    def _write_two_cds(dirp, name, gene_id, offset=0):
        """One 2-CDS gene at 100-900 (+offset) on chr1."""
        s = 100 + offset
        path = dirp / name
        path.write_text(
            f"chr1\tpred\tgene\t{s}\t{s + 800}\t.\t+\t.\tID={gene_id}\n"
            f"chr1\tpred\tmRNA\t{s}\t{s + 800}\t.\t+\t.\tID=m_{gene_id};Parent={gene_id}\n"
            f"chr1\tpred\tCDS\t{s}\t{s + 200}\t.\t+\t0\tID=c1_{gene_id};Parent=m_{gene_id}\n"
            f"chr1\tpred\tCDS\t{s + 400}\t{s + 800}\t.\t+\t0\tID=c2_{gene_id};Parent=m_{gene_id}\n"
        )
        return path

    def test_locus_concordant_in_all_three_sources_kept_once(self, tmp_path, caplog):
        """One identical locus in all three sources, under each source's own
        IDs: (genemark,prot) and (genemark,orfs) return gm1, (orfs,prot)
        returns orf1. The training set must hold one gene, not gm1 + orf1."""
        import logging
        gm = self._write_two_cds(tmp_path, "genemark.gff3", "gm1")
        prot = self._write_two_cds(tmp_path, "prot.gff3", "pr1")
        orfs = self._write_two_cds(tmp_path, "orfs.gff3", "orf1")

        with caplog.at_level(logging.WARNING, logger="eukan.gff.concordance"):
            result = extract_supported_models(gm, prot, orfs, output_dir=tmp_path)

        assert [g.id for g in result.features_of_type("gene")] == ["gm1"]
        written = (tmp_path / "training_set.gff3").read_text()
        assert written.count("\tgene\t") == 1
        assert "Only 1 concordant" in caplog.text

    def test_locus_supported_only_by_third_pair_is_kept(self, tmp_path):
        """A locus where only the third and second sources agree (no
        first-source model there) still enters the training set."""
        gm = self._write_two_cds(tmp_path, "genemark.gff3", "gm1")
        prot = self._write_two_cds(tmp_path, "prot.gff3", "pr1", offset=5000)
        orfs = self._write_two_cds(tmp_path, "orfs.gff3", "orf1", offset=5000)

        result = extract_supported_models(gm, prot, orfs, output_dir=tmp_path)

        assert [g.id for g in result.features_of_type("gene")] == ["orf1"]

    def test_two_paths_unchanged(self, tmp_path):
        """The 2-path branch: one concordance pass."""
        a = self._write_simple_pair(tmp_path, "a.gff3", "ga")
        b = self._write_simple_pair(tmp_path, "b.gff3", "ga")

        result = extract_supported_models(a, b, output_dir=tmp_path)
        genes = list(result.features_of_type("gene"))
        assert len(genes) == 1
        assert genes[0].id == "ga"

    def test_no_concordant_models_returns_empty_db(self, tmp_path):
        """Zero concordant models: an empty training_set.gff3 and an empty
        DB (re-reading the empty file raised gffutils' EmptyInputError)."""
        a = self._write_simple_pair(tmp_path, "a.gff3", "ga")
        b = self._write_simple_pair(tmp_path, "b.gff3", "gb", start=5000, end=5300)

        result = extract_supported_models(a, b, output_dir=tmp_path)

        assert result.count_features_of_type("gene") == 0
        assert (tmp_path / "training_set.gff3").read_text() == ""

    def test_three_paths_logs_per_source_counts(self, tmp_path, caplog):
        """3-way concordance: per-source gene counts always go to INFO,
        labelled by file stem so users see ``genemark`` / ``spaln`` /
        ``transcript_orfs`` in the eukan annotate output."""
        import logging
        a = self._write_simple_pair(tmp_path, "genemark.gff3", "ga")
        b = self._write_simple_pair(tmp_path, "spaln.gff3", "ga")
        c = self._write_simple_pair(tmp_path, "transcript_orfs.gff3", "ga")

        with caplog.at_level(logging.INFO, logger="eukan.gff.concordance"):
            extract_supported_models(a, b, c, output_dir=tmp_path)

        assert "Evidence model counts" in caplog.text
        assert "genemark=1" in caplog.text
        assert "spaln=1" in caplog.text
        assert "transcript_orfs=1" in caplog.text

    def test_three_paths_warns_on_weak_concordance(self, tmp_path, caplog):
        """Below the WEAK_CONCORDANCE_THRESHOLD a WARNING fires for 3-way input."""
        import logging
        a = self._write_simple_pair(tmp_path, "genemark.gff3", "ga")
        b = self._write_simple_pair(tmp_path, "spaln.gff3", "ga")
        c = self._write_simple_pair(tmp_path, "transcript_orfs.gff3", "ga")

        with caplog.at_level(logging.WARNING, logger="eukan.gff.concordance"):
            extract_supported_models(a, b, c, output_dir=tmp_path)

        assert "weak concordance" in caplog.text.lower()
        assert "Only 1 concordant" in caplog.text

    def test_three_paths_relabels_prot_as_spaln(self, tmp_path, caplog):
        """``prot.gff3`` is the spaln output; the log should call it ``spaln``."""
        import logging
        a = self._write_simple_pair(tmp_path, "genemark.gff3", "ga")
        b = self._write_simple_pair(tmp_path, "prot.gff3", "ga")
        c = self._write_simple_pair(tmp_path, "transcript_orfs.gff3", "ga")

        with caplog.at_level(logging.INFO, logger="eukan.gff.concordance"):
            extract_supported_models(a, b, c, output_dir=tmp_path)

        assert "spaln=1" in caplog.text
        assert "prot=" not in caplog.text

    def test_snap_triple_uses_augustus_not_genemark(self, tmp_path, caplog):
        """SNAP trains from AUGUSTUS predictions, not GeneMark — the log
        should reflect that the first source is augustus."""
        import logging
        a = self._write_simple_pair(tmp_path, "augustus.gff3", "ga")
        b = self._write_simple_pair(tmp_path, "prot.gff3", "ga")
        c = self._write_simple_pair(tmp_path, "transcript_orfs.gff3", "ga")

        with caplog.at_level(logging.INFO, logger="eukan.gff.concordance"):
            extract_supported_models(a, b, c, output_dir=tmp_path)

        assert "augustus=1" in caplog.text
        assert "spaln=1" in caplog.text
        assert "transcript_orfs=1" in caplog.text
        assert "genemark" not in caplog.text

    def test_two_paths_does_not_warn_on_low_concordance(self, tmp_path, caplog):
        """The warning only applies to the 3-way path."""
        import logging
        a = self._write_simple_pair(tmp_path, "genemark.gff3", "ga")
        b = self._write_simple_pair(tmp_path, "spaln.gff3", "ga")

        with caplog.at_level(logging.WARNING, logger="eukan.gff.concordance"):
            extract_supported_models(a, b, output_dir=tmp_path)

        assert "weak concordance" not in caplog.text.lower()
