"""Tests for eukan.gff.intersecter — interval operations with pure gffutils."""

import gffutils

from eukan.gff.intersecter import (
    _cds_boundaries_match,
    _features_overlap,
    _find_overlapping_genes,
    combine_nonredundant_models,
    extract_supported_models,
    find_concordant_models,
    find_nonoverlapping_genes,
    merge_fully_overlapping_transcript_genes,
)


def _db(gff_string: str) -> gffutils.FeatureDB:
    return gffutils.create_db(
        gff_string, ":memory:", from_string=True, merge_strategy="create_unique"
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


class TestFeaturesOverlap:
    def test_overlapping(self):
        a = _feature("chr1", 100, 500, "+")
        b = _feature("chr1", 400, 800, "+")
        assert _features_overlap(a, b)

    def test_non_overlapping(self):
        a = _feature("chr1", 100, 300, "+")
        b = _feature("chr1", 400, 800, "+")
        assert not _features_overlap(a, b)

    def test_different_strand(self):
        a = _feature("chr1", 100, 500, "+")
        b = _feature("chr1", 100, 500, "-")
        assert not _features_overlap(a, b)

    def test_different_chrom(self):
        a = _feature("chr1", 100, 500, "+")
        b = _feature("chr2", 100, 500, "+")
        assert not _features_overlap(a, b)

    def test_adjacent_not_overlapping(self):
        a = _feature("chr1", 100, 200, "+")
        b = _feature("chr1", 201, 300, "+")
        assert not _features_overlap(a, b)

    def test_touching(self):
        a = _feature("chr1", 100, 200, "+")
        b = _feature("chr1", 200, 300, "+")
        assert _features_overlap(a, b)


class TestFindOverlappingGenes:
    def test_two_overlapping(self):
        db = _db(
            "chr1\ttest\tgene\t100\t500\t.\t+\t.\tID=g1\n"
            "chr1\ttest\tgene\t300\t700\t.\t+\t.\tID=g2\n"
        )
        clusters = _find_overlapping_genes(db)
        assert len(clusters) == 1
        assert len(clusters[0]) == 2

    def test_no_overlap(self):
        db = _db(
            "chr1\ttest\tgene\t100\t200\t.\t+\t.\tID=g1\n"
            "chr1\ttest\tgene\t500\t700\t.\t+\t.\tID=g2\n"
        )
        clusters = _find_overlapping_genes(db)
        assert len(clusters) == 0

    def test_different_strands_not_merged(self):
        db = _db(
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
    def test_reparents_mrna(self):
        db = _db(
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
    def test_finds_non_overlapping_with_orf(self):
        db1 = _db(
            "chr1\ttest\tgene\t100\t300\t.\t+\t.\tID=g1\n"
            "chr1\ttest\tmRNA\t100\t300\t.\t+\t.\tID=m1;Parent=g1\n"
            "chr1\ttest\tCDS\t100\t300\t.\t+\t0\tID=c1;Parent=m1\n"
            "chr1\ttest\tgene\t800\t1000\t.\t+\t.\tID=g2\n"
            "chr1\ttest\tmRNA\t800\t1000\t.\t+\t.\tID=m2;Parent=g2\n"
            "chr1\ttest\tCDS\t800\t1000\t.\t+\t0\tID=c2;Parent=m2\n"
        )
        db2 = _db(
            "chr1\ttest\tgene\t100\t300\t.\t+\t.\tID=t1\n"
        )
        result = find_nonoverlapping_genes(db1, db2)
        gene_ids = [f.id for f in result if f.featuretype == "gene"]
        # g1 overlaps t1, g2 does not
        assert "g2" in gene_ids
        assert "g1" not in gene_ids

    def test_skips_genes_without_cds(self):
        db1 = _db(
            "chr1\ttest\tgene\t800\t1000\t.\t+\t.\tID=g1\n"
            "chr1\ttest\tmRNA\t800\t1000\t.\t+\t.\tID=m1;Parent=g1\n"
            "chr1\ttest\texon\t800\t1000\t.\t+\t.\tID=e1;Parent=m1\n"
        )
        db2 = _db(
            "chr1\ttest\tgene\t100\t300\t.\t+\t.\tID=t1\n"
        )
        result = find_nonoverlapping_genes(db1, db2)
        assert len([f for f in result if f.featuretype == "gene"]) == 0

    def test_long_spanning_target_engulfs_source(self):
        """A long-spanning target before short non-overlapping ones must
        still be detected as overlapping. Exercises the prefix-max-end
        bound used by the index."""
        db_source = _db(
            "chr1\ttest\tgene\t5000\t5100\t.\t+\t.\tID=src\n"
            "chr1\ttest\tmRNA\t5000\t5100\t.\t+\t.\tID=m;Parent=src\n"
            "chr1\ttest\tCDS\t5000\t5100\t.\t+\t0\tID=c;Parent=m\n"
        )
        # Long-spanning gene starts well before src and engulfs it.
        # Several short genes with start > src.end serve as decoys.
        db_target = _db(
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
    def test_deduplicates_by_id(self):
        db1 = _db(
            "chr1\ttest\tgene\t100\t300\t.\t+\t.\tID=g1\n"
            "chr1\ttest\tmRNA\t100\t300\t.\t+\t.\tID=m1;Parent=g1\n"
        )
        db2 = _db(
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
# Training-set extraction (3-path concurrent branch)
# ---------------------------------------------------------------------------


class TestExtractSupportedModels:
    """Regression: the 3-path branch fans concordance passes across a
    ``ThreadPoolExecutor``.  Returning ``gffutils.FeatureDB`` from worker
    threads broke at the next main-thread query with
    ``sqlite3.ProgrammingError: SQLite objects created in a thread can
    only be used in that same thread``.  The fix returns Feature lists
    from workers and builds dbs on the main thread.
    """

    def _write_simple_pair(self, dirp, name, gene_id, start=100, end=400):
        path = dirp / name
        path.write_text(
            f"chr1\tpred\tgene\t{start}\t{end}\t.\t+\t.\tID={gene_id}\n"
            f"chr1\tpred\tmRNA\t{start}\t{end}\t.\t+\t.\tID=m_{gene_id};Parent={gene_id}\n"
            f"chr1\tpred\tCDS\t{start}\t{end}\t.\t+\t0\tID=c_{gene_id};Parent=m_{gene_id}\n"
        )
        return path

    def test_three_paths_does_not_leak_db_across_threads(self, tmp_path):
        """Calling extract_supported_models with 3 paths must produce a
        usable FeatureDB on the main thread, with no SQLite cross-thread
        error."""
        a = self._write_simple_pair(tmp_path, "a.gff3", "ga")
        b = self._write_simple_pair(tmp_path, "b.gff3", "ga")
        c = self._write_simple_pair(tmp_path, "c.gff3", "ga")

        result = extract_supported_models(a, b, c, output_dir=tmp_path)

        # The query that previously raised ProgrammingError in the failing run.
        genes = list(result.features_of_type("gene"))
        assert len(genes) == 1
        assert genes[0].id == "ga"
        assert (tmp_path / "training_set.gff3").exists()

    def test_three_paths_partial_concordance_queryable_main_thread(self, tmp_path):
        """When some pairs concord and others don't, the per-worker
        Feature lists vary in length — but the assembled db must still be
        queryable on the main thread without a cross-thread error."""
        # ga concordant in (a,b) and (a,c); (b,c) finds gc only when c==b's gene.
        # We arrange: a==b==c on ga; (a,c) concord ga; (c,b) finds nothing.
        a = self._write_simple_pair(tmp_path, "a.gff3", "ga", start=100, end=400)
        b = self._write_simple_pair(tmp_path, "b.gff3", "ga", start=105, end=395)
        # c shares ga's locus with a, but different ID — concord (a,c) finds ga,
        # (c,b) compares c's gc-named gene against b's ga, which won't match by ID
        # but find_concordant_models doesn't compare IDs, only structure.
        c_path = tmp_path / "c.gff3"
        c_path.write_text(
            "chr1\tpred\tgene\t102\t398\t.\t+\t.\tID=gc\n"
            "chr1\tpred\tmRNA\t102\t398\t.\t+\t.\tID=mc;Parent=gc\n"
            "chr1\tpred\tCDS\t102\t398\t.\t+\t0\tID=cc;Parent=mc\n"
        )

        result = extract_supported_models(a, b, c_path, output_dir=tmp_path)

        # Must not raise — and must be queryable on the main thread.
        genes = list(result.features_of_type("gene"))
        assert {g.id for g in genes} <= {"ga", "gc"}
        assert len(genes) >= 1

    def test_two_paths_unchanged(self, tmp_path):
        """The 2-path branch is sequential; verify it still works."""
        a = self._write_simple_pair(tmp_path, "a.gff3", "ga")
        b = self._write_simple_pair(tmp_path, "b.gff3", "ga")

        result = extract_supported_models(a, b, output_dir=tmp_path)
        genes = list(result.features_of_type("gene"))
        assert len(genes) == 1
        assert genes[0].id == "ga"
