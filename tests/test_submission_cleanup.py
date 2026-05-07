"""Tests for submission/cleanup.py — GFF3 preprocessing for table2asn."""

from __future__ import annotations

from pathlib import Path

import pytest

from eukan.gff import create_gff_db
from eukan.submission.cleanup import (
    DEFAULT_PRODUCT,
    INFERENCE_CAP,
    _clean_product,
    _identify_drop_set,
    clean_gff3_for_submission,
)


# ---------------------------------------------------------------------------
# _clean_product
# ---------------------------------------------------------------------------


class TestCleanProduct:
    def test_strips_uniprot_metadata(self):
        cleaned, os_strip, frag = _clean_product(
            "Malic acid transport protein OS=Schizosaccharomyces pombe "
            "OX=4896 GN=mae1 PE=1 SV=1"
        )
        assert cleaned == "Malic acid transport protein"
        assert os_strip is True
        assert frag is False

    def test_strips_fragment_qualifier(self):
        cleaned, os_strip, frag = _clean_product(
            "ATP-dependent DNA helicase tlh1 (Fragment) OS=Foo OX=1 GN=x PE=3 SV=1"
        )
        assert cleaned == "ATP-dependent DNA helicase tlh1"
        assert os_strip is True
        assert frag is True

    def test_fragment_without_uniprot(self):
        cleaned, os_strip, frag = _clean_product("Some protein (Fragment)")
        assert cleaned == "Some protein"
        assert os_strip is False
        assert frag is True

    def test_clean_product_unchanged(self):
        cleaned, os_strip, frag = _clean_product("hypothetical protein")
        assert cleaned == "hypothetical protein"
        assert os_strip is False
        assert frag is False

    def test_empty_after_cleanup_falls_back(self):
        # A pathological input that strips down to whitespace
        cleaned, _, _ = _clean_product(" OS=org OX=1")
        assert cleaned == DEFAULT_PRODUCT


# ---------------------------------------------------------------------------
# _identify_drop_set
# ---------------------------------------------------------------------------


class TestIdentifyDropSet:
    def _write(self, tmp_path: Path, lines: list[str]) -> Path:
        gff = tmp_path / "in.gff3"
        gff.write_text("##gff-version 3\n" + "\n".join(lines) + "\n")
        return gff

    def test_drops_mrna_without_cds(self, tmp_path: Path):
        gff = self._write(tmp_path, [
            "chr1\teukan\tgene\t1\t100\t.\t+\t.\tID=g1",
            "chr1\teukan\tmRNA\t1\t100\t.\t+\t.\tID=m1;Parent=g1",
            "chr1\teukan\texon\t1\t100\t.\t+\t.\tID=e1;Parent=m1",
            # m2 has CDS, m1 does not
            "chr1\teukan\tmRNA\t1\t100\t.\t+\t.\tID=m2;Parent=g1",
            "chr1\teukan\texon\t1\t100\t.\t+\t.\tID=e2;Parent=m2",
            "chr1\teukan\tCDS\t10\t90\t.\t+\t0\tID=c2;Parent=m2",
        ])
        db = create_gff_db(gff)
        mrnas, genes = _identify_drop_set(db)
        assert mrnas == {"m1"}
        assert genes == set()

    def test_drops_orphan_gene_when_all_mrnas_lack_cds(self, tmp_path: Path):
        gff = self._write(tmp_path, [
            "chr1\teukan\tgene\t1\t100\t.\t+\t.\tID=g1",
            "chr1\teukan\tmRNA\t1\t100\t.\t+\t.\tID=m1;Parent=g1",
            "chr1\teukan\texon\t1\t100\t.\t+\t.\tID=e1;Parent=m1",
        ])
        db = create_gff_db(gff)
        mrnas, genes = _identify_drop_set(db)
        assert mrnas == {"m1"}
        assert genes == {"g1"}

    def test_keeps_everything_when_all_mrnas_have_cds(self, tmp_path: Path):
        gff = self._write(tmp_path, [
            "chr1\teukan\tgene\t1\t100\t.\t+\t.\tID=g1",
            "chr1\teukan\tmRNA\t1\t100\t.\t+\t.\tID=m1;Parent=g1",
            "chr1\teukan\tCDS\t10\t90\t.\t+\t0\tID=c1;Parent=m1",
        ])
        db = create_gff_db(gff)
        mrnas, genes = _identify_drop_set(db)
        assert mrnas == set()
        assert genes == set()


# ---------------------------------------------------------------------------
# clean_gff3_for_submission (integration)
# ---------------------------------------------------------------------------


class TestCleanGff3:
    @pytest.fixture
    def messy_gff(self, tmp_path: Path) -> Path:
        gff = tmp_path / "in.gff3"
        gff.write_text(
            "##gff-version 3\n"
            # Gene 1: kept (has CDS-bearing mRNA), product is UniProt-style
            "chr1\teukan\tgene\t1\t1000\t.\t+\t.\tID=g1\n"
            "chr1\teukan\tmRNA\t1\t1000\t.\t+\t.\tID=m1;Parent=g1;"
            "product=Malic acid transport protein OS%3DSchizosaccharomyces "
            "pombe OX%3D4896 GN%3Dmae1 PE%3D1 SV%3D1;"
            "inference=similar to AA sequence:UniProtKB:sp|A|B,"
            "protein motif:PFAM:X,protein motif:PFAM:Y,protein motif:PFAM:Z,"
            "protein motif:PFAM:W\n"
            "chr1\teukan\texon\t1\t1000\t.\t+\t.\tID=e1;Parent=m1\n"
            "chr1\teukan\tCDS\t100\t900\t.\t+\t0\tID=c1;Parent=m1\n"
            # m1b: same gene, no CDS — should be dropped, gene survives
            "chr1\teukan\tmRNA\t1\t1000\t.\t+\t.\tID=m1b;Parent=g1\n"
            "chr1\teukan\texon\t1\t1000\t.\t+\t.\tID=e1b;Parent=m1b\n"
            # Gene 2: only mRNA without CDS — gene + mRNA dropped
            "chr2\teukan\tgene\t1\t500\t.\t+\t.\tID=g2\n"
            "chr2\teukan\tmRNA\t1\t500\t.\t+\t.\tID=m2;Parent=g2\n"
            "chr2\teukan\texon\t1\t500\t.\t+\t.\tID=e2;Parent=m2\n"
            # Gene 3: has fragment in name
            "chr3\teukan\tgene\t1\t300\t.\t+\t.\tID=g3\n"
            "chr3\teukan\tmRNA\t1\t300\t.\t+\t.\tID=m3;Parent=g3;"
            "product=Helicase tlh1 (Fragment)\n"
            "chr3\teukan\texon\t1\t300\t.\t+\t.\tID=e3;Parent=m3\n"
            "chr3\teukan\tCDS\t10\t290\t.\t+\t0\tID=c3;Parent=m3\n"
        )
        return gff

    def test_full_cleanup(self, messy_gff: Path, tmp_path: Path):
        out = tmp_path / "cleaned.gff3"
        report = clean_gff3_for_submission(messy_gff, out)

        assert report.products_cleaned == 1
        assert report.fragments_stripped == 1
        assert report.inferences_capped == 1
        assert report.mrnas_dropped == 2  # m1b and m2
        assert report.genes_dropped == 1  # g2

        text = out.read_text()
        # UniProt cruft gone
        assert "OS=" not in text
        assert "GN=" not in text
        assert "Malic acid transport protein" in text
        # Fragment stripped
        assert "(Fragment)" not in text
        assert "Helicase tlh1" in text
        # Dropped features absent
        assert "ID=g2" not in text
        assert "ID=m2" not in text
        assert "ID=m1b" not in text
        # Surviving features present
        assert "ID=g1" in text
        assert "ID=g3" in text

    def test_inference_cap(self, messy_gff: Path, tmp_path: Path):
        out = tmp_path / "cleaned.gff3"
        clean_gff3_for_submission(messy_gff, out)
        cleaned_db = create_gff_db(out)
        m1 = cleaned_db["m1"]
        assert len(m1.attributes["inference"]) == INFERENCE_CAP

    def test_idempotent(self, messy_gff: Path, tmp_path: Path):
        first = tmp_path / "first.gff3"
        second = tmp_path / "second.gff3"
        clean_gff3_for_submission(messy_gff, first)
        report2 = clean_gff3_for_submission(first, second)
        # Running twice does nothing on the second pass
        assert report2.products_cleaned == 0
        assert report2.fragments_stripped == 0
        assert report2.inferences_capped == 0
        assert report2.mrnas_dropped == 0
        assert report2.genes_dropped == 0


# ---------------------------------------------------------------------------
# CleanupReport.summary
# ---------------------------------------------------------------------------


def test_report_summary():
    from eukan.submission.cleanup import CleanupReport

    r = CleanupReport(
        products_cleaned=10, fragments_stripped=2, inferences_capped=5,
        mrnas_dropped=3, genes_dropped=1,
    )
    s = r.summary()
    assert "10" in s and "products" in s
    assert "2" in s and "fragments" in s
    assert "5" in s and "inferences" in s
    assert "3" in s and "1" in s
