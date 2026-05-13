"""Tests for KOfam parsing, EC extraction, and GFF3 attribute emission."""

from __future__ import annotations

from pathlib import Path

from eukan.functional.kofam import (
    KOEntry,
    extract_ec_numbers,
    parse_ko_list,
)
from eukan.functional.search import _add_func_info, annotate_gff3
from eukan.gff import create_gff_db

# ---------------------------------------------------------------------------
# parse_ko_list
# ---------------------------------------------------------------------------


_HEADER = "knum\tthreshold\tscore_type\tprofile_type\tF-measure\tnseq\tnseq_used\talen\tmlen\teff_nseq\tre/pos\tdefinition\n"


def _write_ko_list(path: Path, rows: list[str]) -> Path:
    path.write_text(_HEADER + "".join(r if r.endswith("\n") else r + "\n" for r in rows))
    return path


class TestParseKoList:
    def test_basic_row(self, tmp_path: Path):
        ko_list = _write_ko_list(tmp_path / "ko_list", [
            "K28880\t443.33\tfull\tall\t1.000000\t812\t773\t1077\t529\t2.48\t0.590\t"
            "peptidoglycan beta-N-acetylmuramidase [EC:3.2.1.92]",
        ])
        entries = parse_ko_list(ko_list)
        assert "K28880" in entries
        e = entries["K28880"]
        assert e == KOEntry(
            k_number="K28880",
            threshold=443.33,
            score_type="full",
            definition="peptidoglycan beta-N-acetylmuramidase [EC:3.2.1.92]",
        )

    def test_missing_threshold_becomes_none(self, tmp_path: Path):
        ko_list = _write_ko_list(tmp_path / "ko_list", [
            "K99999\t-\tfull\tall\t-\t-\t-\t-\t-\t-\t-\tplaceholder KO without threshold",
        ])
        e = parse_ko_list(ko_list)["K99999"]
        assert e.threshold is None

    def test_domain_score_type(self, tmp_path: Path):
        ko_list = _write_ko_list(tmp_path / "ko_list", [
            "K28881\t275.30\tdomain\tall\t0.996650\t308\t300\t1189\t387\t1.94\t0.590\t"
            "potassium large conductance calcium-activated channel auxiliary gamma subunit LRRC26",
        ])
        e = parse_ko_list(ko_list)["K28881"]
        assert e.score_type == "domain"
        assert e.threshold == 275.30

    def test_skips_malformed_rows(self, tmp_path: Path):
        ko_list = _write_ko_list(tmp_path / "ko_list", [
            "K28880\t443.33\tfull\tall\t1.000000\t812\t773\t1077\t529\t2.48\t0.590\tnormal row",
            "garbage line with too few columns",
            "K12345\t100.0\tfull\tall\t1.0\t10\t10\t100\t50\t1.0\t0.5\tanother good row",
        ])
        entries = parse_ko_list(ko_list)
        assert set(entries) == {"K28880", "K12345"}


# ---------------------------------------------------------------------------
# extract_ec_numbers
# ---------------------------------------------------------------------------


class TestExtractEcNumbers:
    def test_single_ec(self):
        cleaned, ecs = extract_ec_numbers(
            "peptidoglycan beta-N-acetylmuramidase [EC:3.2.1.92]"
        )
        assert cleaned == "peptidoglycan beta-N-acetylmuramidase"
        assert ecs == ["3.2.1.92"]

    def test_no_ec(self):
        cleaned, ecs = extract_ec_numbers(
            "potassium large conductance calcium-activated channel auxiliary gamma subunit LRRC26"
        )
        assert cleaned.endswith("LRRC26")
        assert ecs == []

    def test_multiple_ec_space_separated(self):
        # KEGG sometimes packs multiple ECs inside one bracket
        cleaned, ecs = extract_ec_numbers(
            "bifunctional widget [EC:1.1.1.1 2.7.1.1]"
        )
        assert cleaned == "bifunctional widget"
        assert ecs == ["1.1.1.1", "2.7.1.1"]

    def test_multiple_ec_comma_separated(self):
        cleaned, ecs = extract_ec_numbers(
            "bifunctional widget [EC:1.1.1.1,2.7.1.1]"
        )
        assert ecs == ["1.1.1.1", "2.7.1.1"]

    def test_partial_ec_with_dash_kept(self):
        cleaned, ecs = extract_ec_numbers(
            "putative dehydrogenase [EC:1.1.1.-]"
        )
        assert cleaned == "putative dehydrogenase"
        assert ecs == ["1.1.1.-"]

    def test_preliminary_ec_with_n_digit_kept(self):
        # Preliminary EC like 'n6' is allowed by NCBI ec_number qualifier
        cleaned, ecs = extract_ec_numbers("preliminary enzyme [EC:3.4.21.n6]")
        assert ecs == ["3.4.21.n6"]

    def test_placeholder_n_n_n_n_dropped(self):
        # The KEGG 'no EC assigned' placeholder must not become an NCBI EC qualifier
        cleaned, ecs = extract_ec_numbers("uncharacterized protein [EC:n.n.n.n]")
        assert cleaned == "uncharacterized protein"
        assert ecs == []

    def test_dedupes_repeated_codes(self):
        _cleaned, ecs = extract_ec_numbers("foo [EC:1.1.1.1] bar [EC:1.1.1.1]")
        assert ecs == ["1.1.1.1"]


# ---------------------------------------------------------------------------
# KOfam attribute emission via _add_func_info / _apply_kofam_hits
# ---------------------------------------------------------------------------


def _minimal_gff(tmp_path: Path) -> Path:
    gff = tmp_path / "in.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\teukan\tgene\t1\t300\t.\t+\t.\tID=g1\n"
        "chr1\teukan\tmRNA\t1\t300\t.\t+\t.\tID=m1;Parent=g1\n"
        "chr1\teukan\tCDS\t1\t300\t.\t+\t0\tID=c1;Parent=m1\n"
    )
    return gff


class TestApplyKofamHits:
    def test_above_threshold_emits_product_ec_dbxref_inference(self, tmp_path: Path):
        db = create_gff_db(_minimal_gff(tmp_path))
        kofam_res = {
            "m1": {
                "K28880": {
                    "description": "peptidoglycan beta-N-acetylmuramidase",
                    "evalue": 1e-30,
                    "score": 500.0,
                    "threshold": 443.33,
                    "score_type": "full",
                    "above_threshold": True,
                    "ec_numbers": ["3.2.1.92"],
                }
            }
        }
        feats = list(_add_func_info(db, kofam_res, {}, homology_db="kofam"))
        mrna = next(f for f in feats if f.featuretype == "mRNA")
        assert mrna.attributes["product"] == ["peptidoglycan beta-N-acetylmuramidase"]
        assert mrna.attributes["ec_number"] == ["3.2.1.92"]
        assert mrna.attributes["Dbxref"] == ["KEGG:K28880"]
        assert mrna.attributes["inference"] == ["protein motif:KOFAM:K28880"]

    def test_below_threshold_falls_back_to_hypothetical(self, tmp_path: Path):
        db = create_gff_db(_minimal_gff(tmp_path))
        kofam_res = {
            "m1": {
                "K28880": {
                    "description": "peptidoglycan beta-N-acetylmuramidase",
                    "evalue": 1e-5,
                    "score": 300.0,
                    "threshold": 443.33,
                    "score_type": "full",
                    "above_threshold": False,
                    "ec_numbers": ["3.2.1.92"],
                }
            }
        }
        feats = list(_add_func_info(db, kofam_res, {}, homology_db="kofam"))
        mrna = next(f for f in feats if f.featuretype == "mRNA")
        assert mrna.attributes["product"] == ["hypothetical protein"]
        # No KOfam attributes emitted when nothing is above threshold
        assert "ec_number" not in mrna.attributes
        assert "Dbxref" not in mrna.attributes

    def test_multiple_above_threshold_top_score_drives_product(self, tmp_path: Path):
        db = create_gff_db(_minimal_gff(tmp_path))
        kofam_res = {
            "m1": {
                "K00001": {
                    "description": "weaker hit",
                    "evalue": 1e-10, "score": 250.0, "threshold": 200.0,
                    "score_type": "full", "above_threshold": True,
                    "ec_numbers": [],
                },
                "K00002": {
                    "description": "stronger hit",
                    "evalue": 1e-50, "score": 800.0, "threshold": 200.0,
                    "score_type": "full", "above_threshold": True,
                    "ec_numbers": ["1.1.1.1"],
                },
            }
        }
        feats = list(_add_func_info(db, kofam_res, {}, homology_db="kofam"))
        mrna = next(f for f in feats if f.featuretype == "mRNA")
        assert mrna.attributes["product"] == ["stronger hit"]
        # Both KOs make it into Dbxref / inference, top-score first
        assert mrna.attributes["Dbxref"] == ["KEGG:K00002", "KEGG:K00001"]
        assert mrna.attributes["inference"] == [
            "protein motif:KOFAM:K00002",
            "protein motif:KOFAM:K00001",
        ]
        assert mrna.attributes["ec_number"] == ["1.1.1.1"]

    def test_dedupes_ec_across_hits(self, tmp_path: Path):
        db = create_gff_db(_minimal_gff(tmp_path))
        kofam_res = {
            "m1": {
                "K00001": {
                    "description": "iso A", "evalue": 1e-50, "score": 500.0,
                    "threshold": 100.0, "score_type": "full",
                    "above_threshold": True, "ec_numbers": ["1.1.1.1"],
                },
                "K00002": {
                    "description": "iso B", "evalue": 1e-40, "score": 400.0,
                    "threshold": 100.0, "score_type": "full",
                    "above_threshold": True, "ec_numbers": ["1.1.1.1", "2.7.1.1"],
                },
            }
        }
        feats = list(_add_func_info(db, kofam_res, {}, homology_db="kofam"))
        mrna = next(f for f in feats if f.featuretype == "mRNA")
        assert mrna.attributes["ec_number"] == ["1.1.1.1", "2.7.1.1"]


class TestKofamPlusPfamCombined:
    def test_kofam_and_pfam_inferences_both_appear(self, tmp_path: Path):
        db = create_gff_db(_minimal_gff(tmp_path))
        kofam_res = {
            "m1": {
                "K28880": {
                    "description": "peptidoglycan beta-N-acetylmuramidase",
                    "evalue": 1e-30, "score": 500.0, "threshold": 443.33,
                    "score_type": "full", "above_threshold": True,
                    "ec_numbers": ["3.2.1.92"],
                }
            }
        }
        pfam_res = {
            "m1": {
                "Glyco_hydro": {
                    "description": "Glycoside hydrolase domain",
                    "evalue": 1e-20,
                    "accession": "PF00710.15",
                }
            }
        }
        feats = list(_add_func_info(db, kofam_res, pfam_res, homology_db="kofam"))
        mrna = next(f for f in feats if f.featuretype == "mRNA")
        assert mrna.attributes["inference"] == [
            "protein motif:KOFAM:K28880",
            "protein motif:PFAM:PF00710",
        ]


# ---------------------------------------------------------------------------
# Round-trip through annotate_gff3 → submission cleanup keeps ec_number
# ---------------------------------------------------------------------------


class TestEcNumberThroughCleanup:
    """ec_number must survive cleanup_gff3_for_submission untouched."""

    def test_ec_number_preserved(self, tmp_path: Path):
        from eukan.submission.cleanup import clean_gff3_for_submission

        gff = _minimal_gff(tmp_path)
        kofam_res = {
            "m1": {
                "K28880": {
                    "description": "peptidoglycan beta-N-acetylmuramidase",
                    "evalue": 1e-30, "score": 500.0, "threshold": 443.33,
                    "score_type": "full", "above_threshold": True,
                    "ec_numbers": ["3.2.1.92"],
                }
            }
        }
        annotated = annotate_gff3(gff, kofam_res, {}, output_dir=tmp_path, homology_db="kofam")
        cleaned = tmp_path / "cleaned.gff3"
        clean_gff3_for_submission(annotated, cleaned)
        text = cleaned.read_text()
        assert "ec_number=3.2.1.92" in text
        assert "Dbxref=KEGG:K28880" in text
