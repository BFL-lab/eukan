"""Tests for functional/search.py — focus on PFAM accession handling."""

from __future__ import annotations

from pathlib import Path

from eukan.functional.search import _add_func_info, _pfam_accession
from eukan.gff import create_gff_db


# ---------------------------------------------------------------------------
# _pfam_accession
# ---------------------------------------------------------------------------


class TestPfamAccession:
    def test_strips_version_suffix(self):
        info = {"description": "...", "evalue": 1e-20, "accession": "PF00710.15"}
        assert _pfam_accession(info, "Asparaginase") == "PF00710"

    def test_no_version_passthrough(self):
        info = {"description": "...", "evalue": 1e-20, "accession": "PF00710"}
        assert _pfam_accession(info, "Asparaginase") == "PF00710"

    def test_falls_back_to_name_when_missing(self):
        # Legacy cache without accession key
        info = {"description": "...", "evalue": 1e-20}
        assert _pfam_accession(info, "Asparaginase") == "Asparaginase"

    def test_falls_back_to_name_when_empty(self):
        info = {"description": "...", "evalue": 1e-20, "accession": ""}
        assert _pfam_accession(info, "Asparaginase") == "Asparaginase"


# ---------------------------------------------------------------------------
# _add_func_info — PFAM inference uses accession, not family name
# ---------------------------------------------------------------------------


class TestAddFuncInfoPfam:
    def _gff(self, tmp_path: Path) -> Path:
        gff = tmp_path / "in.gff3"
        gff.write_text(
            "##gff-version 3\n"
            "chr1\teukan\tgene\t1\t300\t.\t+\t.\tID=g1\n"
            "chr1\teukan\tmRNA\t1\t300\t.\t+\t.\tID=m1;Parent=g1\n"
            "chr1\teukan\tCDS\t1\t300\t.\t+\t0\tID=c1;Parent=m1\n"
        )
        return gff

    def test_pfam_inference_uses_accession(self, tmp_path: Path):
        db = create_gff_db(self._gff(tmp_path))
        phmmer_res = {}
        hmmscan_res = {
            "m1": {
                "Asparaginase": {
                    "description": "Asparaginase domain",
                    "evalue": 1e-30,
                    "accession": "PF00710.15",
                }
            }
        }
        feats = list(_add_func_info(db, phmmer_res, hmmscan_res))
        # mRNA carries the inference attribute
        mrna = next(f for f in feats if f.featuretype == "mRNA")
        assert mrna.attributes["inference"] == ["protein motif:PFAM:PF00710"]

    def test_legacy_cache_falls_back_to_name(self, tmp_path: Path):
        """Old hmmscan.json without accession field still works."""
        db = create_gff_db(self._gff(tmp_path))
        phmmer_res = {}
        hmmscan_res = {
            "m1": {
                "Asparaginase": {
                    "description": "Asparaginase domain",
                    "evalue": 1e-30,
                }
            }
        }
        feats = list(_add_func_info(db, phmmer_res, hmmscan_res))
        mrna = next(f for f in feats if f.featuretype == "mRNA")
        # Falls back to name when accession is missing
        assert mrna.attributes["inference"] == ["protein motif:PFAM:Asparaginase"]

    def test_combines_phmmer_and_pfam_inferences(self, tmp_path: Path):
        db = create_gff_db(self._gff(tmp_path))
        phmmer_res = {
            "m1": {
                "sp|Q9UTS7|ASPG2_SCHPO": {
                    "description": "Probable L-asparaginase 2",
                    "evalue": 1e-50,
                }
            }
        }
        hmmscan_res = {
            "m1": {
                "Asparaginase": {
                    "description": "Asparaginase domain",
                    "evalue": 1e-30,
                    "accession": "PF00710.15",
                }
            }
        }
        feats = list(_add_func_info(db, phmmer_res, hmmscan_res))
        mrna = next(f for f in feats if f.featuretype == "mRNA")
        assert mrna.attributes["inference"] == [
            "similar to AA sequence:UniProtKB:sp|Q9UTS7|ASPG2_SCHPO",
            "protein motif:PFAM:PF00710",
        ]
