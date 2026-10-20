"""Tests for functional/search.py — focus on PFAM accession handling."""

from __future__ import annotations

import itertools
from pathlib import Path

import pyhmmer
import pytest

from eukan.functional.search import (
    _add_func_info,
    _open_hmm_db,
    _pfam_accession,
    annotate_gff3,
    apply_uniprot_hits,
    run_hmmscan_search,
)
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
        feats = list(_add_func_info(db, phmmer_res, hmmscan_res, apply_uniprot_hits))
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
        feats = list(_add_func_info(db, phmmer_res, hmmscan_res, apply_uniprot_hits))
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
        feats = list(_add_func_info(db, phmmer_res, hmmscan_res, apply_uniprot_hits))
        mrna = next(f for f in feats if f.featuretype == "mRNA")
        assert mrna.attributes["inference"] == [
            "similar to AA sequence:UniProtKB:sp|Q9UTS7|ASPG2_SCHPO",
            "protein motif:PFAM:PF00710",
        ]


# ---------------------------------------------------------------------------
# annotate_gff3 — output routing
# ---------------------------------------------------------------------------


class TestAnnotateGff3OutputRouting:
    """The .mod.gff3 must land in the func-annot work_dir, not next to the
    input. ``Artifact.FINAL_FUNC_GFF3`` declares ``func-annot`` as the
    producer; auto-discovery in prep-submission and the test runner both
    look there. Writing next to the input (e.g. ``annotate/final.gff3``)
    leaves the file in the wrong place and downstream steps skip it.
    """

    def _write_minimal_gff(self, path: Path) -> None:
        path.write_text(
            "##gff-version 3\n"
            "chr1\ttest\tgene\t100\t300\t.\t+\t.\tID=g1\n"
            "chr1\ttest\tmRNA\t100\t300\t.\t+\t.\tID=m1;Parent=g1\n"
            "chr1\ttest\tCDS\t100\t300\t.\t+\t0\tID=c1;Parent=m1\n"
        )

    def test_output_dir_overrides_input_parent(self, tmp_path):
        annotate_dir = tmp_path / "annotate"
        func_dir = tmp_path / "func-annot"
        annotate_dir.mkdir()
        func_dir.mkdir()
        gff3 = annotate_dir / "final.gff3"
        self._write_minimal_gff(gff3)

        out = annotate_gff3(gff3, {}, {}, apply_uniprot_hits, output_dir=func_dir)

        assert out == func_dir / "final.mod.gff3"
        assert out.exists()
        assert not (annotate_dir / "final.mod.gff3").exists()

    def test_default_writes_next_to_input(self, tmp_path):
        gff3 = tmp_path / "final.gff3"
        self._write_minimal_gff(gff3)

        out = annotate_gff3(gff3, {}, {}, apply_uniprot_hits)

        assert out == tmp_path / "final.mod.gff3"
        assert out.exists()


# ---------------------------------------------------------------------------
# run_hmmscan_search — non-overlapping domain hits (real pyhmmer run)
# ---------------------------------------------------------------------------

# The PyPI pyhmmer wheel ships its own test data, including these RRE family HMMs.
RREFAM = Path(pyhmmer.__file__).parent / "tests" / "data" / "hmms" / "txt" / "RREFam.hmm"


@pytest.fixture(scope="module")
def rrefam() -> list[pyhmmer.plan7.HMM]:
    if not RREFAM.exists():
        pytest.skip("pyhmmer installed without its bundled test HMMs")
    with pyhmmer.plan7.HMMFile(str(RREFAM)) as hf:
        return list(hf)


class TestHmmscanNonOverlap:
    """A protein made of two RRE domains (the HMMs' consensus sequences)
    joined by a linker. Thiopeptide_F_RRE has the best E-value in either
    order, and both families also match the other domain's region weakly,
    so a hit's first domain is not always its best one."""

    @pytest.mark.parametrize("lasso_first", [True, False], ids=["lasso_first", "lasso_last"])
    def test_both_domains_reported_in_either_order(self, rrefam, lasso_first):
        by_name = {hmm.name: hmm for hmm in rrefam}
        lasso = by_name["Stand_Alone_Lasso_RRE"].consensus.upper()
        thiopeptide = by_name["Thiopeptide_F_RRE"].consensus.upper()
        linker = "GGSGGS" * 4
        domains = (lasso, thiopeptide) if lasso_first else (thiopeptide, lasso)
        query = pyhmmer.easel.TextSequence(name="q", sequence=linker.join(domains)).digitize(
            pyhmmer.easel.Alphabet.amino(),
        )

        res = run_hmmscan_search([query], rrefam, num_cpu=1, evalue_threshold=1e-3)

        assert set(res["q"]) == {"Stand_Alone_Lasso_RRE", "Thiopeptide_F_RRE"}


# The same RRE family HMMs, pressed (.h3m/.h3f/.h3p/.h3i beside the .hmm).
PRESSED_RREFAM = RREFAM.parents[1] / "db" / "RREFam.hmm"


class TestPressedDatabase:
    def test_pressed_hits_match_loaded_hits(self, rrefam):
        """A pressed database's optimized profiles are read straight into one
        block; the hits are those of the HMMs loaded whole, and of the same
        models in the unpressed file."""
        linker = "GGSGGS" * 4
        consensus = [hmm.consensus.upper() for hmm in rrefam[:4]]
        queries = [
            pyhmmer.easel.TextSequence(name=f"q{i}", sequence=a + linker + b).digitize(
                pyhmmer.easel.Alphabet.amino(),
            )
            for i, (a, b) in enumerate(itertools.permutations(consensus, 2))
        ]

        with _open_hmm_db(PRESSED_RREFAM) as profiles:
            assert isinstance(profiles, pyhmmer.plan7.OptimizedProfileBlock)
            assert len(profiles) == len(rrefam)
            pressed = run_hmmscan_search(queries, profiles, num_cpu=2, evalue_threshold=1e-3)
        with _open_hmm_db(RREFAM) as profiles:
            unpressed = run_hmmscan_search(queries, profiles, num_cpu=2, evalue_threshold=1e-3)
        loaded = run_hmmscan_search(queries, rrefam, num_cpu=2, evalue_threshold=1e-3)

        assert len(pressed) == len(queries)
        assert pressed == unpressed == loaded
