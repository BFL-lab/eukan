"""Tests for the functional annotation pipeline driver.

Two layered resume concerns:

* Manifest step keys are scoped per proteome+mode (``_step_scope``) so several
  proteomes can be annotated in one work_dir without the second inheriting the
  first's "search complete" record (which used to crash reading a cache file
  that was never written). See TestStepScope / TestMultiProteomeResume.
* Within a scope, each step declares input fingerprints
  (``_func_step_fingerprints``) so the same proteome edited in place re-runs and
  a byte-identical re-run is a friendly no-op pointing at ``-f``. See
  TestFuncAnnotInputIdentity / TestFuncStepFingerprints.

The heavy homology search + annotation are monkeypatched: the scoping/resume
logic under test lives in the driver, and the real steps need multi-GB DBs.
"""

from __future__ import annotations

import json
import logging

import pytest

from eukan.functional import pipeline as fpipe
from eukan.settings import FunctionalConfig


def _config(proteins, work_dir, *, homology_db="uniprot", gff3_path=None):
    """Build a FunctionalConfig. The caller writes *proteins* first."""
    return FunctionalConfig(
        proteins=proteins, work_dir=work_dir, manifest_dir=work_dir,
        homology_db=homology_db, gff3_path=gff3_path, num_cpu=1,
    )


# --- stubs --------------------------------------------------------------------

def _fake_search(calls):
    """A _search_and_cache stand-in recording proteome stems; writes caches."""
    def _search_and_cache(config, homology_json, hmmscan_json):
        calls.append(config.proteins.stem)
        homology_json.write_text(json.dumps({}))
        hmmscan_json.write_text(json.dumps({}))
        return homology_json
    return _search_and_cache


def _fake_annotate_fasta(proteins, homology_res, hmmscan_res, homology_db="uniprot"):
    out = proteins.parent / f"{proteins.stem}.mod{proteins.suffix}"
    out.write_text(">x\nMA\n")
    return out


@pytest.fixture
def calls(monkeypatch):
    """Record the proteome stem each search runs for (search + fasta stubbed)."""
    recorded: list[str] = []
    monkeypatch.setattr(fpipe, "_search_and_cache", _fake_search(recorded))
    monkeypatch.setattr(fpipe, "annotate_fasta", _fake_annotate_fasta)
    return recorded


def _install_counts(monkeypatch) -> dict[str, int]:
    """Stub all three steps with per-step invocation counters.

    The stubs write the same output files the real steps would (JSON caches,
    ``.mod.faa``, ``.mod.gff3``) so the manifest records valid, integrity-
    checkable outputs.
    """
    counts = {"search": 0, "annotate_fasta": 0, "annotate_gff3": 0}

    def fake_search(config, homology_json, hmmscan_json):
        counts["search"] += 1
        homology_json.write_text("{}")
        hmmscan_json.write_text("{}")
        return homology_json

    def fake_annotate_fasta(proteins, homology_res, hmmscan_res, homology_db="uniprot"):
        counts["annotate_fasta"] += 1
        out = proteins.parent / f"{proteins.stem}.mod{proteins.suffix}"
        out.write_text(">s1 hypothetical protein\nMACE\n")
        return out

    def fake_annotate_gff3(
        gff3_path, homology_res, hmmscan_res, output_dir=None, homology_db="uniprot"
    ):
        counts["annotate_gff3"] += 1
        out = (output_dir or gff3_path.parent) / f"{gff3_path.stem}.mod{gff3_path.suffix}"
        out.write_text("##gff-version 3\nchr1\teukan\tgene\t1\t100\t.\t+\t.\tID=g1\n")
        return out

    monkeypatch.setattr(fpipe, "_search_and_cache", fake_search)
    monkeypatch.setattr(fpipe, "annotate_fasta", fake_annotate_fasta)
    monkeypatch.setattr(fpipe, "annotate_gff3", fake_annotate_gff3)
    return counts


# --- per-proteome+mode scoping ------------------------------------------------

class TestStepScope:
    def test_distinct_per_proteome_mode_and_dir(self, tmp_path):
        a = tmp_path / "a.faa"
        a.touch()
        b = tmp_path / "b.faa"
        b.touch()
        sub = tmp_path / "sub"
        sub.mkdir()
        a_other_dir = sub / "a.faa"
        a_other_dir.touch()

        scope = fpipe._step_scope
        # Different proteomes -> different scopes.
        assert scope(_config(a, tmp_path)) != scope(_config(b, tmp_path))
        # Same proteome, different homology DB -> different scopes.
        assert scope(_config(a, tmp_path, homology_db="uniprot")) != scope(
            _config(a, tmp_path, homology_db="kofam")
        )
        # Same stem in a different directory -> different scopes (distinct caches).
        assert scope(_config(a, tmp_path)) != scope(_config(a_other_dir, tmp_path))
        # Deterministic for the same inputs.
        assert scope(_config(a, tmp_path)) == scope(_config(a, tmp_path))


class TestMultiProteomeResume:
    def test_second_proteome_same_workdir_does_not_crash(self, tmp_path, calls):
        work_dir = tmp_path / "func-annot"
        prot_a = tmp_path / "speciesA.faa"
        prot_a.write_text(">x\nMA\n")
        fpipe.run_functional_annotation(_config(prot_a, work_dir))

        prot_b = tmp_path / "speciesB.faa"
        prot_b.write_text(">y\nMB\n")
        # Before the scoping fix this raised FileNotFoundError: the search was
        # skipped (manifest said "complete") so speciesB's cache was never
        # written, yet the driver then tried to read it.
        fpipe.run_functional_annotation(_config(prot_b, work_dir))

        assert calls == ["speciesA", "speciesB"]  # each ran its own search
        assert (tmp_path / "speciesA.phmmer.json").exists()
        assert (tmp_path / "speciesB.phmmer.json").exists()

    def test_same_proteome_rerun_is_skipped(self, tmp_path, calls):
        work_dir = tmp_path / "func-annot"
        prot = tmp_path / "speciesA.faa"
        prot.write_text(">x\nMA\n")
        fpipe.run_functional_annotation(_config(prot, work_dir))
        fpipe.run_functional_annotation(_config(prot, work_dir))
        assert calls == ["speciesA"]  # resume: search ran once, re-run skipped

    def test_mode_switch_reruns_search(self, tmp_path, calls):
        work_dir = tmp_path / "func-annot"
        prot = tmp_path / "speciesA.faa"
        prot.write_text(">x\nMA\n")
        fpipe.run_functional_annotation(_config(prot, work_dir, homology_db="uniprot"))
        fpipe.run_functional_annotation(_config(prot, work_dir, homology_db="kofam"))
        assert calls == ["speciesA", "speciesA"]  # each mode writes its own cache
        assert (tmp_path / "speciesA.phmmer.json").exists()
        assert (tmp_path / "speciesA.kofam.json").exists()


# --- content-sensitivity within a scope + friendly no-op ----------------------

class TestFuncAnnotInputIdentity:
    def test_different_input_reruns_without_force(self, tmp_path, monkeypatch):
        counts = _install_counts(monkeypatch)
        a = tmp_path / "prot_a.faa"
        a.write_text(">a\nMAAA\n")
        b = tmp_path / "prot_b.faa"
        b.write_text(">b\nMBBB\n")

        fpipe.run_functional_annotation(_config(a, tmp_path))
        assert (counts["search"], counts["annotate_fasta"]) == (1, 1)
        fpipe.run_functional_annotation(_config(b, tmp_path))  # different proteome
        assert (counts["search"], counts["annotate_fasta"]) == (2, 2)

    def test_same_name_different_content_reruns(self, tmp_path, monkeypatch):
        counts = _install_counts(monkeypatch)
        prot = tmp_path / "proteins.faa"
        prot.write_text(">a\nMAAA\n")
        fpipe.run_functional_annotation(_config(prot, tmp_path))
        prot.write_text(">a\nMCCC\n")  # same path+scope, edited in place
        fpipe.run_functional_annotation(_config(prot, tmp_path))
        # Scope alone wouldn't catch this (same path); the content fingerprint does.
        assert counts["search"] == 2

    def test_identical_input_is_noop(self, tmp_path, monkeypatch, caplog):
        counts = _install_counts(monkeypatch)
        prot = tmp_path / "proteins.faa"
        prot.write_text(">s1\nMACE\n")

        fpipe.run_functional_annotation(_config(prot, tmp_path))
        with caplog.at_level(logging.INFO, logger="eukan.functional.pipeline"):
            fpipe.run_functional_annotation(_config(prot, tmp_path))  # identical -> no-op

        assert counts["search"] == 1
        assert counts["annotate_fasta"] == 1
        assert "Already annotated" in caplog.text
        assert "-f" in caplog.text

    def test_identical_input_with_force_reruns(self, tmp_path, monkeypatch):
        counts = _install_counts(monkeypatch)
        prot = tmp_path / "proteins.faa"
        prot.write_text(">s1\nMACE\n")
        fpipe.run_functional_annotation(_config(prot, tmp_path))
        fpipe.run_functional_annotation(_config(prot, tmp_path), force=True)
        assert counts["search"] == 2
        assert counts["annotate_fasta"] == 2

    def test_added_gff3_runs_without_force(self, tmp_path, monkeypatch):
        """Re-running the same proteome but newly asking for GFF3 annotation
        runs only the new step, no -f needed."""
        counts = _install_counts(monkeypatch)
        prot = tmp_path / "proteins.faa"
        prot.write_text(">s1\nMACE\n")
        fpipe.run_functional_annotation(_config(prot, tmp_path))
        assert counts["annotate_gff3"] == 0

        gff3 = tmp_path / "genes.gff3"
        gff3.write_text("##gff-version 3\nchr1\teukan\tgene\t1\t100\t.\t+\t.\tID=g1\n")
        fpipe.run_functional_annotation(_config(prot, tmp_path, gff3_path=gff3))
        # search + fasta cached (identical proteome), gff3 newly run.
        assert counts["search"] == 1
        assert counts["annotate_fasta"] == 1
        assert counts["annotate_gff3"] == 1


class TestFuncStepFingerprints:
    """``_func_step_fingerprints`` declares the right inputs/scalars per step."""

    def test_uniprot_scalars_and_inputs(self, tmp_path):
        proteins = tmp_path / "p.faa"
        proteins.write_text(">s\nM\n")
        cfg = _config(proteins, tmp_path)
        hj, kj = tmp_path / "p.phmmer.json", tmp_path / "p.hmmscan.json"

        fp_map = fpipe._func_step_fingerprints(cfg, hj, kj)

        search_files, search_scalars = fp_map["search"]
        assert search_files == [proteins]
        assert "homology_db=uniprot" in search_scalars
        assert any(s.startswith("evalue=") for s in search_scalars)
        assert any(s.startswith("uniprot_db=") for s in search_scalars)
        # thread count must never enter the fingerprint
        assert not any("num_cpu" in s for s in search_scalars)

        fasta_files, _ = fp_map["annotate_fasta"]
        assert fasta_files == [proteins, hj, kj]  # depends on both JSON caches
        assert "annotate_gff3" not in fp_map  # no gff3_path

    def test_gff3_entry_present_when_set(self, tmp_path):
        proteins = tmp_path / "p.faa"
        proteins.write_text(">s\nM\n")
        gff3 = tmp_path / "genes.gff3"
        gff3.write_text("##gff-version 3\n")
        cfg = _config(proteins, tmp_path, gff3_path=gff3)
        hj, kj = tmp_path / "p.phmmer.json", tmp_path / "p.hmmscan.json"

        gff3_files, _ = fpipe._func_step_fingerprints(cfg, hj, kj)["annotate_gff3"]
        assert gff3_files == [gff3, hj, kj]

    def test_kofam_mode_scalars(self, tmp_path):
        proteins = tmp_path / "p.faa"
        proteins.write_text(">s\nM\n")
        cfg = _config(proteins, tmp_path, homology_db="kofam")
        hj, kj = tmp_path / "p.kofam.json", tmp_path / "p.hmmscan.json"

        _, scalars = fpipe._func_step_fingerprints(cfg, hj, kj)["search"]
        assert any(s.startswith("kofam_db=") for s in scalars)
        assert any(s.startswith("ko_list=") for s in scalars)
        assert not any(s.startswith("uniprot_db=") for s in scalars)
