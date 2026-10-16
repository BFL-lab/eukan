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
import shutil

import pytest

from eukan.functional import pipeline as fpipe
from eukan.infra.manifest import (
    MANIFEST_FILE,
    RunManifest,
    StepRecord,
    StepStatus,
    load_manifest,
    save_manifest,
)
from eukan.infra.steps import fingerprint_inputs
from eukan.infra.utils import md5_file
from eukan.settings import FunctionalConfig


def _config(proteins, work_dir, *, homology_db="uniprot", gff3_path=None):
    """Build a FunctionalConfig. The caller writes *proteins* first."""
    return FunctionalConfig(
        proteins=proteins, work_dir=work_dir, manifest_dir=work_dir,
        homology_db=homology_db, gff3_path=gff3_path, num_cpu=1,
    )


# --- stubs --------------------------------------------------------------------

def _fake_search(calls):
    """A _search_and_cache stand-in recording proteome stems; writes the cache."""
    def _search_and_cache(config, cache):
        calls.append(config.proteins.stem)
        results = {"homology": {}, "pfam": {}}
        cache.write_text(json.dumps(results))
        return results
    return _search_and_cache


def _fake_annotate_fasta(proteins, homology_res, hmmscan_res, format_hit):
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

    The stubs write the same output files the real steps would (JSON cache,
    ``.mod.faa``, ``.mod.gff3``) so the manifest records valid, integrity-
    checkable outputs.
    """
    counts = {"search": 0, "annotate_fasta": 0, "annotate_gff3": 0}

    def fake_search(config, cache):
        counts["search"] += 1
        results = {"homology": {}, "pfam": {}}
        cache.write_text(json.dumps(results))
        return results

    def fake_annotate_fasta(proteins, homology_res, hmmscan_res, format_hit):
        counts["annotate_fasta"] += 1
        out = proteins.parent / f"{proteins.stem}.mod{proteins.suffix}"
        out.write_text(">s1 hypothetical protein\nMACE\n")
        return out

    def fake_annotate_gff3(gff3_path, homology_res, hmmscan_res, apply_hits, output_dir=None):
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
        assert (tmp_path / "speciesA.phmmer.search.json").exists()
        assert (tmp_path / "speciesB.phmmer.search.json").exists()

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
        assert (tmp_path / "speciesA.phmmer.search.json").exists()
        assert (tmp_path / "speciesA.kofam.search.json").exists()


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
        cache = tmp_path / "p.phmmer.search.json"

        fp_map = fpipe._func_step_fingerprints(cfg, cache)

        search_files, search_scalars = fp_map["search"]
        assert search_files == [proteins]
        assert "homology_db=uniprot" in search_scalars
        assert any(s.startswith("evalue=") for s in search_scalars)
        assert any(s.startswith("uniprot_db=") for s in search_scalars)
        # thread count must never enter the fingerprint
        assert not any("num_cpu" in s for s in search_scalars)

        fasta_files, _ = fp_map["annotate_fasta"]
        assert fasta_files == [proteins, cache]  # depends on the search cache
        assert "annotate_gff3" not in fp_map  # no gff3_path

    def test_gff3_entry_present_when_set(self, tmp_path):
        proteins = tmp_path / "p.faa"
        proteins.write_text(">s\nM\n")
        gff3 = tmp_path / "genes.gff3"
        gff3.write_text("##gff-version 3\n")
        cfg = _config(proteins, tmp_path, gff3_path=gff3)
        cache = tmp_path / "p.phmmer.search.json"

        gff3_files, _ = fpipe._func_step_fingerprints(cfg, cache)["annotate_gff3"]
        assert gff3_files == [gff3, cache]

    def test_kofam_mode_scalars(self, tmp_path):
        proteins = tmp_path / "p.faa"
        proteins.write_text(">s\nM\n")
        cfg = _config(proteins, tmp_path, homology_db="kofam")
        cache = tmp_path / "p.kofam.search.json"

        _, scalars = fpipe._func_step_fingerprints(cfg, cache)["search"]
        assert any(s.startswith("kofam_db=") for s in scalars)
        assert any(s.startswith("ko_list=") for s in scalars)
        assert not any(s.startswith("uniprot_db=") for s in scalars)


# --- the search cache ---------------------------------------------------------


_HOMOLOGY = {"s1": {"P12345": {"description": "Kinase X", "evalue": 1e-30}}}
_PFAM = {"s1": {"Pkinase": {
    "description": "Protein kinase domain", "evalue": 1e-20, "accession": "PF00069.1",
}}}


class TestSearchCache:
    """The search caches its homology and Pfam hits together in one file, its
    declared output, so a resume never reads a search file it does not track and
    a lost cache re-runs the search. Runs the real _search_and_cache with the
    database searches stubbed."""

    @staticmethod
    def _stub_searches(monkeypatch) -> list[str]:
        searches: list[str] = []
        monkeypatch.setattr(
            fpipe, "_run_uniprot_phmmer",
            lambda *a: searches.append("phmmer") or _HOMOLOGY,
        )
        monkeypatch.setattr(
            fpipe, "_run_pfam_hmmscan", lambda *a: searches.append("hmmscan") or _PFAM,
        )
        return searches

    @staticmethod
    def _search_record(manifest):
        (record,) = [r for k, r in manifest.steps.items() if k.startswith("functional/search")]
        return record

    def test_resume_reads_no_untracked_search_file(self, tmp_path, monkeypatch):
        # The Pfam hits used to sit in an untracked <stem>.hmmscan.json beside the
        # tracked homology cache: with it removed the search still counted as
        # complete and the resume crashed reading it.
        self._stub_searches(monkeypatch)
        prot = tmp_path / "speciesA.faa"
        prot.write_text(">s1\nMACE\n")
        work_dir = tmp_path / "func-annot"
        cfg = _config(prot, work_dir)
        fpipe.run_functional_annotation(cfg)
        tracked = work_dir / self._search_record(load_manifest(work_dir)).output_file
        for path in tmp_path.glob("speciesA.*.json"):
            if path != tracked:
                path.unlink()
        fpipe.run_functional_annotation(cfg)
        assert "Protein kinase domain" in (tmp_path / "speciesA.mod.faa").read_text()

    def test_lost_cache_reruns_search(self, tmp_path, monkeypatch):
        searches = self._stub_searches(monkeypatch)
        prot = tmp_path / "speciesA.faa"
        prot.write_text(">s1\nMACE\n")
        cfg = _config(prot, tmp_path / "func-annot")
        fpipe.run_functional_annotation(cfg)
        (tmp_path / "speciesA.phmmer.search.json").unlink()
        fpipe.run_functional_annotation(cfg)
        assert searches == ["phmmer", "hmmscan"] * 2
        assert "Protein kinase domain" in (tmp_path / "speciesA.mod.faa").read_text()

    def test_record_from_the_two_file_layout_reruns_search(self, tmp_path, monkeypatch):
        # A run dir from before the single cache: the search record points at the
        # homology-only <stem>.phmmer.json under a fingerprint without the cache
        # name. Trusting it would skip the search and then miss the cache.
        searches = self._stub_searches(monkeypatch)
        prot = tmp_path / "speciesA.faa"
        prot.write_text(">s1\nMACE\n")
        work_dir = tmp_path / "func-annot"
        cfg = _config(prot, work_dir)
        fpipe.run_functional_annotation(cfg)
        cache = fpipe._search_cache_path(cfg)
        cache.unlink()
        old = tmp_path / "speciesA.phmmer.json"
        old.write_text(json.dumps(_HOMOLOGY))
        files, scalars = fpipe._func_step_fingerprints(cfg, cache)["search"]
        manifest = load_manifest(work_dir)
        record = self._search_record(manifest)
        record.output_file, record.output_md5 = str(old), md5_file(old)
        record.input_md5 = fingerprint_inputs(
            files, [s for s in scalars if not s.startswith("cache=")],
        )
        save_manifest(work_dir, manifest)

        fpipe.run_functional_annotation(cfg)

        assert searches == ["phmmer", "hmmscan"] * 2
        assert cache.exists()


class TestMovedRunDir:
    def test_moved_run_dir_resumes(self, tmp_path, monkeypatch):
        """A proteome inside a run dir moved elsewhere keeps its step keys, cache
        and fingerprints: the resume is a no-op."""
        counts = _install_counts(monkeypatch)

        def config(run_dir):
            return FunctionalConfig(
                proteins=run_dir / "speciesA.faa", work_dir=run_dir / "func-annot",
                manifest_dir=run_dir, num_cpu=1,
            )

        run_a, run_b = tmp_path / "runA", tmp_path / "runB"
        run_a.mkdir()
        (run_a / "speciesA.faa").write_text(">s1\nMACE\n")
        fpipe.run_functional_annotation(config(run_a))
        shutil.move(run_a, run_b)

        fpipe.run_functional_annotation(config(run_b))

        assert counts == {"search": 1, "annotate_fasta": 1, "annotate_gff3": 0}


# --- manifest location --------------------------------------------------------


class TestManifestLocation:
    """The functional steps record into the run-dir manifest (``manifest_dir``),
    shared with the other pipelines; their outputs stay under ``work_dir``."""

    def test_records_into_shared_run_dir_manifest(self, tmp_path, monkeypatch):
        counts = _install_counts(monkeypatch)
        run_dir = tmp_path / "run"
        work_dir = run_dir / "func-annot"
        run_dir.mkdir()
        # An earlier pipeline's record already in the run-dir manifest.
        save_manifest(run_dir, RunManifest(steps={
            "annotation/final": StepRecord(
                name="annotation/final", status=StepStatus.completed,
            ),
        }))
        prot = tmp_path / "speciesA.faa"
        prot.write_text(">s1\nMACE\n")
        gff3 = tmp_path / "genes.gff3"
        gff3.write_text("##gff-version 3\nchr1\teukan\tgene\t1\t100\t.\t+\t.\tID=g1\n")
        cfg = FunctionalConfig(
            proteins=prot, work_dir=work_dir, manifest_dir=run_dir,
            gff3_path=gff3, num_cpu=1,
        )

        fpipe.run_functional_annotation(cfg)

        assert not (work_dir / MANIFEST_FILE).exists()
        manifest = load_manifest(run_dir)
        assert manifest is not None
        assert "annotation/final" in manifest.steps  # added to, not replaced
        assert {k.split(".")[0] for k in manifest.steps if k.startswith("functional/")} == {
            "functional/search", "functional/annotate_fasta", "functional/annotate_gff3",
        }
        assert (work_dir / "genes.mod.gff3").exists()
        # A resume reads the same manifest back: an identical re-run is a no-op.
        fpipe.run_functional_annotation(cfg)
        assert counts == {"search": 1, "annotate_fasta": 1, "annotate_gff3": 1}
