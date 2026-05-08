"""Tests for infra/layout.py and the artifacts cross-step fallback."""

from __future__ import annotations

from pathlib import Path

import pytest

from eukan.infra import artifacts
from eukan.infra.artifacts import Artifact, masked_genome
from eukan.infra.layout import (
    PIPELINE_SUBDIRS,
    sibling_step_dir,
    step_work_dir,
)


# ---------------------------------------------------------------------------
# Layout helpers
# ---------------------------------------------------------------------------


class TestStepWorkDir:
    def test_known_steps_resolve(self, tmp_path: Path):
        for step in ("annotate", "assemble", "mask-repeats", "func-annot", "prep-submission"):
            wd = step_work_dir(step, tmp_path)
            assert wd.parent == tmp_path
            assert wd.name == PIPELINE_SUBDIRS[step]

    def test_unknown_step_raises(self, tmp_path: Path):
        with pytest.raises(KeyError):
            step_work_dir("nonexistent", tmp_path)

    def test_default_base_is_cwd(self, tmp_path: Path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        assert step_work_dir("annotate") == tmp_path / "annotate"


class TestSiblingStepDir:
    def test_resolves_sibling_under_same_root(self, tmp_path: Path):
        annotate_wd = step_work_dir("annotate", tmp_path)
        assert sibling_step_dir(annotate_wd, "assemble") == tmp_path / "assemble"
        assert sibling_step_dir(annotate_wd, "func-annot") == tmp_path / "func-annot"


# ---------------------------------------------------------------------------
# Artifact resolution — flat layout vs. step layout
# ---------------------------------------------------------------------------


class TestFindArtifact:
    def test_flat_layout_finds_in_work_dir(self, tmp_path: Path):
        """Legacy: artifact written directly into the caller's work_dir."""
        (tmp_path / Artifact.FINAL_GFF3.value).write_text("##gff-version 3\n")
        assert artifacts.find(tmp_path, Artifact.FINAL_GFF3) is not None

    def test_step_layout_falls_back_to_producer(self, tmp_path: Path):
        """Caller's work_dir is empty; artifact is in the producer's sibling dir."""
        annotate_dir = step_work_dir("annotate", tmp_path)
        annotate_dir.mkdir()
        (annotate_dir / Artifact.FINAL_GFF3.value).write_text("##gff-version 3\n")

        # Submission, looking from its own step dir, should still find it.
        submission_dir = step_work_dir("prep-submission", tmp_path)
        submission_dir.mkdir()
        found = artifacts.find(submission_dir, Artifact.FINAL_GFF3)
        assert found == annotate_dir / Artifact.FINAL_GFF3.value

    def test_own_dir_wins_over_sibling(self, tmp_path: Path):
        """If both exist, the caller's own work_dir takes precedence."""
        annotate_dir = step_work_dir("annotate", tmp_path)
        annotate_dir.mkdir()
        (annotate_dir / Artifact.FINAL_GFF3.value).write_text("from-annotate")

        # final.gff3 also present in submission/ for some reason
        submission_dir = step_work_dir("prep-submission", tmp_path)
        submission_dir.mkdir()
        (submission_dir / Artifact.FINAL_GFF3.value).write_text("from-submission")

        found = artifacts.find(submission_dir, Artifact.FINAL_GFF3)
        assert found.read_text() == "from-submission"

    def test_returns_none_when_absent(self, tmp_path: Path):
        annotate_dir = step_work_dir("annotate", tmp_path)
        annotate_dir.mkdir()
        # nothing written
        submission_dir = step_work_dir("prep-submission", tmp_path)
        submission_dir.mkdir()
        assert artifacts.find(submission_dir, Artifact.FINAL_GFF3) is None


class TestMaskedGenome:
    def test_falls_back_to_repeats_dir(self, tmp_path: Path):
        repeats_dir = step_work_dir("mask-repeats", tmp_path)
        repeats_dir.mkdir()
        (repeats_dir / "genome.masked.fasta").write_text(">chr\nACGT\n")

        annotate_dir = step_work_dir("annotate", tmp_path)
        annotate_dir.mkdir()
        found = masked_genome(annotate_dir, "genome")
        assert found == repeats_dir / "genome.masked.fasta"
        assert found.exists()

    def test_returns_in_step_path_when_missing(self, tmp_path: Path):
        annotate_dir = step_work_dir("annotate", tmp_path)
        annotate_dir.mkdir()
        # nothing written anywhere
        result = masked_genome(annotate_dir, "genome")
        # Caller can use this path to create the file; just verify it's a Path.
        assert isinstance(result, Path)
        assert result.name == "genome.masked.fasta"
