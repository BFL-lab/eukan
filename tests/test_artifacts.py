"""Tests for eukan.infra.artifacts — find_or_warn warn-on-loss semantics."""

from __future__ import annotations

import logging

from eukan.infra.artifacts import Artifact, find_or_warn
from eukan.infra.manifest import RunManifest, StepRecord, StepStatus


def _manifest_with(step_key: str, status: StepStatus = StepStatus.completed) -> RunManifest:
    manifest = RunManifest()
    manifest.steps[step_key] = StepRecord(name=step_key, status=status)
    return manifest


class TestFindOrWarn:
    """A missing optional artifact warns only when its producer actually ran."""

    def test_present_returns_path(self, tmp_path):
        f = tmp_path / Artifact.REPEATMASK_HINTS.value
        f.write_text("data")
        assert find_or_warn(tmp_path, Artifact.REPEATMASK_HINTS, None) == f

    def test_missing_producer_completed_warns(self, tmp_path, caplog):
        manifest = _manifest_with("repeats/masker")  # mask-repeats producer ran
        with caplog.at_level(logging.WARNING, logger="eukan.infra.artifacts"):
            result = find_or_warn(tmp_path, Artifact.REPEATMASK_HINTS, manifest)
        assert result is None
        assert "now missing" in caplog.text
        assert "mask-repeats" in caplog.text

    def test_missing_producer_never_ran_is_silent(self, tmp_path, caplog):
        with caplog.at_level(logging.WARNING, logger="eukan.infra.artifacts"):
            result = find_or_warn(tmp_path, Artifact.REPEATMASK_HINTS, RunManifest())
        assert result is None
        assert "now missing" not in caplog.text

    def test_missing_no_manifest_is_silent(self, tmp_path, caplog):
        with caplog.at_level(logging.WARNING, logger="eukan.infra.artifacts"):
            result = find_or_warn(tmp_path, Artifact.SPLICE_SUMMARY, None)
        assert result is None
        assert "now missing" not in caplog.text

    def test_producer_failed_not_completed_is_silent(self, tmp_path, caplog):
        # Producer step is in the manifest but failed -> not a genuine loss.
        manifest = _manifest_with("repeats/masker", status=StepStatus.failed)
        with caplog.at_level(logging.WARNING, logger="eukan.infra.artifacts"):
            result = find_or_warn(tmp_path, Artifact.REPEATMASK_HINTS, manifest)
        assert result is None
        assert "now missing" not in caplog.text
