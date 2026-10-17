"""Tests for the run manifest data model and step lifecycle."""

import json
import subprocess

import pytest

from eukan.cli.status import format_status
from eukan.infra import manifest as manifest_mod
from eukan.infra.manifest import (
    RunManifest,
    StepRecord,
    StepStatus,
    load_manifest,
    save_manifest,
)
from eukan.infra.steps import (
    SENTINEL,
    _should_verify_md5,
    clean_interrupted_step,
    is_step_complete,
    is_step_interrupted,
    pipeline_step,
)
from eukan.infra.utils import md5_file


class TestManifestIO:
    def test_round_trip(self, tmp_path):
        """Save and load should preserve all fields."""
        manifest = RunManifest(
            started_at="2026-01-01T00:00:00+00:00",
            genome="/data/genome.fa",
            kingdom="protist",
            genetic_code="6",
            num_cpu=8,
        )
        manifest.steps["genemark"] = StepRecord(
            name="genemark", status=StepStatus.completed,
            duration_seconds=120.5,
        )
        save_manifest(tmp_path, manifest)
        loaded = load_manifest(tmp_path)

        assert loaded is not None
        assert loaded.genome == "/data/genome.fa"
        assert loaded.kingdom == "protist"
        assert loaded.num_cpu == 8
        assert "genemark" in loaded.steps
        assert loaded.steps["genemark"].duration_seconds == 120.5

    def test_load_missing(self, tmp_path):
        assert load_manifest(tmp_path) is None

    def test_load_corrupt(self, tmp_path):
        (tmp_path / "eukan-run.json").write_text("{bad json")
        assert load_manifest(tmp_path) is None

    def test_pydantic_serialization(self):
        manifest = RunManifest(started_at="2026-01-01", genome="g.fa")
        text = manifest.model_dump_json()
        loaded = RunManifest.model_validate_json(text)
        assert loaded.genome == "g.fa"


class TestSaveManifestMerge:
    """Copies of one run dir's manifest (two eukan invocations, say) each save
    only what they changed, so neither loses the other's records."""

    @staticmethod
    def _two_copies(tmp_path, **fields):
        save_manifest(tmp_path, RunManifest(**fields))
        return load_manifest(tmp_path), load_manifest(tmp_path)

    @staticmethod
    def _completed(name):
        return StepRecord(name=name, status=StepStatus.completed)

    def test_copies_saving_different_steps_keep_both(self, tmp_path):
        a, b = self._two_copies(tmp_path)
        b.steps["repeats/modeler"] = self._completed("repeats/modeler")
        save_manifest(tmp_path, b)
        a.steps["assembly/minimap2"] = self._completed("assembly/minimap2")
        save_manifest(tmp_path, a)

        assert set(load_manifest(tmp_path).steps) == {"repeats/modeler", "assembly/minimap2"}

    def test_unchanged_stale_record_does_not_overwrite_a_newer_one(self, tmp_path):
        a, b = self._two_copies(tmp_path, steps={
            "repeats/modeler": StepRecord(name="repeats/modeler", status=StepStatus.failed),
        })
        b.steps["repeats/modeler"].status = StepStatus.completed  # b re-ran it
        save_manifest(tmp_path, b)
        a.steps["assembly/minimap2"] = self._completed("assembly/minimap2")
        save_manifest(tmp_path, a)  # a still holds the failed record

        on_disk = load_manifest(tmp_path)
        assert on_disk.steps["repeats/modeler"].status == StepStatus.completed
        assert "assembly/minimap2" in on_disk.steps

    def test_run_level_change_is_kept(self, tmp_path):
        a, b = self._two_copies(tmp_path)
        b.status = "completed"
        save_manifest(tmp_path, b)
        a.steps["assembly/minimap2"] = self._completed("assembly/minimap2")
        save_manifest(tmp_path, a)

        assert load_manifest(tmp_path).status == "completed"


class TestPipelineStepContextManager:
    def test_success_creates_and_removes_sentinel(self, tmp_path):
        manifest = RunManifest()

        with pipeline_step(tmp_path, manifest, "genemark") as step:
            # Sentinel should exist during execution
            assert (tmp_path / "genemark" / SENTINEL).exists()
            assert manifest.steps["genemark"].status == StepStatus.running
            step.output_file = str(tmp_path / "genemark" / "out.gff3")
            # Create the fake output
            (tmp_path / "genemark" / "out.gff3").write_text("fake")

        # After success: sentinel removed, status completed
        assert not (tmp_path / "genemark" / SENTINEL).exists()
        assert manifest.steps["genemark"].status == StepStatus.completed
        assert manifest.steps["genemark"].output_md5 is not None
        assert manifest.steps["genemark"].duration_seconds is not None

    def test_failure_records_error(self, tmp_path):
        manifest = RunManifest()

        try:
            with pipeline_step(tmp_path, manifest, "genemark"):
                raise RuntimeError("tool crashed")
        except RuntimeError:
            pass

        assert not (tmp_path / "genemark" / SENTINEL).exists()
        assert manifest.steps["genemark"].status == StepStatus.failed
        assert "tool crashed" in manifest.steps["genemark"].error

    def test_interrupt_keeps_sentinel(self, tmp_path):
        """Ctrl-C mid-step must not leave the step "running" with its sentinel
        gone: the sentinel stays so the next run cleans up and re-runs it."""
        manifest = RunManifest()

        with pytest.raises(KeyboardInterrupt):
            with pipeline_step(tmp_path, manifest, "genemark"):
                (tmp_path / "genemark" / "partial.gtf").write_text("half")
                raise KeyboardInterrupt

        assert is_step_interrupted(tmp_path, "genemark")
        record = load_manifest(tmp_path).steps["genemark"]
        assert record.status == StepStatus.failed
        assert record.error == "interrupted"


class TestIsStepComplete:
    def test_completed_step(self, tmp_path):
        manifest = RunManifest()
        output = tmp_path / "test_output.txt"
        output.write_text("data")
        manifest.steps["test"] = StepRecord(
            name="test", status=StepStatus.completed,
            output_file=str(output),
        )
        assert is_step_complete(tmp_path, manifest, "test") == output

    def test_incomplete_step(self, tmp_path):
        manifest = RunManifest()
        manifest.steps["test"] = StepRecord(name="test", status=StepStatus.running)
        assert is_step_complete(tmp_path, manifest, "test") is None

    def test_missing_step(self, tmp_path):
        manifest = RunManifest()
        assert is_step_complete(tmp_path, manifest, "nonexistent") is None

    def test_missing_output_file(self, tmp_path):
        manifest = RunManifest()
        manifest.steps["test"] = StepRecord(
            name="test", status=StepStatus.completed,
            output_file=str(tmp_path / "deleted.gff3"),
        )
        assert is_step_complete(tmp_path, manifest, "test") is None


class TestIsStepCompleteIntegrity:
    """A completed step whose output went missing/empty/corrupt re-runs (None)."""

    @staticmethod
    def _completed(tmp_path, name="out.txt", content="data", md5=None):
        manifest = RunManifest()
        output = tmp_path / name
        output.write_text(content)
        manifest.steps["test"] = StepRecord(
            name="test", status=StepStatus.completed,
            output_file=str(output), output_md5=md5,
        )
        return manifest, output

    def test_empty_output_rebuilds(self, tmp_path):
        manifest, output = self._completed(tmp_path, content="")
        assert is_step_complete(tmp_path, manifest, "test") is None

    def test_checksum_mismatch_rebuilds(self, tmp_path):
        manifest, output = self._completed(tmp_path, content="v1")
        manifest.steps["test"].output_md5 = md5_file(output)
        output.write_text("v2-corrupted")  # same path, different bytes
        assert is_step_complete(tmp_path, manifest, "test") is None

    def test_checksum_match_is_reused(self, tmp_path):
        manifest, output = self._completed(tmp_path, content="v1")
        manifest.steps["test"].output_md5 = md5_file(output)
        assert is_step_complete(tmp_path, manifest, "test") == output

    def test_verify_md5_false_skips_checksum(self, tmp_path):
        manifest, output = self._completed(tmp_path, content="v1")
        manifest.steps["test"].output_md5 = md5_file(output)
        output.write_text("v2-corrupted")
        # Opt-out: existence/non-empty still pass, checksum is not checked.
        assert is_step_complete(tmp_path, manifest, "test", verify_md5=False) == output

    def test_unparseable_gff_rebuilds(self, tmp_path):
        manifest, _ = self._completed(tmp_path, name="out.gff3", content="not gff")
        assert is_step_complete(tmp_path, manifest, "test") is None

    def test_valid_gff_is_reused(self, tmp_path):
        gff = "##gff-version 3\nchr1\teukan\tgene\t1\t100\t.\t+\t.\tID=g1\n"
        manifest, output = self._completed(tmp_path, name="out.gff3", content=gff)
        assert is_step_complete(tmp_path, manifest, "test") == output

    def test_md5_skipped_when_policy_declines(self, tmp_path, monkeypatch):
        # When _should_verify_md5 declines (e.g. a large BAM over the size
        # ceiling), a recorded md5 that no longer matches does NOT force a
        # rebuild — existence + non-empty suffice (cost policy).
        import eukan.infra.steps as steps_mod
        out = tmp_path / "aln.bam"
        out.write_bytes(b"data")
        manifest = RunManifest()
        manifest.steps["test"] = StepRecord(
            name="test", status=StepStatus.completed,
            output_file=str(out), output_md5="not-the-real-md5",
        )
        monkeypatch.setattr(steps_mod, "_should_verify_md5", lambda p: False)
        assert is_step_complete(tmp_path, manifest, "test") == out


class TestShouldVerifyMd5:
    def test_text_and_gff_suffixes_verify(self, tmp_path):
        for name in ("a.gff3", "a.fasta", "a.json", "a.txt"):
            assert _should_verify_md5(tmp_path / name)

    def test_small_unknown_suffix_verifies(self, tmp_path):
        p = tmp_path / "a.weird"
        p.write_bytes(b"small")
        assert _should_verify_md5(p)

    def test_large_binary_skips(self, tmp_path, monkeypatch):
        p = tmp_path / "a.bam"
        p.write_bytes(b"x")

        class _BigStat:
            st_size = 300 * 1024 * 1024  # over the 256 MiB ceiling

        monkeypatch.setattr(type(p), "stat", lambda self: _BigStat())
        assert not _should_verify_md5(p)


class TestInterruptionDetection:
    def test_detects_sentinel(self, tmp_path):
        sdir = tmp_path / "genemark"
        sdir.mkdir()
        (sdir / SENTINEL).write_text("started")
        assert is_step_interrupted(tmp_path, "genemark")

    def test_no_false_positive(self, tmp_path):
        sdir = tmp_path / "genemark"
        sdir.mkdir()
        assert not is_step_interrupted(tmp_path, "genemark")

    def test_clean_removes_dir(self, tmp_path):
        sdir = tmp_path / "genemark"
        sdir.mkdir()
        (sdir / SENTINEL).write_text("started")
        (sdir / "partial_output.gff3").write_text("partial")

        clean_interrupted_step(tmp_path, "genemark")
        assert not sdir.exists()


class TestFormatStatus:
    def test_basic_format(self):
        manifest = RunManifest(
            started_at="2026-01-01T00:00:00+00:00",
            genome="/data/genome.fa",
            kingdom="protist",
            num_cpu=8,
        )
        manifest.steps["genemark"] = StepRecord(
            name="genemark", status=StepStatus.completed,
            duration_seconds=120.5,
        )
        manifest.steps["augustus"] = StepRecord(
            name="augustus", status=StepStatus.failed,
            error="exit code 1",
        )

        output = format_status(manifest)
        assert "genemark" in output
        assert "completed" in output
        assert "120.5s" in output
        assert "augustus" in output
        assert "failed" in output
        assert "exit code 1" in output


class TestToolVersionSnapshot:
    _CONDA_LIST = (
        "# packages in environment at /env:\n"
        "#\n"
        "# Name    Version    Build    Channel\n"
        "diamond    2.1.11    h43eeafb_0    bioconda\n"
        "hmmer    3.4    hdbdd923_1    bioconda\n"
        "samtools    1.21    h50ea8bc_0    bioconda\n"
        "trinity    2.15.2    pl5321h077b44d_6    bioconda\n"
        "zlib    1.3.1    hb9d3cd8_2    conda-forge\n"
    )

    @pytest.fixture
    def conda_calls(self, tmp_path, monkeypatch):
        """A conda prefix and cache dir under tmp_path, and a recording fake conda."""
        (tmp_path / "env").mkdir()
        monkeypatch.setenv("CONDA_PREFIX", str(tmp_path / "env"))
        monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
        calls = []

        def fake_run(cmd, **_kw):
            calls.append(cmd)
            return subprocess.CompletedProcess(cmd, 0, stdout=self._CONDA_LIST, stderr="")

        monkeypatch.setattr(manifest_mod.subprocess, "run", fake_run)
        return calls

    def test_records_the_conda_packages_tools_toml_declares(self, conda_calls):
        """Not a hard-coded six that missed trinity and diamond and still listed
        hmmer, which pyhmmer replaced."""
        assert manifest_mod._snapshot_tool_versions() == {
            "diamond": "2.1.11", "samtools": "1.21", "trinity": "2.15.2",
        }

    def test_cache_written_for_another_package_list_is_refreshed(self, tmp_path, conda_calls):
        env = tmp_path / "env"
        cache = tmp_path / "cache" / "eukan" / "tool-versions.json"
        cache.parent.mkdir(parents=True)
        cache.write_text(json.dumps({
            "prefix": str(env), "mtime": env.stat().st_mtime, "versions": {"hmmer": "3.4"},
        }))

        assert manifest_mod._snapshot_tool_versions()["trinity"] == "2.15.2"
        assert manifest_mod._snapshot_tool_versions()["trinity"] == "2.15.2"
        assert len(conda_calls) == 1  # the refreshed entry is reused
