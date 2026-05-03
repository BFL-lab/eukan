"""Tests for eukan.manifest — run tracking and step lifecycle."""


from eukan.infra.manifest import (
    SENTINEL,
    RunManifest,
    StepRecord,
    StepStatus,
    clean_interrupted_step,
    format_status,
    is_step_complete,
    is_step_interrupted,
    load_manifest,
    pipeline_step,
    save_manifest,
)


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


class TestIsStepComplete:
    def test_completed_step(self, tmp_path):
        manifest = RunManifest()
        output = tmp_path / "test_output.gff3"
        output.write_text("data")
        manifest.steps["test"] = StepRecord(
            name="test", status=StepStatus.completed,
            output_file=str(output),
        )
        assert is_step_complete(manifest, "test") == output

    def test_incomplete_step(self, tmp_path):
        manifest = RunManifest()
        manifest.steps["test"] = StepRecord(name="test", status=StepStatus.running)
        assert is_step_complete(manifest, "test") is None

    def test_missing_step(self):
        manifest = RunManifest()
        assert is_step_complete(manifest, "nonexistent") is None

    def test_missing_output_file(self, tmp_path):
        manifest = RunManifest()
        manifest.steps["test"] = StepRecord(
            name="test", status=StepStatus.completed,
            output_file=str(tmp_path / "deleted.gff3"),
        )
        assert is_step_complete(manifest, "test") is None


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
