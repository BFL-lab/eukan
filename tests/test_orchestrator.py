"""Tests for eukan.{annotation,assembly}.pipeline — CLI flag → step translation."""

from __future__ import annotations

from typing import ClassVar

from eukan.annotation.pipeline import force_steps_from_run_flags
from eukan.assembly.pipeline import (
    _steps_for,
)
from eukan.assembly.pipeline import (
    force_steps_from_run_flags as assembly_force_steps_from_run_flags,
)
from eukan.infra.manifest import RunManifest, StepRecord, StepStatus
from eukan.infra.pipeline import run_orchestrated_step
from eukan.infra.steps import fingerprint_inputs
from eukan.settings import AssemblyConfig


class TestForceStepsFromRunFlags:
    """``--run-*`` CLI flags translate to the right manifest keys.

    The CLI surface is the dict of booleans accepted as kwargs here;
    the pipeline surface is the list of ``annotation/<step>`` keys
    that get popped from the manifest before re-execution.
    """

    def test_no_flags_returns_empty(self):
        assert force_steps_from_run_flags() == []

    def test_run_genemark_groups_orf_finder(self):
        """--run-genemark forces both genemark and orf_finder (shared flag)."""
        result = force_steps_from_run_flags(run_genemark=True)
        assert set(result) == {"annotation/genemark", "annotation/orf_finder"}

    def test_run_snap_groups_codingquarry(self):
        """--run-snap forces both snap and codingquarry (shared flag)."""
        result = force_steps_from_run_flags(run_snap=True)
        assert set(result) == {"annotation/snap", "annotation/codingquarry"}

    def test_run_augustus_alone(self):
        assert force_steps_from_run_flags(run_augustus=True) == ["annotation/augustus"]

    def test_run_consensus_alone(self):
        assert force_steps_from_run_flags(run_consensus=True) == [
            "annotation/consensus_models"
        ]

    def test_run_prot_align_default_picks_non_ssp(self):
        """Without spaln_ssp, --run-prot-align forces prot_align only."""
        result = force_steps_from_run_flags(run_prot_align=True)
        assert result == ["annotation/prot_align"]

    def test_run_prot_align_with_spsp_picks_ssp(self):
        """With spaln_ssp=True, --run-prot-align forces prot_align_ssp only."""
        result = force_steps_from_run_flags(run_prot_align=True, spaln_ssp=True)
        assert result == ["annotation/prot_align_ssp"]

    def test_spsp_alone_does_not_force(self):
        """spaln_ssp gates which prot-align step gets forced; alone it's a no-op."""
        assert force_steps_from_run_flags(spaln_ssp=True) == []

    def test_all_flags_together(self):
        """Every --run-* flag set forces every step exactly once."""
        result = force_steps_from_run_flags(
            spaln_ssp=False,
            run_genemark=True,
            run_prot_align=True,
            run_augustus=True,
            run_snap=True,
            run_consensus=True,
        )
        assert set(result) == {
            "annotation/genemark",
            "annotation/orf_finder",
            "annotation/prot_align",
            "annotation/augustus",
            "annotation/snap",
            "annotation/codingquarry",
            "annotation/consensus_models",
        }
        assert "annotation/prot_align_ssp" not in result

    def test_all_flags_with_spsp(self):
        """Same as above but spaln_ssp swaps prot_align → prot_align_ssp."""
        result = force_steps_from_run_flags(
            spaln_ssp=True,
            run_genemark=True,
            run_prot_align=True,
            run_augustus=True,
            run_snap=True,
            run_consensus=True,
        )
        assert "annotation/prot_align_ssp" in result
        assert "annotation/prot_align" not in result

    def test_returned_keys_are_prefixed(self):
        """Every returned key carries the ``annotation/`` prefix."""
        result = force_steps_from_run_flags(
            run_genemark=True, run_augustus=True, run_snap=True,
        )
        assert all(k.startswith("annotation/") for k in result)


class TestAssemblyForceStepsFromRunFlags:
    """``--run-X`` on assemble narrows the step list AND forces re-run.

    Returns full ``assembly/<step>`` keys, harmonized with the annotation
    pipeline. Empty list = "run all pending, force nothing".
    """

    _ALL_KEYS: ClassVar[list[str]] = [
        "assembly/minimap2", "assembly/trinity",
        "assembly/jaccard", "assembly/map_transcripts",
        "assembly/strand_correct", "assembly/defuse", "assembly/max_intron_split",
        "assembly/sl_detect", "assembly/sl_cut", "assembly/combinr",
    ]

    def test_no_flags_returns_empty(self):
        """No flags → empty list → run all pending, force nothing."""
        assert assembly_force_steps_from_run_flags() == []

    def test_force_alone_returns_all_keys(self):
        """--force alone → re-run every step from scratch."""
        assert assembly_force_steps_from_run_flags(force=True) == self._ALL_KEYS

    def test_run_minimap2_cascades_to_genome_guided_consumers(self):
        """Re-mapping reads invalidates Trinity (both modes read the aligner BAM) and
        the SL read-side, which together cascade through the whole transcript chain —
        so re-mapping reads forces every step."""
        assert assembly_force_steps_from_run_flags(run_minimap2=True) == self._ALL_KEYS

    def test_run_combinr_alone_forces_combinr_only(self):
        assert assembly_force_steps_from_run_flags(run_combinr=True) == ["assembly/combinr"]

    def test_run_trinity_cascades_through_transcript_chain(self):
        """Re-assembling Trinity (both de novo + genome-guided FASTAs) invalidates the
        whole downstream transcript chain: jaccard re-clips → map_transcripts re-maps →
        strand_correct/defuse/max_intron_split/sl_detect/sl_cut/combinr all re-run."""
        assert assembly_force_steps_from_run_flags(run_trinity=True) == [
            "assembly/trinity", "assembly/jaccard", "assembly/map_transcripts",
            "assembly/strand_correct", "assembly/defuse", "assembly/max_intron_split",
            "assembly/sl_detect", "assembly/sl_cut", "assembly/combinr",
        ]

    def test_run_sl_steps_cascade_downstream(self):
        # max_intron_split is UPSTREAM of the SL steps, so it does not appear here.
        assert assembly_force_steps_from_run_flags(run_sl_detect=True) == [
            "assembly/sl_detect", "assembly/sl_cut", "assembly/combinr"
        ]
        assert assembly_force_steps_from_run_flags(run_sl_cut=True) == [
            "assembly/sl_cut", "assembly/combinr"
        ]

    def test_run_max_intron_split_cascades_downstream(self):
        assert assembly_force_steps_from_run_flags(run_max_intron_split=True) == [
            "assembly/max_intron_split", "assembly/sl_cut", "assembly/combinr"
        ]

    def test_run_strand_correct_cascades_downstream(self):
        assert assembly_force_steps_from_run_flags(run_strand_correct=True) == [
            "assembly/strand_correct", "assembly/defuse", "assembly/max_intron_split",
            "assembly/sl_cut", "assembly/combinr",
        ]

    def test_run_defuse_cascades_downstream(self):
        assert assembly_force_steps_from_run_flags(run_defuse=True) == [
            "assembly/defuse", "assembly/max_intron_split", "assembly/sl_cut",
            "assembly/combinr",
        ]

    def test_run_map_transcripts_cascades_to_sl_and_combinr(self):
        """The new spliced BAM invalidates every step that reads it, directly or not."""
        assert assembly_force_steps_from_run_flags(run_map_transcripts=True) == [
            "assembly/map_transcripts", "assembly/strand_correct", "assembly/defuse",
            "assembly/max_intron_split", "assembly/sl_detect", "assembly/sl_cut",
            "assembly/combinr",
        ]

    def test_run_minimap2_with_force_takes_run_flag(self):
        """--run-minimap2 --force scopes to minimap2's cascade; --run-X takes precedence
        over --force. Its cascade spans the whole pipeline (Trinity reads its BAM)."""
        assert (
            assembly_force_steps_from_run_flags(run_minimap2=True, force=True)
            == self._ALL_KEYS
        )

    def test_multiple_run_flags(self):
        """--run-minimap2 (cascades to everything) plus --run-combinr; union in order."""
        result = assembly_force_steps_from_run_flags(run_minimap2=True, run_combinr=True)
        assert result == self._ALL_KEYS

    def test_step_order_is_pipeline_order(self):
        """Returned keys follow pipeline order regardless of kwarg order."""
        result = assembly_force_steps_from_run_flags(
            run_combinr=True, run_sl_cut=True, run_sl_detect=True,
            run_map_transcripts=True, run_jaccard=True, run_trinity=True,
            run_minimap2=True,
        )
        assert result == self._ALL_KEYS

    def test_returned_keys_are_prefixed(self):
        result = assembly_force_steps_from_run_flags(run_minimap2=True, run_combinr=True)
        assert all(k.startswith("assembly/") for k in result)


class TestRunOrchestratedStepOutput:
    """A declared output_file is recorded in the manifest even when missing,
    so is_step_complete flags it on resume (and the step rebuilds) instead of
    it being lost."""

    def test_completed_step_with_deleted_output_reruns(self, tmp_path):
        """Deleting a completed step's output makes the next run rebuild it,
        no force flag needed (auto-rebuild, not a hard stop)."""
        manifest = RunManifest()
        out = tmp_path / "step" / "expected.txt"
        runs: list[int] = []

        def writes_output():
            runs.append(1)
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text("done")
            return None

        run_orchestrated_step(
            tmp_path, manifest, "annotation/thing", writes_output,
            step_dir=tmp_path / "step", output_file=out,
        )
        out.unlink()  # output lost between runs
        run_orchestrated_step(
            tmp_path, manifest, "annotation/thing", writes_output,
            step_dir=tmp_path / "step", output_file=out,
        )
        assert runs == [1, 1]  # rebuilt rather than trusting the manifest
        assert out.exists()

    def test_missing_declared_output_is_still_recorded(self, tmp_path):
        manifest = RunManifest()
        out = tmp_path / "step" / "expected.gff3"

        def writes_nothing():
            return None

        result = run_orchestrated_step(
            tmp_path, manifest, "annotation/thing",
            writes_nothing,
            step_dir=tmp_path / "step",
            output_file=out,
        )
        assert result == out
        record = manifest.steps["annotation/thing"]
        assert record.output_file == str(out)
        assert record.output_md5 is None  # nothing to checksum

    def test_existing_declared_output_is_recorded_and_checksummed(self, tmp_path):
        manifest = RunManifest()
        out = tmp_path / "step" / "expected.gff3"

        def writes_output():
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text("data")
            return None

        result = run_orchestrated_step(
            tmp_path, manifest, "annotation/thing2",
            writes_output,
            step_dir=tmp_path / "step",
            output_file=out,
        )
        assert result == out
        record = manifest.steps["annotation/thing2"]
        assert record.output_file == str(out)
        assert record.output_md5 is not None


class TestStepDirCleanup:
    """Per-step dirs that the step never wrote into are removed, not left empty."""

    def test_empty_step_dir_removed_after_completion(self, tmp_path):
        manifest = RunManifest()
        sdir = tmp_path / "mystep"
        out = tmp_path / "out.txt"  # written to work_dir, not the step dir

        def fn():
            out.write_text("data")
            return None

        run_orchestrated_step(
            tmp_path, manifest, "assembly/mystep", fn,
            step_dir=sdir, output_file=out,
        )
        assert out.exists()
        assert not sdir.exists()  # empty step dir cleaned up

    def test_used_step_dir_is_kept(self, tmp_path):
        manifest = RunManifest()
        sdir = tmp_path / "mystep"

        def fn():
            (sdir / "work.txt").write_text("data")  # step writes into its own dir
            return None

        run_orchestrated_step(
            tmp_path, manifest, "assembly/mystep", fn, step_dir=sdir,
        )
        assert (sdir / "work.txt").exists()  # non-empty → kept


class TestFingerprintInputs:
    """``fingerprint_inputs`` is a stable, content + presence sensitive digest."""

    def test_none_for_no_inputs(self):
        assert fingerprint_inputs(None) is None
        assert fingerprint_inputs([]) is None

    def test_changes_with_content(self, tmp_path):
        f = tmp_path / "a"
        f.write_text("1")
        before = fingerprint_inputs([f])
        f.write_text("2")
        assert fingerprint_inputs([f]) != before

    def test_missing_differs_from_present(self, tmp_path):
        f = tmp_path / "a"
        missing = fingerprint_inputs([f])
        f.write_text("x")
        assert fingerprint_inputs([f]) != missing

    def test_order_independent(self, tmp_path):
        a = tmp_path / "a"
        a.write_text("1")
        b = tmp_path / "b"
        b.write_text("2")
        assert fingerprint_inputs([a, b]) == fingerprint_inputs([b, a])

    def test_extra_salt_changes_digest(self, tmp_path):
        f = tmp_path / "a"
        f.write_text("1")
        base = fingerprint_inputs([f])
        assert fingerprint_inputs([f], ["max_intron_len=5000"]) != base
        assert fingerprint_inputs([f], ["max_intron_len=2000"]) != fingerprint_inputs(
            [f], ["max_intron_len=5000"]
        )

    def test_extra_only_is_non_none(self):
        assert fingerprint_inputs(None, ["max_intron_len=5000"]) is not None
        assert fingerprint_inputs(None, ["a"]) != fingerprint_inputs(None, ["b"])


class TestInputFingerprintResume:
    """A completed step re-runs on resume only when its declared inputs changed."""

    @staticmethod
    def _run(tmp_path, manifest, inp, runs):
        def fn():
            runs.append(1)
            (tmp_path / "out.txt").write_text("done")
            return None

        return run_orchestrated_step(
            tmp_path, manifest, "assembly/x", fn,
            step_dir=tmp_path / "x",
            output_file=tmp_path / "out.txt",
            input_files=[inp],
        )

    def test_unchanged_inputs_skip_rerun(self, tmp_path):
        manifest = RunManifest()
        inp = tmp_path / "in.txt"
        inp.write_text("v1")
        runs: list[int] = []
        self._run(tmp_path, manifest, inp, runs)
        self._run(tmp_path, manifest, inp, runs)
        assert runs == [1]  # second invocation reused the cached output

    def test_changed_inputs_force_rerun(self, tmp_path):
        manifest = RunManifest()
        inp = tmp_path / "in.txt"
        inp.write_text("v1")
        runs: list[int] = []
        self._run(tmp_path, manifest, inp, runs)
        inp.write_text("v2")  # upstream output changed
        self._run(tmp_path, manifest, inp, runs)
        assert runs == [1, 1]  # re-ran on the changed input

    def test_legacy_record_without_fingerprint_is_reused(self, tmp_path):
        """A completed step recorded before input_md5 existed is not treated stale."""
        manifest = RunManifest()
        out = tmp_path / "out.txt"
        out.write_text("done")
        manifest.steps["assembly/x"] = StepRecord(
            name="assembly/x", status=StepStatus.completed,
            output_file=str(out), input_md5=None,
        )
        runs: list[int] = []
        self._run(tmp_path, manifest, tmp_path / "in.txt", runs)
        assert runs == []  # reused despite a (missing) declared input


class TestAssemblyStepScalars:
    """The enforcement steps fold ``max_intron_len`` into their resume fingerprint
    so a changed ``-M`` re-runs them (but not the expensive mappers)."""

    @staticmethod
    def _scalars(tmp_path, name, max_intron=5000):
        cfg = AssemblyConfig(
            genome=tmp_path / "g.fa", work_dir=tmp_path, num_cpu=2,
            max_intron_len=max_intron,
        )
        spec = next(s for s in _steps_for() if s.name == name)
        return spec.scalars(cfg) if spec.scalars else None

    def test_jaccard_has_no_scalar(self, tmp_path):
        # jaccard has no declared output, so it always re-runs on resume and a
        # scalar would never be consulted (is_step_complete short-circuits first);
        # a changed greediness re-clips via the always-re-run + map_transcripts.
        assert self._scalars(tmp_path, "jaccard") is None

    def test_combinr_tracks_max_intron_and_stringent_overlap(self, tmp_path):
        scalars = self._scalars(tmp_path, "combinr")
        assert "max_intron_len=5000" in scalars
        assert any(s.startswith("combinr_stringent_overlap=") for s in scalars)

    def test_max_intron_split_tracks_max_intron_and_min_fragment(self, tmp_path):
        # The model-level -M enforcement lives here now (a declared output makes the
        # scalar live), so changing -M re-runs max_intron_split and cascades on.
        scalars = self._scalars(tmp_path, "max_intron_split")
        assert "max_intron_len=5000" in scalars
        assert any(s.startswith("min_sl_fragment=") for s in scalars)

    def test_sl_cut_tracks_min_fragment_not_max_intron(self, tmp_path):
        # sl_cut no longer owns the max-intron split, so it tracks only min_sl_fragment.
        scalars = self._scalars(tmp_path, "sl_cut")
        assert any(s.startswith("min_sl_fragment=") for s in scalars)
        assert not any(s.startswith("max_intron_len=") for s in scalars)

    def test_defuse_tracks_its_knobs(self, tmp_path):
        scalars = self._scalars(tmp_path, "defuse")
        assert any(s.startswith("defuse=") for s in scalars)
        assert any(s.startswith("defuse_overlap_tolerance=") for s in scalars)

    def test_mappers_do_not_track_max_intron(self, tmp_path):
        # minimap2 bounds introns at map time via -G, not a resume scalar; tightening
        # -M re-runs only the post-mapping steps (max_intron_split/combinr), not a
        # multi-hour re-map. Trinity likewise emits FASTAs with no -M enforcement.
        for name in ("minimap2", "trinity", "map_transcripts"):
            assert self._scalars(tmp_path, name) is None


class TestScalarFingerprintResume:
    """A completed step re-runs on resume when a tracked scalar input changed,
    even though its files are byte-identical (e.g. ``-M/--max-intron``)."""

    @staticmethod
    def _run(tmp_path, manifest, scalars, runs):
        def fn():
            runs.append(1)
            (tmp_path / "out.txt").write_text("done")
            return None

        return run_orchestrated_step(
            tmp_path, manifest, "assembly/x", fn,
            step_dir=tmp_path / "x",
            output_file=tmp_path / "out.txt",
            input_scalars=scalars,
        )

    def test_unchanged_scalar_skips_rerun(self, tmp_path):
        manifest = RunManifest()
        runs: list[int] = []
        self._run(tmp_path, manifest, ["max_intron_len=5000"], runs)
        self._run(tmp_path, manifest, ["max_intron_len=5000"], runs)
        assert runs == [1]

    def test_changed_scalar_forces_rerun(self, tmp_path):
        manifest = RunManifest()
        runs: list[int] = []
        self._run(tmp_path, manifest, ["max_intron_len=5000"], runs)
        self._run(tmp_path, manifest, ["max_intron_len=2000"], runs)
        assert runs == [1, 1]


class TestRerunClearsStaleOutput:
    """A re-run of a completed step clears its declared output first, so step
    functions that reuse a complete output as their own resume check (the
    transcript mappers) actually regenerate it rather than reuse a stale file."""

    @staticmethod
    def _mapper_like(out, runs):
        """A step fn that, like the mappers, reuses a complete output if present."""

        def fn():
            if out.exists():
                return None
            runs.append(1)
            out.write_text("mapped")
            return None

        return fn

    def test_forced_rerun_clears_output(self, tmp_path):
        manifest = RunManifest()
        out = tmp_path / "out.bam"
        runs: list[int] = []
        fn = self._mapper_like(out, runs)
        run_orchestrated_step(
            tmp_path, manifest, "assembly/x", fn,
            step_dir=tmp_path / "x", output_file=out,
        )
        run_orchestrated_step(
            tmp_path, manifest, "assembly/x", fn,
            step_dir=tmp_path / "x", output_file=out, force=True,
        )
        assert runs == [1, 1]  # forced re-run cleared the BAM and re-mapped

    def test_stale_input_rerun_clears_output(self, tmp_path):
        manifest = RunManifest()
        inp = tmp_path / "in.fasta"
        inp.write_text("v1")
        out = tmp_path / "out.bam"
        runs: list[int] = []
        fn = self._mapper_like(out, runs)
        run_orchestrated_step(
            tmp_path, manifest, "assembly/x", fn,
            step_dir=tmp_path / "x", output_file=out, input_files=[inp],
        )
        inp.write_text("v2")  # the query the mapper reads changed
        run_orchestrated_step(
            tmp_path, manifest, "assembly/x", fn,
            step_dir=tmp_path / "x", output_file=out, input_files=[inp],
        )
        assert runs == [1, 1]  # stale BAM cleared, so the mapper re-ran
