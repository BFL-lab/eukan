"""Tests for eukan.{annotation,assembly}.pipeline — CLI flag → step translation."""

from __future__ import annotations

import itertools
import logging
import shutil
from dataclasses import replace
from pathlib import Path
from typing import ClassVar

import pytest

from eukan.annotation import pipeline as annotation_pipeline
from eukan.annotation.pipeline import force_steps_from_run_flags
from eukan.assembly import pipeline as assembly_pipeline
from eukan.assembly.pipeline import (
    STEPS,
)
from eukan.assembly.pipeline import (
    force_steps_from_run_flags as assembly_force_steps_from_run_flags,
)
from eukan.cli.status import format_status
from eukan.infra.manifest import (
    ANNOTATION,
    ASSEMBLY,
    MANIFEST_FILE,
    REPEATS,
    RunManifest,
    StepRecord,
    StepStatus,
    load_manifest,
    step_key,
)
from eukan.infra.pipeline import (
    StepSpec,
    run_orchestrated_step,
    run_simple_pipeline,
)
from eukan.infra.pipeline import (
    force_steps_from_run_flags as generic_force_steps_from_run_flags,
)
from eukan.infra.steps import fingerprint_inputs, is_step_complete
from eukan.infra.utils import md5_file
from eukan.settings import AssemblyConfig, Kingdom, PipelineConfig, RepeatsConfig


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


class TestGenericForceStepsFromRunFlags:
    """The shared helper forces a flagged step and, through *downstream*, every
    step that reads its output, transitively, in pipeline order."""

    STEPS: ClassVar[list[StepSpec]] = [
        StepSpec(name, lambda config: None) for name in ("a", "b", "c", "d")
    ]

    def test_downstream_is_followed_transitively(self):
        downstream = {"b": ("c",), "c": ("d",)}
        assert generic_force_steps_from_run_flags(
            ASSEMBLY, self.STEPS, downstream=downstream, run_b=True,
        ) == ["assembly/b", "assembly/c", "assembly/d"]

    def test_without_downstream_only_flagged_steps(self):
        assert generic_force_steps_from_run_flags(
            ASSEMBLY, self.STEPS, run_d=True, run_a=True,
        ) == ["assembly/a", "assembly/d"]


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


class TestAssemblyForceStepsMatchFlagDict:
    """For every combination of the ten ``--run-X`` flags and ``--force``, the
    generic driver with ``downstream=`` forces the steps the hand-written flag dict
    and ``_expand_downstream`` forced before (f7c164c), frozen here as the oracle."""

    ORDER = (
        "minimap2", "trinity", "jaccard", "map_transcripts", "strand_correct",
        "defuse", "max_intron_split", "sl_detect", "sl_cut", "combinr",
    )
    OLD_DOWNSTREAM: ClassVar[dict[str, tuple[str, ...]]] = {
        "minimap2": ("trinity", "sl_detect"),
        "trinity": ("jaccard",),
        "jaccard": ("map_transcripts",),
        "map_transcripts": ("strand_correct", "sl_detect"),
        "strand_correct": ("defuse",),
        "defuse": ("max_intron_split",),
        "max_intron_split": ("sl_cut",),
        "sl_detect": ("sl_cut",),
        "sl_cut": ("combinr",),
    }

    def _old_answer(self, flags: dict[str, bool], force: bool) -> list[str]:
        selected = {name for name, on in flags.items() if on}
        if selected:
            expanded, stack = set(selected), list(selected)
            while stack:
                for dep in self.OLD_DOWNSTREAM.get(stack.pop(), ()):
                    if dep not in expanded:
                        expanded.add(dep)
                        stack.append(dep)
            return [f"assembly/{n}" for n in self.ORDER if n in expanded]
        if force:
            return [f"assembly/{n}" for n in self.ORDER]
        return []

    def test_every_flag_combination(self):
        assert tuple(s.name for s in STEPS) == self.ORDER
        for *on, force in itertools.product((False, True), repeat=len(self.ORDER) + 1):
            flags = dict(zip(self.ORDER, on, strict=True))
            got = assembly_force_steps_from_run_flags(
                force=force, **{f"run_{name}": v for name, v in flags.items()},
            )
            assert got == self._old_answer(flags, force), (flags, force)


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
        assert record.output_file == "step/expected.gff3"  # relative to the run dir
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
        assert record.output_file == "step/expected.gff3"
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


class TestOutputlessStepResume:
    """A step with no declared output resumes on its recorded input fingerprint:
    it is skipped while that matches and re-runs when it does not; one that
    recorded no fingerprint re-runs every time."""

    @staticmethod
    def _run(tmp_path, manifest, fn, **fingerprint):
        return run_orchestrated_step(
            tmp_path, manifest, "assembly/x", fn, step_dir=tmp_path / "x", **fingerprint,
        )

    def test_unchanged_fingerprint_skips(self, tmp_path):
        manifest = RunManifest()
        inp = tmp_path / "in.txt"
        inp.write_text("v1")
        runs: list[int] = []
        for _ in range(2):
            result = self._run(
                tmp_path, manifest, lambda: runs.append(1),
                input_files=[inp], input_scalars=["knob=1"],
            )
            assert result is None
        assert runs == [1]

    def test_changed_input_or_scalar_reruns(self, tmp_path):
        manifest = RunManifest()
        inp = tmp_path / "in.txt"
        inp.write_text("v1")
        runs: list[int] = []
        for content, knob in (("v1", 1), ("v2", 1), ("v2", 2)):
            inp.write_text(content)
            self._run(
                tmp_path, manifest, lambda: runs.append(1),
                input_files=[inp], input_scalars=[f"knob={knob}"],
            )
        assert runs == [1, 1, 1]

    def test_no_fingerprint_always_reruns(self, tmp_path):
        manifest = RunManifest()
        runs: list[int] = []
        for _ in range(2):
            self._run(tmp_path, manifest, lambda: runs.append(1))
        assert runs == [1, 1]

    def test_removed_written_file_reruns(self, tmp_path):
        # The step names the file it writes among its inputs; the fingerprint is
        # taken after it runs, so removing that file re-runs the step.
        manifest = RunManifest()
        written = tmp_path / "written.txt"
        runs: list[int] = []

        def fn():
            runs.append(1)
            written.write_text("models")

        for _ in range(2):
            self._run(tmp_path, manifest, fn, input_files=[written])
        written.unlink()
        self._run(tmp_path, manifest, fn, input_files=[written])
        assert runs == [1, 1]
        assert written.exists()


class TestAssemblyStepScalars:
    """Every step that consumes ``max_intron_len`` folds it into its resume
    fingerprint, so a changed ``-M`` re-runs it."""

    @staticmethod
    def _scalars(tmp_path, name, max_intron=5000):
        cfg = AssemblyConfig(
            genome=tmp_path / "g.fa", work_dir=tmp_path, num_cpu=2,
            max_intron_len=max_intron,
        )
        spec = next(s for s in STEPS if s.name == name)
        return spec.scalars(cfg) if spec.scalars else None

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

    def test_mappers_and_trinity_track_max_intron(self, tmp_path):
        # minimap2 maps reads and transcripts with -G and Trinity genome-guided takes
        # -M as --genome_guided_max_intron, so all three track it.
        for name in ("minimap2", "trinity", "map_transcripts"):
            assert self._scalars(tmp_path, name) == ["max_intron_len=5000"]


class TestMaxIntronConsumersRerun:
    """Through the driver, a changed ``-M`` on resume re-runs the read mapping, the
    Trinity assembly and the transcript mapping; an unchanged one reuses them."""

    @pytest.mark.parametrize("name", ["minimap2", "trinity", "map_transcripts"])
    def test_changed_max_intron_reruns(self, tmp_path, name):
        spec = next(s for s in STEPS if s.name == name)
        runs: list[int] = []

        def fake(config):
            runs.append(config.max_intron_len)
            (config.work_dir / spec.output).write_text("out")

        for max_intron in (5000, 5000, 2000):
            cfg = AssemblyConfig(
                genome=tmp_path / "g.fa", work_dir=tmp_path, num_cpu=2,
                max_intron_len=max_intron,
            )
            run_simple_pipeline(ASSEMBLY, [replace(spec, fn=fake)], cfg)
        assert runs == [5000, 2000]


class TestDefuseTracksProteinDb:
    """defuse blasts every transcript against the --uniprot DB, so a changed DB
    re-runs it on resume, as it re-runs strand_correct. On a stranded library
    strand_correct's models do not change with the DB, so nothing else would."""

    def test_changed_uniprot_db_reruns_defuse(self, tmp_path):
        spec = next(s for s in STEPS if s.name == "defuse")
        db = tmp_path / "uniprot.dmnd"
        db.write_text("v1")
        runs: list[int] = []

        def fake(config):
            runs.append(1)
            (config.work_dir / spec.output).write_text(_GFF3)

        cfg = AssemblyConfig(
            genome=tmp_path / "g.fa", work_dir=tmp_path, num_cpu=2,
            uniprot_db=db, defuse=True,
        )
        for _ in range(2):
            run_simple_pipeline(ASSEMBLY, [replace(spec, fn=fake)], cfg)
        db.write_text("v2")
        run_simple_pipeline(ASSEMBLY, [replace(spec, fn=fake)], cfg)
        assert runs == [1, 1]


class TestOutputlessAssemblyStepsResume:
    """Through the driver, the assemble steps with no declared output (jaccard,
    strand_correct, sl_detect, sl_cut) are skipped on an unchanged resume and re-run
    when a knob they read, a file they read, or a file they write has changed."""

    # step -> (a file it reads, a file it writes)
    FILES: ClassVar[dict[str, tuple[str, str]]] = {
        "jaccard": ("trinity-gg.fasta", "trinity-denovo.jaccard.fasta"),
        "strand_correct": ("trinity-gg.genome.bam", "trinity-denovo.genome.stranded.gff3"),
        "sl_detect": ("minimap2_Aligned.sortedByCoord.out.bam", "sl_acceptors.gff3"),
        "sl_cut": ("trinity-gg.genome.maxintron.gff3", "trinity-denovo.genome.cut.gff3"),
    }
    # step -> every config knob the step function reads, with a changed value
    KNOBS: ClassVar[dict[str, dict[str, object]]] = {
        "jaccard": {
            "jaccard_greediness": 2.0, "jaccard_max_trough": 0.1,
            "jaccard_min_delta": 0.5, "jaccard_max_adaptive_trough": 0.4,
            "min_sl_fragment": 30,
        },
        "strand_correct": {
            "strand_specific": "RF", "genetic_code": "6",
            "strand_blastx_evalue": 1e-3, "min_strand_consensus": 10,
        },
        "sl_detect": {
            "sl_sequence": "GTTTCTGTACTTTATT", "sl_adapter_filter": False,
            "adapter_sequences": ["ACGTACGTACGT"], "min_sl_clip_len": 12,
            "min_sl_insertion_len": 12, "sl_cluster_window": 10,
        },
        "sl_cut": {"min_sl_fragment": 30},
    }

    def _driver(self, tmp_path, name):
        spec = next(s for s in STEPS if s.name == name)
        read, written = self.FILES[name]
        (tmp_path / read).write_text("v1")
        runs: list[int] = []

        def fake(config):
            runs.append(1)
            (config.work_dir / written).write_text("out")

        def run(**knobs):
            cfg = AssemblyConfig(
                genome=tmp_path / "g.fa", work_dir=tmp_path, num_cpu=2, **knobs,
            )
            run_simple_pipeline(ASSEMBLY, [replace(spec, fn=fake)], cfg)
            return len(runs)

        return run, tmp_path / read, tmp_path / written

    @pytest.mark.parametrize("name", list(FILES))
    def test_resume_follows_files_read_and_written(self, tmp_path, name):
        run, read, written = self._driver(tmp_path, name)
        assert run() == 1
        assert run() == 1  # unchanged: skipped
        read.write_text("v2")
        assert run() == 2  # a file it reads changed
        written.unlink()
        assert run() == 3  # a file it writes was removed
        assert run() == 3

    @pytest.mark.parametrize(
        ("name", "knob"), [(n, k) for n, knobs in KNOBS.items() for k in knobs],
    )
    def test_changed_knob_reruns(self, tmp_path, name, knob):
        run, _, _ = self._driver(tmp_path, name)
        assert run() == 1
        assert run() == 1
        assert run(**{knob: self.KNOBS[name][knob]}) == 2

    def test_strand_correct_follows_uniprot_db(self, tmp_path):
        run, _, _ = self._driver(tmp_path, "strand_correct")
        db = tmp_path / "sprot.faa"
        db.write_text(">a\nMK\n")
        assert run() == 1
        assert run(uniprot_db=db) == 2  # DB given
        assert run(uniprot_db=db) == 2
        db.write_text(">a\nMKV\n")
        assert run(uniprot_db=db) == 3  # DB content changed


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

    def test_rerun_clears_sibling_outputs(self, tmp_path):
        # A step fn that skips each sub-output already present (Trinity's per-mode
        # FASTAs, the per-set transcript BAMs) regenerates them all on a re-run.
        manifest = RunManifest()
        out, sibling = tmp_path / "out.bam", tmp_path / "sibling.bam"
        made: list[str] = []

        def fn():
            for p in (sibling, out):
                if not p.exists():
                    made.append(p.name)
                    p.write_text("mapped")

        for force in (False, True):
            run_orchestrated_step(
                tmp_path, manifest, "assembly/x", fn, step_dir=tmp_path / "x",
                output_file=out, clear_files=[sibling], force=force,
            )
        assert made == ["sibling.bam", "out.bam", "sibling.bam", "out.bam"]

    def test_failed_step_resume_keeps_sibling_outputs(self, tmp_path):
        # Only a re-run of a *completed* step clears: a step that failed part-way
        # keeps the sub-output it finished, so the resume does not redo it.
        manifest = RunManifest()
        out, sibling = tmp_path / "out.bam", tmp_path / "sibling.bam"
        made: list[str] = []

        def fn(fail):
            for p in (sibling, out):
                if not p.exists():
                    if fail and p == out:
                        raise RuntimeError("killed")
                    made.append(p.name)
                    p.write_text("mapped")

        for fail in (True, False):
            try:
                run_orchestrated_step(
                    tmp_path, manifest, "assembly/x", fn, fail, step_dir=tmp_path / "x",
                    output_file=out, clear_files=[sibling],
                )
            except RuntimeError:
                assert fail
        assert made == ["sibling.bam", "out.bam"]


_GFF3 = "##gff-version 3\nchr1\tt\tgene\t1\t10\t.\t+\t.\tID=g1\n"


class TestRerunFromScratch:
    """A forced step (whatever its record says) or a completed one whose inputs
    changed re-runs from scratch: its declared output, clear_files and step dir go
    first, so a cache a step function keeps in its step dir (GeneMark's
    genemark.gtf, the AUGUSTUS training set) cannot short-circuit it."""

    @staticmethod
    def _cached_step(sdir, ran):
        """A step fn that, like GeneMark, reuses a cache in its step dir."""

        def fn():
            cache = sdir / "genemark.gtf"
            ran.append("reused" if cache.exists() else "fresh")
            cache.write_text("gtf")
            (sdir / "genemark.gff3").write_text(_GFF3)
            return sdir / "genemark.gff3"

        return fn

    def test_forced_rerun_clears_step_dir(self, tmp_path):
        manifest = RunManifest()
        sdir = tmp_path / "genemark"
        ran: list[str] = []
        fn = self._cached_step(sdir, ran)
        for force in (False, True):
            run_orchestrated_step(
                tmp_path, manifest, "annotation/genemark", fn, step_dir=sdir, force=force,
            )
        assert ran == ["fresh", "fresh"]

    def test_changed_inputs_rerun_clears_step_dir(self, tmp_path):
        manifest = RunManifest()
        sdir = tmp_path / "genemark"
        hints = tmp_path / "hints.gff"
        ran: list[str] = []
        fn = self._cached_step(sdir, ran)
        for content in ("v1", "v2"):
            hints.write_text(content)
            run_orchestrated_step(
                tmp_path, manifest, "annotation/genemark", fn, step_dir=sdir,
                input_files=[hints],
            )
        assert ran == ["fresh", "fresh"]

    def test_missing_output_rerun_keeps_step_dir(self, tmp_path):
        # The output was lost, not invalidated: the cache still matches the inputs.
        manifest = RunManifest()
        sdir = tmp_path / "genemark"
        ran: list[str] = []
        fn = self._cached_step(sdir, ran)
        run_orchestrated_step(tmp_path, manifest, "annotation/genemark", fn, step_dir=sdir)
        (sdir / "genemark.gff3").unlink()
        run_orchestrated_step(tmp_path, manifest, "annotation/genemark", fn, step_dir=sdir)
        assert ran == ["fresh", "reused"]

    def test_forced_rerun_keeps_manifest_dir(self, tmp_path):
        manifest = RunManifest()
        keep = tmp_path / "keep.txt"
        keep.write_text("x")
        for force in (False, True):
            run_orchestrated_step(
                tmp_path, manifest, "repeats/x", lambda: None, step_dir=tmp_path, force=force,
            )
        assert keep.exists()

    def test_forced_rerun_of_failed_step_clears_outputs(self, tmp_path):
        manifest = RunManifest()
        out, sibling = tmp_path / "out.bam", tmp_path / "sibling.bam"
        made: list[str] = []

        def fn(fail):
            for p in (sibling, out):
                if not p.exists():
                    if fail and p == out:
                        raise RuntimeError("killed")
                    made.append(p.name)
                    p.write_text("mapped")

        with pytest.raises(RuntimeError):
            run_orchestrated_step(
                tmp_path, manifest, "assembly/x", fn, True, step_dir=tmp_path / "x",
                output_file=out, clear_files=[sibling],
            )
        run_orchestrated_step(
            tmp_path, manifest, "assembly/x", fn, False, step_dir=tmp_path / "x",
            output_file=out, clear_files=[sibling], force=True,
        )
        assert made == ["sibling.bam", "sibling.bam", "out.bam"]

    def test_forced_annotation_step_clears_its_cache(self, tmp_path, monkeypatch):
        # --run-genemark: run_annotation_pipeline forces the step through
        # run_orchestrated_step rather than dropping its manifest record, so the
        # genemark.gtf GeneMark keeps in its step dir is cleared too.
        genome = tmp_path / "genome.fa"
        genome.write_text(">chr1\nACGTACGTAC\n")
        config = PipelineConfig(
            genome=genome, proteins=[tmp_path / "proteins.faa"], work_dir=tmp_path,
            kingdom="animal",
        )
        ran: list[str] = []

        def other_step(name):
            def fn(config, *args, **kwargs):
                out = config.work_dir / name / f"{name}.gff3"
                out.write_text(_GFF3)
                return out

            return fn

        monkeypatch.setattr(
            annotation_pipeline, "run_genemark",
            lambda config: self._cached_step(config.work_dir / "genemark", ran)(),
        )
        for attr, name in (
            ("align_proteins", "prot_align"), ("run_augustus", "augustus"),
            ("run_snap", "snap"), ("build_consensus_models", "consensus_models"),
        ):
            monkeypatch.setattr(annotation_pipeline, attr, other_step(name))

        annotation_pipeline.run_annotation_pipeline(config)
        annotation_pipeline.run_annotation_pipeline(
            config, force_steps=[step_key(ANNOTATION, "genemark")],
        )
        assert ran == ["fresh", "fresh"]


class TestMovedRunDir:
    """A run dir moved elsewhere resumes where it left off: the manifest records
    outputs, and keys input fingerprints, relative to the run dir."""

    @staticmethod
    def _step(run_dir, manifest, runs):
        out = run_dir / "x" / "out.txt"

        def fn():
            runs.append(1)
            out.write_text("done")

        return run_orchestrated_step(
            run_dir, manifest, "assembly/x", fn,
            step_dir=run_dir / "x", output_file=out,
            input_files=[run_dir / "in.txt"], input_scalars=["k=1"],
        )

    def test_moved_run_dir_resumes(self, tmp_path):
        run_a, run_b = tmp_path / "runA", tmp_path / "runB"
        run_a.mkdir()
        (run_a / "in.txt").write_text("v1")
        runs: list[int] = []
        self._step(run_a, RunManifest(), runs)
        shutil.move(run_a, run_b)

        assert self._step(run_b, load_manifest(run_b), runs) == run_b / "x" / "out.txt"
        assert runs == [1]
        assert (run_b / "x" / "out.txt").read_text() == "done"

    def test_moved_assembly_run_dir_reruns_nothing(self, tmp_path, monkeypatch):
        # Before, every assemble step re-ran after a move, its outputs deleted first.
        ran: list[str] = []
        specs = [replace(s, fn=lambda config, _n=s.name: ran.append(_n)) for s in STEPS]
        monkeypatch.setattr(assembly_pipeline, "STEPS", tuple(specs))

        def config(run_dir):
            return AssemblyConfig(
                genome=run_dir / "g.fa", left_reads=run_dir / "l.fq",
                right_reads=run_dir / "r.fq", work_dir=run_dir / "assemble",
                manifest_dir=run_dir, num_cpu=2,
            )

        run_a, run_b = tmp_path / "runA", tmp_path / "runB"
        cfg = config(run_a)
        cfg.ensure_dirs()
        for name in ("g.fa", "l.fq", "r.fq"):
            (run_a / name).write_text("x")
        for spec in specs:
            for name in [*([spec.output] if spec.output else []), *spec.clears]:
                (cfg.work_dir / name).write_text(_GFF3 if name.endswith(".gff3") else "x")
        assembly_pipeline.run_assembly(cfg)
        assert ran
        files = [p.relative_to(run_a) for p in run_a.rglob("*") if p.name != MANIFEST_FILE]
        shutil.move(run_a, run_b)
        ran.clear()

        assembly_pipeline.run_assembly(config(run_b))

        assert ran == []
        assert all((run_b / f).exists() for f in files)

    def test_record_from_before_run_dir_relative_paths_is_reused(self, tmp_path):
        # Older manifests hold the absolute output path and a fingerprint keyed by
        # absolute input path; in place, the step is still complete.
        inp = tmp_path / "in.txt"
        inp.write_text("v1")
        out = tmp_path / "x" / "out.txt"
        out.parent.mkdir()
        out.write_text("done")
        manifest = RunManifest(steps={
            "assembly/x": StepRecord(
                name="assembly/x", status=StepStatus.completed, output_file=str(out),
                output_md5=md5_file(out), input_md5=fingerprint_inputs([inp], ["k=1"]),
            ),
            "assembly/y": StepRecord(
                name="assembly/y", status=StepStatus.completed,
                input_md5=fingerprint_inputs([inp]),
            ),
        })
        runs: list[int] = []

        assert self._step(tmp_path, manifest, runs) == out
        assert runs == []
        assert out.read_text() == "done"
        assert is_step_complete(tmp_path, manifest, "assembly/y", [inp]) is True


class _FakeAnnotation:
    """Stub the annotation step functions of a protein-only animal run (GeneMark,
    protein alignment, AUGUSTUS, SNAP, consensus): each records its step name in
    ``ran`` and writes a GFF3 into its step dir — ``content[name]`` if set — or
    raises if its name is in ``fail``."""

    STEPS: ClassVar[dict[str, str]] = {
        "run_genemark": "genemark", "align_proteins": "prot_align",
        "run_augustus": "augustus", "run_snap": "snap",
        "build_consensus_models": "consensus_models",
    }

    def __init__(self, monkeypatch, tmp_path):
        self.ran: list[str] = []
        self.content: dict[str, str] = {}
        self.fail: set[str] = set()
        genome = tmp_path / "genome.fa"
        genome.write_text(">chr1\nACGTACGTAC\n")
        self.proteins = tmp_path / "proteins.faa"
        self.proteins.write_text(">p\nMK\n")
        self.config = PipelineConfig(
            genome=genome, proteins=[self.proteins], work_dir=tmp_path, kingdom="animal",
        )
        for attr, name in self.STEPS.items():
            monkeypatch.setattr(annotation_pipeline, attr, self._step(name))

    def _step(self, name):
        def fn(config, *args, **kwargs):
            self.ran.append(name)
            if name in self.fail:
                raise RuntimeError(f"{name} failed")
            out = config.work_dir / name / f"{name}.gff3"
            out.write_text(self.content.get(name, _GFF3))
            return out

        return fn

    def run(self, *forced: str) -> list[str]:
        """Run the pipeline, forcing the named steps; return the steps that ran."""
        self.ran.clear()
        annotation_pipeline.run_annotation_pipeline(
            self.config, force_steps=[step_key(ANNOTATION, n) for n in forced],
        )
        return list(self.ran)


class TestAnnotationResume:
    """Each annotation step decides on its own record whether it re-runs."""

    def test_nothing_left_runs_no_step(self, tmp_path, monkeypatch):
        fake = _FakeAnnotation(monkeypatch, tmp_path)
        assert fake.run() == ["genemark", "prot_align", "augustus", "snap", "consensus_models"]
        assert fake.run() == []

    def test_failed_step_reruns_on_resume(self, tmp_path, monkeypatch):
        # A forced re-run that failed leaves the step's record failed while the
        # consensus record from the earlier run is still complete; the resume must
        # re-run the failed step rather than report every step complete.
        fake = _FakeAnnotation(monkeypatch, tmp_path)
        fake.run()
        fake.fail.add("augustus")
        with pytest.raises(RuntimeError):
            fake.run("augustus")
        fake.fail.clear()
        assert fake.run() == ["augustus"]


class TestAnnotationInputFingerprints:
    """Each annotation step fingerprints the files it reads, so one that changed
    since the step completed re-runs it, and the re-run cascades on through the
    steps whose inputs it rewrites."""

    def test_rerun_genemark_reruns_its_consumers(self, tmp_path, monkeypatch):
        fake = _FakeAnnotation(monkeypatch, tmp_path)
        fake.run()
        fake.content["genemark"] = _GFF3.replace("ID=g1", "ID=g2")
        # SNAP and the consensus read the AUGUSTUS and alignment outputs, which came
        # out unchanged; combinr is never given the GeneMark models.
        assert fake.run("genemark") == ["genemark", "prot_align", "augustus"]
        assert fake.run() == []

    def test_unchanged_genemark_output_reruns_nothing_else(self, tmp_path, monkeypatch):
        fake = _FakeAnnotation(monkeypatch, tmp_path)
        fake.run()
        assert fake.run("genemark") == ["genemark"]

    def test_edited_proteins_rerun_the_alignment(self, tmp_path, monkeypatch):
        fake = _FakeAnnotation(monkeypatch, tmp_path)
        fake.run()
        fake.proteins.write_text(">p\nMKV\n")
        assert fake.run() == ["prot_align"]

    def test_new_repeat_hints_rerun_augustus(self, tmp_path, monkeypatch):
        from eukan.infra.artifacts import Artifact

        fake = _FakeAnnotation(monkeypatch, tmp_path)
        fake.run()
        (tmp_path / Artifact.REPEATMASK_HINTS.value).write_text("hints\n")
        assert fake.run() == ["augustus"]

    @pytest.mark.parametrize(("field", "value", "reruns"), [
        ("genetic_code", "4", ["genemark", "prot_align", "augustus", "consensus_models"]),
        ("allow_noncanonical_splice", True, ["augustus"]),
        # GeneMark reads only whether the kingdom is fungus; a fungus run adds CodingQuarry.
        ("kingdom", Kingdom.fungus, ["genemark", "codingquarry", "consensus_models"]),
        ("kingdom", Kingdom.plant, []),
        ("combinr_path", Path("/opt/combinr/combinr"), ["consensus_models"]),
        ("orf_patch", "all", ["consensus_models"]),
    ])
    def test_changed_scalar_reruns_the_steps_that_read_it(
        self, tmp_path, monkeypatch, field, value, reruns,
    ):
        fake = _FakeAnnotation(monkeypatch, tmp_path)
        monkeypatch.setattr(annotation_pipeline, "run_codingquarry", fake._step("codingquarry"))
        fake.run()
        fake.config = fake.config.model_copy(update={field: value})
        assert fake.run() == reruns

    def test_transcript_orfs_rerun_the_consensus_only_under_orf_patch(
        self, tmp_path, monkeypatch,
    ):
        fake = _FakeAnnotation(monkeypatch, tmp_path)
        orfs = tmp_path / "orf_finder" / "transcript_orfs.gff3"
        orfs.parent.mkdir()
        fake.run()
        orfs.write_text(_GFF3)
        assert fake.run() == []  # the patch is off: the consensus does not read them
        fake.config = fake.config.model_copy(update={"orf_patch": "all"})
        assert fake.run() == ["consensus_models"]
        orfs.write_text(_GFF3.replace("ID=g1", "ID=g2"))
        assert fake.run() == ["consensus_models"]


class TestTrinityTrackRerun:
    """A re-run of the completed trinity / map_transcripts steps regenerates the
    genome-guided output too, not only the declared de novo one: run_trinity skips
    a mode whose FASTA exists and the mapper reuses a per-set BAM that passes
    quickcheck, so the driver must clear both before re-running."""

    @staticmethod
    def _stub_tools(monkeypatch, wd, ran):
        from eukan.assembly import minimap2, trinity

        def fake_trinity(cmd, **kw):
            prefix = cmd[cmd.index("--output") + 1]
            ran.append(prefix)
            name = "Trinity-GG.fasta" if prefix == "trinity-gg" else "Trinity.fasta"
            (wd / f"{prefix}.{name}").write_text(">t\nACGT\n")

        def fake_piped(cmd1, cmd2, **kw):
            if cmd1[0] == "minimap2":
                ran.append(cmd1[-1].rsplit("/", 1)[-1])  # the transcript FASTA mapped
            (wd / cmd2[cmd2.index("-o") + 1]).write_text("bam")
            return ""

        monkeypatch.setattr(trinity, "run_cmd", fake_trinity)
        monkeypatch.setattr(minimap2, "run_piped", fake_piped)
        monkeypatch.setattr(minimap2, "run_cmd", lambda *a, **k: None)
        monkeypatch.setattr(minimap2, "_write_unmapped_fasta", lambda u, o: o.write_text(">u\nA\n"))
        monkeypatch.setattr(minimap2, "_finalize_transcript_diagnostics", lambda c: None)
        # A BAM from a completed run passes samtools quickcheck.
        monkeypatch.setattr(minimap2, "_bam_is_complete", lambda b: b.exists() and b.stat().st_size > 0)

    def test_forced_rerun_regenerates_both_tracks(self, tmp_path, monkeypatch):
        ran: list[str] = []
        self._stub_tools(monkeypatch, tmp_path, ran)
        for name in ("genome.fa", "R1.fq", "R2.fq"):
            (tmp_path / name).write_text("x\n")
        cfg = AssemblyConfig(
            genome=tmp_path / "genome.fa", work_dir=tmp_path, num_cpu=2,
            left_reads=tmp_path / "R1.fq", right_reads=tmp_path / "R2.fq",
        )
        steps = [s for s in STEPS if s.name in ("trinity", "map_transcripts")]

        run_simple_pipeline(ASSEMBLY, steps, cfg)
        ran.clear()
        # --run-trinity (which cascades to map_transcripts) on the completed run
        run_simple_pipeline(
            ASSEMBLY, steps, cfg, force_steps=[step_key(ASSEMBLY, s.name) for s in steps],
        )
        assert ran == [
            "trinity-gg", "trinity-denovo", "trinity-denovo.fasta", "trinity-gg.fasta",
        ]


class TestSkipAppliesToForcedSteps:
    """A step the config disables stays skipped when it is forced: ``-f`` and the
    ``--run-X`` cascades name jaccard and defuse whatever --no-jaccard-clip,
    --defuse and --uniprot say, and a forced defuse without --uniprot crashed."""

    def test_forced_step_matching_skip_does_not_run(self, tmp_path, caplog):
        ran: list[str] = []
        steps = [StepSpec(n, lambda config, n=n: ran.append(n)) for n in ("a", "b", "c")]
        cfg = AssemblyConfig(genome=tmp_path / "g.fa", work_dir=tmp_path, num_cpu=2)
        with caplog.at_level(logging.WARNING, logger="eukan.infra.pipeline"):
            run_simple_pipeline(
                ASSEMBLY, steps, cfg,
                force_steps=[step_key(ASSEMBLY, n) for n in ("a", "b", "c")],
                skip=lambda s: s.name == "b",
            )
        assert ran == ["a", "c"]
        assert "[assembly/b] skipped: disabled by configuration" in caplog.text


class TestPipelineRecords:
    """Each driver keeps its own record of its run in the shared manifest."""

    def test_concurrent_pipelines_keep_each_others_records(self, tmp_path):
        # mask-repeats starts and finishes while assemble is running: assemble's
        # copy of the manifest never sees the repeats record, and saves after it.
        assembly_cfg = AssemblyConfig(
            genome=tmp_path / "asm.fa", work_dir=tmp_path / "assemble", manifest_dir=tmp_path,
            num_cpu=5,
        )
        repeats_cfg = RepeatsConfig(
            genome=tmp_path / "rep.fa", work_dir=tmp_path / "repeats", manifest_dir=tmp_path,
            num_cpu=3,
        )

        def minimap2(config):
            run_simple_pipeline(REPEATS, [StepSpec("modeler", lambda c: None)], repeats_cfg)

        run_simple_pipeline(ASSEMBLY, [StepSpec("minimap2", minimap2)], assembly_cfg)

        status = format_status(load_manifest(tmp_path))
        assert "assembly: completed\n  Started: " in status
        assert f"  Genome: {tmp_path / 'asm.fa'}\n  Genetic code: 1\n  CPUs: 5" in status
        assert "repeats: completed\n  Started: " in status
        assert f"  Genome: {tmp_path / 'rep.fa'}\n  CPUs: 3" in status

    def test_running_while_it_runs_then_failed(self, tmp_path):
        cfg = AssemblyConfig(genome=tmp_path / "g.fa", work_dir=tmp_path, num_cpu=2)
        seen = []

        def step(config):
            seen.append(load_manifest(tmp_path).pipelines[ASSEMBLY].status)
            raise RuntimeError("boom")

        with pytest.raises(RuntimeError):
            run_simple_pipeline(ASSEMBLY, [StepSpec("a", step)], cfg)

        record = load_manifest(tmp_path).pipelines[ASSEMBLY]
        assert seen == ["running"]
        assert record.status == "failed"
        assert record.finished_at is not None
