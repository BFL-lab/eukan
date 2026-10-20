"""Pipeline driver: declarative step specs + the run loop.

A :class:`StepSpec` is a tool-agnostic description of one pipeline step
(name, function to call, output filename, re-run flag display string).
:func:`run_simple_pipeline` consumes a list of them to drive the linear
case (assembly, repeats); :func:`run_orchestrated_step` is the lower-level
primitive used directly by pipelines whose execution graph isn't a
straight line (annotation's phases, functional's JSON-cache dance).
:func:`pipeline_run` opens the manifest for a driver and records the run.
"""

from __future__ import annotations

import shutil
from collections.abc import Callable, Iterator, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from eukan.infra.logging import get_logger
from eukan.infra.manifest import (
    PipelineName,
    PipelineRecord,
    RunManifest,
    StepStatus,
    _now,
    _snapshot_tool_versions,
    get_or_create_manifest,
    save_manifest,
    step_key,
)
from eukan.infra.steps import (
    clean_interrupted_step,
    fingerprint_matches,
    is_step_complete,
    is_step_interrupted,
    pipeline_step,
    run_dir_relative,
)

log = get_logger(__name__)


# ---------------------------------------------------------------------------
# Step specification
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class StepSpec:
    """Declarative description of one pipeline step.

    Attributes:
        name: Bare step name (e.g. ``"minimap2"``). Combined with the
            pipeline prefix to form the manifest key.
        fn: Step function. Called as ``fn(config)``. The function may
            return a Path (used as the step's output) or write to its
            own files and return None.
        output: Filename under ``work_dir`` that the step is known to
            produce. Used for the manifest record's ``output_file`` and
            for stale-output validation. ``None`` if no canonical output;
            such a step resumes on its input fingerprint alone (see ``inputs``).
        flag: Display string for the CLI re-run flag (e.g.
            ``"-A / --run-minimap2"``). Descriptive only: no code reads it,
            since a stale output is rebuilt rather than reported.
            ``None`` if the step has no dedicated re-run flag.
        inputs: Optional ``config -> [Path, ...]`` resolver naming the
            files this step consumes from upstream steps. Their fingerprint
            is recorded on completion; on resume a change re-runs the step
            (instead of reusing an output built from inputs that changed).
            ``None`` if the step has no inputs worth tracking. A step with
            no ``output`` counts complete only while this fingerprint
            matches, so it also names the files it writes: the fingerprint
            is taken when the step completes, so a removed or rewritten
            output then re-runs the step.
        scalars: Optional ``config -> [str, ...]`` resolver naming non-file
            inputs (e.g. ``"max_intron_len=2000"``) folded into the same
            fingerprint, so a change in a config value the step's behaviour
            depends on re-runs it on resume even when its input files are
            byte-identical. ``None`` if the step depends on no such values.
        clears: Other files under ``work_dir`` the step writes beside ``output``,
            deleted with it when a completed step re-runs (forced, or its inputs
            changed). For a step function that skips a sub-output which already
            exists, so an interrupted run resumes without redoing it: without
            this a re-run would regenerate only ``output`` and silently keep the
            other sub-outputs from the previous run.
    """

    name: str
    fn: Callable[..., Any]
    output: str | None = None
    flag: str | None = None
    inputs: Callable[[Any], list[Path]] | None = None
    scalars: Callable[[Any], list[str]] | None = None
    clears: tuple[str, ...] = ()


# ---------------------------------------------------------------------------
# CLI flag → forced-step keys
# ---------------------------------------------------------------------------


def force_steps_from_run_flags(
    pipeline: PipelineName,
    steps: Sequence[StepSpec],
    *,
    force: bool = False,
    downstream: Mapping[str, Sequence[str]] | None = None,
    **run_flags: bool,
) -> list[str]:
    """Translate ``--run-X`` / ``--force`` flags into manifest keys to force.

    For each step ``s``, looks up ``run_flags[f"run_{s.name}"]``. Steps
    whose flag is set are forced, together with every step that consumes
    their output: *downstream* maps a step name to the names of the steps
    that read its output, followed transitively. ``force=True`` (and no
    individual flags) forces every step. Keys come back in pipeline order.
    """
    selected = {s.name for s in steps if run_flags.get(f"run_{s.name}", False)}
    stack = list(selected)
    while stack:
        for dep in (downstream or {}).get(stack.pop(), ()):
            if dep not in selected:
                selected.add(dep)
                stack.append(dep)
    if selected:
        return [step_key(pipeline, s.name) for s in steps if s.name in selected]
    if force:
        return [step_key(pipeline, s.name) for s in steps]
    return []


# ---------------------------------------------------------------------------
# Pipeline run record
# ---------------------------------------------------------------------------


# The inputs each pipeline's record holds: those its config takes from the user.
_RECORDED_INPUTS: dict[PipelineName, Callable[[Any], dict[str, Any]]] = {
    PipelineName.ANNOTATION: lambda c: {
        "genome": c.genome, "proteins": c.proteins, "kingdom": c.kingdom,
        "genetic_code": c.genetic_code, "has_transcripts": c.has_transcripts,
    },
    PipelineName.ASSEMBLY: lambda c: {"genome": c.genome, "genetic_code": c.genetic_code},
    PipelineName.FUNCTIONAL: lambda c: {"proteins": [c.proteins]},
    PipelineName.REPEATS: lambda c: {"genome": c.genome},
}


@contextmanager
def pipeline_run(pipeline: PipelineName, config: Any) -> Iterator[RunManifest]:
    """Yield the run dir's manifest with *pipeline*'s record of this run in it.

    The record -- running, started now, with this run's CPUs, the installed
    tool versions and the pipeline's inputs -- replaces the one from the
    pipeline's previous run and is saved before the body runs. Afterwards it
    is marked completed, or failed if the body raised or was interrupted, and
    saved again.
    """
    manifest = get_or_create_manifest(config.manifest_dir)
    record = PipelineRecord(
        started_at=_now(), num_cpu=config.num_cpu, tool_versions=_snapshot_tool_versions(),
        **_RECORDED_INPUTS[pipeline](config),
    )
    manifest.version = "2"
    manifest.pipelines[pipeline] = record
    save_manifest(config.manifest_dir, manifest)
    try:
        yield manifest
        record.status = "completed"
    except BaseException:
        record.status = "failed"
        raise
    finally:
        record.finished_at = _now()
        save_manifest(config.manifest_dir, manifest)


# ---------------------------------------------------------------------------
# Linear pipeline driver
# ---------------------------------------------------------------------------


def run_simple_pipeline(
    pipeline: PipelineName,
    steps: Sequence[StepSpec],
    config: Any,
    *,
    force_steps: list[str] | None = None,
    skip: Callable[[StepSpec], bool] | None = None,
) -> Path | None:
    """Run a linear pipeline: each step takes ``config`` and runs in order.

    Args:
        pipeline: Manifest-key prefix.
        steps: Step specs, in execution order.
        config: Passed as the first argument to each ``StepSpec.fn``.
        force_steps: Manifest keys to re-run from scratch. ``None`` /
            empty means run all pending steps; cached steps are skipped.
            A non-empty list narrows the active step set to just those
            keys *and* forces re-execution.
        skip: Optional predicate; steps matching it are dropped from
            execution, forced or not (a forced step it drops is logged).

    Returns:
        The output path of the last executed step (or ``None`` if no
        step produced one).
    """
    forced = set(force_steps or ())

    # A missing/corrupt output of a completed step is detected per-step by
    # run_orchestrated_step -> is_step_complete (which returns None and the
    # step rebuilds), so no up-front validation pass is needed here.
    active: list[StepSpec] = []
    for s in steps:
        key = step_key(pipeline, s.name)
        if forced and key not in forced:
            continue
        if skip and skip(s):
            if key in forced:
                log.warning("[%s] skipped: disabled by configuration", key)
            continue
        active.append(s)

    result: Path | None = None
    with pipeline_run(pipeline, config) as manifest:
        for spec in active:
            key = step_key(pipeline, spec.name)
            output_file = (config.work_dir / spec.output) if spec.output else None
            input_files = spec.inputs(config) if spec.inputs else None
            input_scalars = spec.scalars(config) if spec.scalars else None
            result = run_orchestrated_step(
                config.manifest_dir, manifest, key,
                spec.fn, config,
                step_dir=config.work_dir / spec.name,
                force=key in forced,
                output_file=output_file,
                clear_files=[config.work_dir / f for f in spec.clears],
                input_files=input_files,
                input_scalars=input_scalars,
            )
    return result


# ---------------------------------------------------------------------------
# Lower-level primitive: one step with full lifecycle
# ---------------------------------------------------------------------------


def run_orchestrated_step(
    manifest_dir: Path,
    manifest: RunManifest,
    step_name: str,
    fn: Callable[..., Any],
    *args: Any,
    step_dir: Path | None = None,
    force: bool = False,
    output_file: Path | None = None,
    clear_files: list[Path] | None = None,
    input_files: list[Path] | None = None,
    input_scalars: list[str] | None = None,
    **kwargs: Any,
) -> Path | None:
    """Run a pipeline step with full manifest lifecycle management.

    Handles the complete skip-if-complete → clean-interrupted → execute
    dance uniformly across all pipelines.

    Args:
        manifest_dir: Directory containing ``eukan-run.json``.
        manifest: The shared manifest to update.
        step_name: Full manifest key (already prefixed, e.g. ``annotation/genemark``).
        fn: Callable to execute — the caller is responsible for currying
            any config/state needed by the step.
        step_dir: Directory holding the step's sentinel and outputs.
            Defaults to ``manifest_dir / <last segment of step_name>``.
        force: If True, skip the cached-step check and re-run from scratch:
            ``output_file``, ``clear_files`` and the step dir are removed
            first, whatever the step's record says. A completed step whose
            inputs changed re-runs from scratch the same way.
        output_file: If provided, overrides the return value of *fn* as
            the step's output path. Useful for steps that write to a
            fixed filename and return ``None``.
        clear_files: Other outputs deleted along with ``output_file`` when a
            completed step re-runs (see ``StepSpec.clears``).
        input_files: Declared inputs to fingerprint. A completed step whose
            inputs changed since it ran is treated as stale and re-run.
        input_scalars: Non-file inputs (config values) folded into the same
            fingerprint, so a change in one re-runs the step on resume.

    Returns:
        The step's output ``Path`` (cached result, ``output_file``, or
        ``fn``'s return value), or ``None`` if the step has no output.
    """
    prior = manifest.steps.get(step_name)
    if not force:
        cached = is_step_complete(
            manifest_dir, manifest, step_name,
            input_files=input_files, input_scalars=input_scalars,
        )
        if cached:
            log.info("[%s] Already complete, skipping.", step_name)
            return cached if isinstance(cached, Path) else None

    # Re-running a step that previously completed (its inputs changed, or its output
    # is missing or stale): drop the declared output (and clear_files) first, so step
    # functions that reuse an existing output as their own resume check (Trinity's
    # per-mode FASTAs, the transcript mappers' BAM check) actually regenerate it
    # instead of reusing a file built from inputs that have changed. A forced step
    # (whatever its record says) or one whose inputs changed re-runs from scratch, so
    # its step dir goes too: a cache a step function keeps there (GeneMark's
    # genemark.gtf, the AUGUSTUS training set) cannot short-circuit it. A step whose
    # output went missing keeps the step dir: its caches still match its inputs.
    sdir = step_dir if step_dir else manifest_dir / step_name.rsplit("/", 1)[-1]
    reason = None
    if force:
        reason = "forced"
    elif prior is not None and prior.status == StepStatus.completed:
        if (
            (input_files or input_scalars)
            and prior.input_md5 is not None
            and not fingerprint_matches(
                prior.input_md5, input_files, input_scalars, manifest_dir,
            )
        ):
            reason = "inputs changed"
        elif output_file is not None:
            reason = "output missing or stale"
    if reason is not None:
        stale = [
            p for p in (output_file, *(clear_files or ())) if p is not None and p.exists()
        ]
        if reason != "output missing or stale" and sdir != manifest_dir and sdir.exists():
            stale.append(sdir)
        if stale:
            log.warning(
                "[%s] Re-running (%s); clearing stale output %s.",
                step_name, reason, ", ".join(p.name for p in stale),
            )
            for p in stale:
                if p.is_dir():
                    shutil.rmtree(p)
                else:
                    p.unlink()
        else:
            log.info("[%s] Re-running (%s).", step_name, reason)

    if is_step_interrupted(manifest_dir, step_name, step_dir=sdir):
        log.warning("[%s] Cleaning up interrupted previous run...", step_name)
        clean_interrupted_step(manifest_dir, step_name, step_dir=sdir)

    with pipeline_step(
        manifest_dir, manifest, step_name, step_dir=sdir,
        input_files=input_files, input_scalars=input_scalars,
    ) as step:
        result = fn(*args, **kwargs)

        if output_file is not None:
            # Record the declared output even when it's missing: pipeline_step
            # only checksums a path that exists, and is_step_complete flags a
            # completed step whose recorded output is missing/empty on the next
            # run (rebuilding it) — instead of silently recording "no output"
            # and returning the path as if the step had succeeded.
            step.output_file = str(run_dir_relative(output_file, manifest_dir))
            return output_file
        if isinstance(result, (str, Path)):
            step.output_file = str(run_dir_relative(Path(result), manifest_dir))
            return Path(result)
        return None
