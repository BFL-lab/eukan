"""Pipeline driver: declarative step specs + the run loop.

A :class:`StepSpec` is a tool-agnostic description of one pipeline step
(name, function to call, output filename, re-run flag display string).
:func:`run_simple_pipeline` consumes a list of them to drive the linear
case (assembly, repeats); :func:`run_orchestrated_step` is the lower-level
primitive used directly by pipelines whose execution graph isn't a
straight line (annotation's phases, functional's JSON-cache dance).
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from eukan.infra.logging import get_logger
from eukan.infra.manifest import (
    PipelineName,
    RunManifest,
    StepStatus,
    get_or_create_manifest,
    save_manifest,
    step_key,
)
from eukan.infra.steps import (
    clean_interrupted_step,
    fingerprint_inputs,
    is_step_complete,
    is_step_interrupted,
    pipeline_step,
)

log = get_logger(__name__)


# ---------------------------------------------------------------------------
# Step specification
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class StepSpec:
    """Declarative description of one pipeline step.

    Attributes:
        name: Bare step name (e.g. ``"star"``). Combined with the
            pipeline prefix to form the manifest key.
        fn: Step function. Called as ``fn(config)``. The function may
            return a Path (used as the step's output) or write to its
            own files and return None.
        output: Filename under ``work_dir`` that the step is known to
            produce. Used for the manifest record's ``output_file`` and
            for stale-output validation. ``None`` if no canonical output.
        flag: Display string for the CLI re-run flag (e.g.
            ``"-A / --run-star"``). Shown in stale-output error messages.
            ``None`` if the step has no dedicated re-run flag.
        inputs: Optional ``config -> [Path, ...]`` resolver naming the
            files this step consumes from upstream steps. Their fingerprint
            is recorded on completion; on resume a change re-runs the step
            (instead of reusing an output built from inputs that changed).
            ``None`` if the step has no inputs worth tracking.
        scalars: Optional ``config -> [str, ...]`` resolver naming non-file
            inputs (e.g. ``"max_intron_len=2000"``) folded into the same
            fingerprint, so a change in a config value the step's behaviour
            depends on re-runs it on resume even when its input files are
            byte-identical. ``None`` if the step depends on no such values.
    """

    name: str
    fn: Callable[..., Any]
    output: str | None = None
    flag: str | None = None
    inputs: Callable[[Any], list[Path]] | None = None
    scalars: Callable[[Any], list[str]] | None = None


# ---------------------------------------------------------------------------
# CLI flag → forced-step keys
# ---------------------------------------------------------------------------


def force_steps_from_run_flags(
    pipeline: PipelineName,
    steps: Sequence[StepSpec],
    *,
    force: bool = False,
    **run_flags: bool,
) -> list[str]:
    """Translate ``--run-X`` / ``--force`` flags into manifest keys to force.

    For each step ``s``, looks up ``run_flags[f"run_{s.name}"]``. Steps
    whose flag is set are forced individually. ``force=True`` (and no
    individual flags) forces every step.
    """
    selected = [s.name for s in steps if run_flags.get(f"run_{s.name}", False)]
    if selected:
        return [step_key(pipeline, n) for n in selected]
    if force:
        return [step_key(pipeline, s.name) for s in steps]
    return []


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
            execution unless they were explicitly forced.

    Returns:
        The output path of the last executed step (or ``None`` if no
        step produced one).
    """
    manifest = get_or_create_manifest(config.manifest_dir, config)
    forced = set(force_steps or ())

    if forced:
        active = [s for s in steps if step_key(pipeline, s.name) in forced]
    else:
        # A missing/corrupt output of a completed step is detected per-step by
        # run_orchestrated_step -> is_step_complete (which returns None and the
        # step rebuilds), so no up-front validation pass is needed here.
        active = [s for s in steps if not (skip and skip(s))]

    save_manifest(config.manifest_dir, manifest)

    result: Path | None = None
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
        force: If True, skip the cached-step check and re-run.
        output_file: If provided, overrides the return value of *fn* as
            the step's output path. Useful for steps that write to a
            fixed filename and return ``None``.
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
            manifest, step_name, input_files=input_files, input_scalars=input_scalars
        )
        if cached:
            log.info("[%s] Already complete, skipping.", step_name)
            return cached

    # Re-running a step that previously completed (forced, or its inputs changed):
    # drop the declared output first, so step functions that reuse a "complete"
    # output as their own resume check (the transcript mappers' BAM check) actually
    # regenerate it instead of reusing a file built from inputs that have changed.
    if prior is not None and prior.status == StepStatus.completed and output_file is not None:
        if force:
            reason = "forced"
        elif (
            (input_files or input_scalars)
            and prior.input_md5 is not None
            and fingerprint_inputs(input_files, input_scalars) != prior.input_md5
        ):
            reason = "inputs changed"
        else:
            reason = "output missing or stale"
        if output_file.exists():
            log.warning(
                "[%s] Re-running (%s); clearing stale output %s.",
                step_name, reason, output_file.name,
            )
            output_file.unlink()
        else:
            log.info("[%s] Re-running (%s).", step_name, reason)

    sdir = step_dir if step_dir else manifest_dir / step_name.rsplit("/", 1)[-1]
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
            step.output_file = str(output_file)
            return output_file
        if isinstance(result, (str, Path)):
            step.output_file = str(result)
            return Path(result)
        return None
