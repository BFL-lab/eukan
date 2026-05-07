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
    get_or_create_manifest,
    save_manifest,
    step_key,
)
from eukan.infra.steps import (
    clean_interrupted_step,
    is_step_complete,
    is_step_interrupted,
    pipeline_step,
    validate_or_raise,
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
    """

    name: str
    fn: Callable[..., Any]
    output: str | None = None
    flag: str | None = None


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
        active = [s for s in steps if not (skip and skip(s))]
        expected = [step_key(pipeline, s.name) for s in active]
        flag_map = {step_key(pipeline, s.name): s.flag for s in steps if s.flag}
        validate_or_raise(manifest, expected, flag_map)

    save_manifest(config.manifest_dir, manifest)

    result: Path | None = None
    for spec in active:
        key = step_key(pipeline, spec.name)
        output_file = (config.work_dir / spec.output) if spec.output else None
        result = run_orchestrated_step(
            config.manifest_dir, manifest, key,
            spec.fn, config,
            step_dir=config.work_dir / spec.name,
            force=key in forced,
            output_file=output_file,
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

    Returns:
        The step's output ``Path`` (cached result, ``output_file``, or
        ``fn``'s return value), or ``None`` if the step has no output.
    """
    if not force:
        cached = is_step_complete(manifest, step_name)
        if cached:
            log.info("[%s] Already complete, skipping.", step_name)
            return cached

    sdir = step_dir if step_dir else manifest_dir / step_name.rsplit("/", 1)[-1]
    if is_step_interrupted(manifest_dir, step_name, step_dir=sdir):
        log.warning("[%s] Cleaning up interrupted previous run...", step_name)
        clean_interrupted_step(manifest_dir, step_name, step_dir=sdir)

    with pipeline_step(manifest_dir, manifest, step_name, step_dir=sdir) as step:
        result = fn(*args, **kwargs)

        if output_file is not None:
            if output_file.exists():
                step.output_file = str(output_file)
            return output_file
        if isinstance(result, (str, Path)):
            step.output_file = str(result)
            return Path(result)
        return None
