"""Pipeline step lifecycle: directory layout, sentinel handling, validation.

The lifecycle pieces live here (rather than next to the run manifest)
because they're orchestration concerns: a sentinel file marks an
in-flight step on disk, and ``pipeline_step`` is the context manager
every wrapper uses to record start/finish/error in the manifest.
"""

from __future__ import annotations

import re
import shutil
from collections.abc import Iterator
from contextlib import contextmanager
from datetime import UTC, datetime
from pathlib import Path

from eukan.exceptions import ConfigurationError
from eukan.infra.logging import get_logger
from eukan.infra.manifest import (
    RunManifest,
    StepRecord,
    StepStatus,
    save_manifest,
)
from eukan.infra.utils import md5_file
from eukan.validation import validate_gff

log = get_logger(__name__)

SENTINEL = ".running"

_SAFE_STEP_NAME = re.compile(r"^[a-zA-Z0-9_\-]+$")


# ---------------------------------------------------------------------------
# Step directory
# ---------------------------------------------------------------------------


def _validate_step_name(step_name: str) -> None:
    """Ensure step_name is safe for use as a directory component."""
    if not _SAFE_STEP_NAME.match(step_name):
        raise ConfigurationError(
            f"Invalid step name: {step_name!r} — must be alphanumeric, hyphens, or underscores only"
        )


def step_dir(work_dir: Path, step_name: str) -> Path:
    """Create and return the working directory for a pipeline step."""
    _validate_step_name(step_name)
    d = work_dir / step_name
    d.mkdir(parents=True, exist_ok=True)
    return d


# ---------------------------------------------------------------------------
# Step lifecycle context manager
# ---------------------------------------------------------------------------


@contextmanager
def pipeline_step(
    work_dir: Path,
    manifest: RunManifest,
    step_name: str,
    step_dir: Path | None = None,
) -> Iterator[StepRecord]:
    """Context manager for pipeline step lifecycle.

    Usage:
        with pipeline_step(work_dir, manifest, "annotation/genemark") as step:
            result = run_genemark(config)
            step.output_file = str(result)

    Args:
        work_dir: Directory containing eukan-run.json.
        manifest: The shared manifest to update.
        step_name: Unique step identifier (used as manifest key).
        step_dir: Directory for the .running sentinel. Defaults to
            ``work_dir / step_name``.

    On __enter__: marks step as running, writes sentinel, saves manifest.
    On __exit__: marks step as completed/failed, removes sentinel, checksums output.
    """
    record = manifest.steps.get(step_name, StepRecord(name=step_name))
    record.status = StepStatus.running
    record.started_at = _now()
    record.error = None
    manifest.steps[step_name] = record

    sdir = step_dir if step_dir else work_dir / step_name
    sdir.mkdir(parents=True, exist_ok=True)
    (sdir / SENTINEL).write_text(f"started: {record.started_at}\n")
    save_manifest(work_dir, manifest)

    log.info("[%s] Running...", step_name)
    try:
        yield record

        record.status = StepStatus.completed
        record.finished_at = _now()
        _compute_duration(record)

        if record.output_file:
            output_path = Path(record.output_file)
            if output_path.exists():
                record.output_md5 = md5_file(output_path)

        log.info("[%s] Done (%.1fs)", step_name, record.duration_seconds or 0)

    except Exception as e:
        record.status = StepStatus.failed
        record.finished_at = _now()
        record.error = str(e)
        _compute_duration(record)
        log.error("[%s] Failed: %s", step_name, e)
        raise

    finally:
        (sdir / SENTINEL).unlink(missing_ok=True)
        save_manifest(work_dir, manifest)


# ---------------------------------------------------------------------------
# Step status predicates
# ---------------------------------------------------------------------------


def is_step_complete(manifest: RunManifest, step_name: str) -> Path | None:
    """Check if a step was completed in a previous run.

    Returns the output path if complete and file exists, else None.
    Pure predicate: callers are responsible for any user-visible logging
    that follows from the result.
    """
    record = manifest.steps.get(step_name)
    if not record or record.status != StepStatus.completed:
        return None
    if not record.output_file:
        return None
    path = Path(record.output_file)
    if path.exists():
        return path
    return None


def is_step_interrupted(work_dir: Path, step_name: str, step_dir: Path | None = None) -> bool:
    """Check if a step was interrupted (sentinel exists)."""
    sdir = step_dir if step_dir else work_dir / step_name
    return (sdir / SENTINEL).exists()


def clean_interrupted_step(work_dir: Path, step_name: str, step_dir: Path | None = None) -> None:
    """Remove partial output from an interrupted step."""
    sdir = step_dir if step_dir else work_dir / step_name
    if sdir.exists():
        shutil.rmtree(sdir)


# ---------------------------------------------------------------------------
# Manifest output validation
# ---------------------------------------------------------------------------


def validate_or_raise(
    manifest: RunManifest,
    expected_steps: list[str],
    step_to_flag: dict[str, str] | None = None,
) -> None:
    """Validate completed step outputs; raise StaleManifestError on any error.

    Replaces the pipeline-level boilerplate of "log each error then
    raise SystemExit(1)" — that bypasses the CLI's structured error
    handling. Use this from inside library code; the CLI handler renders
    StaleManifestError uniformly.
    """
    from eukan.exceptions import StaleManifestError

    errors = validate_step_outputs(manifest, expected_steps, step_to_flag)
    if errors:
        raise StaleManifestError(errors)


def validate_step_outputs(
    manifest: RunManifest,
    expected_steps: list[str],
    step_to_flag: dict[str, str] | None = None,
) -> list[str]:
    """Validate that completed steps have valid output files.

    Checks each expected step in the manifest: if marked completed,
    verifies the output file exists and is non-empty. For GFF outputs,
    additionally verifies the file is structurally valid. Returns a list
    of error messages (empty if all OK).

    Args:
        manifest: The run manifest to check.
        expected_steps: Manifest step keys to validate.
        step_to_flag: Optional mapping of step key to CLI flag for
            actionable error messages. Falls back to the raw step key.
    """
    from eukan.exceptions import GFFValidationError

    errors: list[str] = []
    flag_map = step_to_flag or {}
    for key in expected_steps:
        record = manifest.steps.get(key)
        if not record or record.status != StepStatus.completed:
            continue
        if not record.output_file:
            continue
        output = Path(record.output_file)
        flag = flag_map.get(key, f"(step: {key})")

        if not output.exists() or output.stat().st_size == 0:
            state = "empty" if output.exists() else "missing"
            errors.append(
                f"Step '{key}' is marked complete but output is "
                f"{state}: {output}. Re-run with: {flag}"
            )
            continue

        if output.suffix in (".gff", ".gff3"):
            try:
                validate_gff(output)
            except GFFValidationError as exc:
                errors.append(
                    f"Step '{key}' is marked complete but output is "
                    f"unparseable: {exc}. Re-run with: {flag}"
                )

    return errors


# ---------------------------------------------------------------------------
# Time helpers (shared with manifest)
# ---------------------------------------------------------------------------


def _now() -> str:
    return datetime.now(UTC).isoformat()


def _compute_duration(record: StepRecord) -> None:
    if record.started_at and record.finished_at:
        t0 = datetime.fromisoformat(record.started_at)
        t1 = datetime.fromisoformat(record.finished_at)
        record.duration_seconds = round((t1 - t0).total_seconds(), 1)
