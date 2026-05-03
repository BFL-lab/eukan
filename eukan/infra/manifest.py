"""Run manifest for pipeline state tracking and reproducibility.

Writes `eukan-run.json` in the working directory, recording:
- Pipeline inputs and configuration
- Tool versions discovered at start
- Per-step status, timing, and output checksums
- Overall run state (running/completed/failed)

This enables:
- Resume: re-run skips completed steps
- Reproducibility: exact record of what ran and with what versions
- Diagnostics: `eukan status` reads the manifest to report progress
"""

from __future__ import annotations

import shutil
import subprocess
import threading
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field

from eukan.infra.logging import get_logger, md5_file, validate_gff

log = get_logger(__name__)

MANIFEST_FILE = "eukan-run.json"
SENTINEL = ".running"

# Manifest step keys are prefixed with a pipeline name so the three
# pipelines can share a single eukan-run.json without colliding.
ANNOTATION = "annotation"
ASSEMBLY = "assembly"
FUNCTIONAL = "functional"


def step_key(pipeline: str, name: str) -> str:
    """Build a prefixed manifest key, e.g. ``annotation/genemark``."""
    return f"{pipeline}/{name}"


# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------


class StepStatus(str, Enum):
    pending = "pending"
    running = "running"
    completed = "completed"
    failed = "failed"
    skipped = "skipped"


class StepRecord(BaseModel):
    """Record of a single pipeline step execution."""

    name: str
    status: StepStatus = StepStatus.pending
    started_at: str | None = None
    finished_at: str | None = None
    duration_seconds: float | None = None
    output_file: str | None = None
    output_md5: str | None = None
    error: str | None = None


class RunManifest(BaseModel):
    """Full pipeline run manifest -- serialized to eukan-run.json."""

    version: str = "1"
    status: str = "running"
    started_at: str = ""
    finished_at: str | None = None
    genome: str = ""
    proteins: list[str] = Field(default_factory=list)
    kingdom: str | None = None
    genetic_code: str = ""
    num_cpu: int = 1
    has_transcripts: bool = False
    tool_versions: dict[str, str] = Field(default_factory=dict)
    steps: dict[str, StepRecord] = Field(default_factory=dict)

    def pipeline_status(self, prefix: str) -> str:
        """Get the aggregate status for a pipeline (assembly, annotation, functional).

        Returns 'completed' if all steps with the prefix are completed,
        'running' if any are running, 'failed' if any failed, else 'pending'.
        """
        pipeline_steps = {k: v for k, v in self.steps.items() if k.startswith(prefix + "/")}
        if not pipeline_steps:
            return "pending"
        statuses = {r.status for r in pipeline_steps.values()}
        if StepStatus.failed in statuses:
            return "failed"
        if StepStatus.running in statuses:
            return "running"
        if all(s == StepStatus.completed for s in statuses):
            return "completed"
        return "running"


# ---------------------------------------------------------------------------
# Manifest I/O
# ---------------------------------------------------------------------------


def load_manifest(work_dir: Path) -> RunManifest | None:
    """Load an existing run manifest, or None if not found/corrupt."""
    path = work_dir / MANIFEST_FILE
    if not path.exists():
        return None
    try:
        return RunManifest.model_validate_json(path.read_text())
    except (ValueError, KeyError, OSError) as exc:
        log.warning("Corrupt manifest at %s (%s), ignoring", path, exc)
        return None


_manifest_lock = threading.Lock()


def save_manifest(work_dir: Path, manifest: RunManifest) -> None:
    """Save the run manifest to disk atomically.

    Writes to a temp file first, then renames -- prevents corruption
    if the process crashes mid-write. Thread-safe via lock.
    """
    with _manifest_lock:
        target = work_dir / MANIFEST_FILE
        tmp = target.with_suffix(".tmp")
        tmp.write_text(manifest.model_dump_json(indent=2) + "\n")
        tmp.replace(target)


def init_manifest(config: Any) -> RunManifest:
    """Create a new manifest from pipeline config, snapshotting tool versions."""
    return RunManifest(
        started_at=_now(),
        genome=str(config.genome),
        proteins=[str(p) for p in config.proteins],
        kingdom=config.kingdom.value if config.kingdom else None,
        genetic_code=config.genetic_code,
        num_cpu=config.num_cpu,
        has_transcripts=config.has_transcripts,
        tool_versions=_snapshot_tool_versions(),
    )


def get_or_create_manifest(work_dir: Path, config: Any = None) -> RunManifest:
    """Load existing manifest or create a fresh one.

    This is the preferred entry point for all pipelines.  A single
    eukan-run.json is shared across assembly, annotation, and functional
    pipelines running from the same work directory.
    """
    manifest = load_manifest(work_dir)
    if manifest:
        return manifest
    if config:
        return init_manifest(config)
    return RunManifest(started_at=_now(), tool_versions=_snapshot_tool_versions())


# ---------------------------------------------------------------------------
# Step lifecycle -- context manager
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

    # Write sentinel
    sdir = step_dir if step_dir else work_dir / step_name
    sdir.mkdir(parents=True, exist_ok=True)
    (sdir / SENTINEL).write_text(f"started: {record.started_at}\n")
    save_manifest(work_dir, manifest)

    log.info("[%s] Running...", step_name)
    try:
        yield record

        # Success path
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
    dance uniformly across all three pipelines.

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


def is_step_complete(manifest: RunManifest, step_name: str) -> Path | None:
    """Check if a step was completed in a previous run.

    Returns the output path if complete and file exists, else None.
    """
    record = manifest.steps.get(step_name)
    if not record or record.status != StepStatus.completed:
        return None
    if not record.output_file:
        return None
    path = Path(record.output_file)
    if path.exists():
        log.info("[%s] Already complete, skipping.", step_name)
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
    for step_key in expected_steps:
        record = manifest.steps.get(step_key)
        if not record or record.status != StepStatus.completed:
            continue
        if not record.output_file:
            continue
        output = Path(record.output_file)
        flag = flag_map.get(step_key, f"(step: {step_key})")

        if not output.exists() or output.stat().st_size == 0:
            state = "empty" if output.exists() else "missing"
            errors.append(
                f"Step '{step_key}' is marked complete but output is "
                f"{state}: {output}. Re-run with: {flag}"
            )
            continue

        if output.suffix in (".gff", ".gff3"):
            try:
                validate_gff(output)
            except GFFValidationError as exc:
                errors.append(
                    f"Step '{step_key}' is marked complete but output is "
                    f"unparseable: {exc}. Re-run with: {flag}"
                )

    return errors


# ---------------------------------------------------------------------------
# Tool version snapshot
# ---------------------------------------------------------------------------


_CONDA_TOOLS = ["samtools", "augustus", "star", "hmmer", "spaln", "snap", "trinity"]


def _snapshot_tool_versions() -> dict[str, str]:
    """Capture version strings for key tools via conda. Best-effort, never fails."""
    import os

    prefix = os.environ.get("CONDA_PREFIX", "")
    if not prefix:
        return {}

    # Use conda list for clean version strings
    try:
        result = subprocess.run(
            ["conda", "list", "-p", prefix, "--no-pip"],
            capture_output=True, text=True, timeout=15,
        )
        if result.returncode != 0:
            return {}
    except Exception:
        return {}

    versions: dict[str, str] = {}
    for line in result.stdout.splitlines():
        if line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 2 and parts[0] in _CONDA_TOOLS:
            versions[parts[0]] = parts[1]
    return versions


# ---------------------------------------------------------------------------
# Status formatting
# ---------------------------------------------------------------------------


_STATUS_ICONS = {
    StepStatus.completed: "\u2713",
    StepStatus.running: "\u25b6",
    StepStatus.failed: "\u2717",
    StepStatus.pending: "\u00b7",
    StepStatus.skipped: "-",
}


def format_status(manifest: RunManifest) -> str:
    """Format the manifest as a human-readable status report."""
    lines = [
        f"Run: {manifest.status}",
        f"Started: {manifest.started_at}",
    ]
    if manifest.finished_at:
        lines.append(f"Finished: {manifest.finished_at}")
    lines += [
        f"Genome: {manifest.genome}",
        f"Kingdom: {manifest.kingdom or 'not set'}",
        f"CPUs: {manifest.num_cpu}",
        "",
    ]

    if manifest.tool_versions:
        lines.append("Tool versions:")
        for tool, ver in manifest.tool_versions.items():
            lines.append(f"  {tool:<20s} {ver}")
        lines.append("")

    lines.append("Steps:")
    for step_name, record in manifest.steps.items():
        icon = _STATUS_ICONS.get(record.status, "?")
        duration = f" ({record.duration_seconds}s)" if record.duration_seconds else ""
        lines.append(f"  {icon} {step_name:<25s} {record.status.value}{duration}")
        if record.error:
            lines.append(f"    error: {record.error[:200]}")
        if record.output_md5:
            lines.append(f"    output: {record.output_file} (md5:{record.output_md5[:12]}...)")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _compute_duration(record: StepRecord) -> None:
    if record.started_at and record.finished_at:
        t0 = datetime.fromisoformat(record.started_at)
        t1 = datetime.fromisoformat(record.finished_at)
        record.duration_seconds = round((t1 - t0).total_seconds(), 1)
