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
from contextlib import contextmanager
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Any, Iterator

from pydantic import BaseModel, Field

from eukan.infra.logging import get_logger, md5_file

log = get_logger(__name__)

MANIFEST_FILE = "eukan-run.json"
SENTINEL = ".running"


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


# ---------------------------------------------------------------------------
# Step lifecycle -- context manager
# ---------------------------------------------------------------------------


@contextmanager
def pipeline_step(
    work_dir: Path,
    manifest: RunManifest,
    step_name: str,
) -> Iterator[StepRecord]:
    """Context manager for pipeline step lifecycle.

    Usage:
        with pipeline_step(work_dir, manifest, "genemark") as step:
            result = run_genemark(config)
            step.output_file = str(result)

    On __enter__: marks step as running, writes sentinel, saves manifest.
    On __exit__: marks step as completed/failed, removes sentinel, checksums output.
    """
    record = manifest.steps.get(step_name, StepRecord(name=step_name))
    record.status = StepStatus.running
    record.started_at = _now()
    record.error = None
    manifest.steps[step_name] = record

    # Write sentinel
    sdir = work_dir / step_name
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


def is_step_interrupted(work_dir: Path, step_name: str) -> bool:
    """Check if a step was interrupted (sentinel exists)."""
    return (work_dir / step_name / SENTINEL).exists()


def clean_interrupted_step(work_dir: Path, step_name: str) -> None:
    """Remove partial output from an interrupted step."""
    sdir = work_dir / step_name
    if sdir.exists():
        shutil.rmtree(sdir)


# ---------------------------------------------------------------------------
# Tool version snapshot
# ---------------------------------------------------------------------------


_VERSION_TOOLS = [
    ("samtools", ["samtools", "--version"]),
    ("augustus", ["augustus", "--version"]),
    ("star", ["STAR", "--version"]),
    ("hmmer", ["hmmscan", "-h"]),
    ("spaln", ["spaln", "-h"]),
    ("snap", ["snap"]),
    ("trinity", ["Trinity", "--version"]),
]


def _snapshot_tool_versions() -> dict[str, str]:
    """Capture version strings for key tools. Best-effort, never fails."""
    versions: dict[str, str] = {}
    for name, cmd in _VERSION_TOOLS:
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=5)
            output = (result.stdout + result.stderr).strip()
            versions[name] = output.split("\n")[0][:100] if output else "unknown"
        except Exception:
            pass
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
