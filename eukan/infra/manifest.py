"""Run manifest data model and I/O.

Defines the ``RunManifest`` pydantic model serialized to ``eukan-run.json``
and the load/save/init/format helpers around it. Pipeline lifecycle (the
context manager that opens/closes a step record) lives in
:mod:`eukan.infra.steps`; the driver loop that consumes step specs lives
in :mod:`eukan.infra.pipeline`.

The manifest enables:
- Resume: re-run skips completed steps
- Reproducibility: exact record of what ran and with what versions
- Diagnostics: ``eukan status`` reads the manifest to report progress
"""

from __future__ import annotations

import subprocess
import threading
from datetime import UTC, datetime
from enum import Enum, StrEnum
from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field

from eukan.infra.logging import get_logger

log = get_logger(__name__)

MANIFEST_FILE = "eukan-run.json"


class PipelineName(StrEnum):
    """Manifest-key prefix for each pipeline.

    Step keys are prefixed (``annotation/genemark`` etc.) so all pipelines
    can share a single ``eukan-run.json`` without colliding. Inherits
    ``StrEnum`` so the value is used in str-context (f-strings, comparisons).
    """

    ANNOTATION = "annotation"
    ASSEMBLY = "assembly"
    FUNCTIONAL = "functional"
    REPEATS = "repeats"


# Re-exports for callers that prefer the bare-name form.
ANNOTATION = PipelineName.ANNOTATION
ASSEMBLY = PipelineName.ASSEMBLY
FUNCTIONAL = PipelineName.FUNCTIONAL
REPEATS = PipelineName.REPEATS


def step_key(pipeline: str | PipelineName, name: str) -> str:
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
    if the process crashes mid-write.  Thread-safe within a process via
    ``_manifest_lock``; cross-process safe via an fcntl lock on the
    target path so two simultaneous ``eukan`` invocations against the
    same work_dir don't tear writes.

    The fcntl import is local because Windows ships without it; on
    Windows we fall back to in-process locking only.
    """
    payload = manifest.model_dump_json(indent=2) + "\n"
    target = work_dir / MANIFEST_FILE
    tmp = target.with_suffix(".tmp")

    with _manifest_lock:
        try:
            import fcntl
        except ImportError:
            tmp.write_text(payload)
            tmp.replace(target)
            return

        # Lock a sidecar file (target itself may not exist yet, and we
        # don't want to hold a fd to the file we're about to replace).
        lock_path = target.with_suffix(".lock")
        with open(lock_path, "w") as lock_fh:
            fcntl.flock(lock_fh.fileno(), fcntl.LOCK_EX)
            try:
                tmp.write_text(payload)
                tmp.replace(target)
            finally:
                fcntl.flock(lock_fh.fileno(), fcntl.LOCK_UN)


def init_manifest(config: Any) -> RunManifest:
    """Create a new manifest from a pipeline config, snapshotting tool versions.

    Defensive against partial configs: only ``PipelineConfig`` carries
    every field below. ``AssemblyConfig`` / ``RepeatsConfig`` /
    ``FunctionalConfig`` populate what they have and leave the rest at
    their defaults — so whichever pipeline runs first in a clean
    work_dir gets a coherent (if sparse) record. Subsequent pipelines
    sharing the manifest are free to add their step records.
    """
    kingdom = getattr(config, "kingdom", None)
    proteins = getattr(config, "proteins", None) or []
    # Normalize scalar paths (FunctionalConfig has a single `proteins: Path`,
    # PipelineConfig has `list[Path]`).
    if isinstance(proteins, (str, Path)):
        proteins = [proteins]
    return RunManifest(
        started_at=_now(),
        genome=str(getattr(config, "genome", "")),
        proteins=[str(p) for p in proteins],
        kingdom=kingdom.value if kingdom else None,
        genetic_code=getattr(config, "genetic_code", "") or "",
        num_cpu=getattr(config, "num_cpu", 1),
        has_transcripts=bool(getattr(config, "has_transcripts", False)),
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
    if config is not None:
        return init_manifest(config)
    return RunManifest(started_at=_now(), tool_versions=_snapshot_tool_versions())


# ---------------------------------------------------------------------------
# Tool version snapshot
# ---------------------------------------------------------------------------


_CONDA_TOOLS = ["samtools", "augustus", "star", "hmmer", "spaln", "snap", "trinity"]


def _tool_versions_cache_path() -> Path:
    """Per-user cache file for the conda tool-version snapshot."""
    import os

    base = os.environ.get("XDG_CACHE_HOME") or str(Path.home() / ".cache")
    return Path(base) / "eukan" / "tool-versions.json"


def _snapshot_tool_versions() -> dict[str, str]:
    """Capture version strings for key tools via conda. Best-effort, never fails.

    Cached in ``$XDG_CACHE_HOME/eukan/tool-versions.json`` (or
    ``~/.cache/eukan/tool-versions.json``) keyed by ``CONDA_PREFIX`` and
    its mtime — re-running ``conda list`` on every fresh manifest costs
    hundreds of milliseconds, and the result only changes when packages
    are installed or upgraded.
    """
    import json
    import os

    prefix = os.environ.get("CONDA_PREFIX", "")
    if not prefix:
        return {}

    try:
        prefix_mtime = Path(prefix).stat().st_mtime
    except OSError:
        prefix_mtime = 0.0

    cache_path = _tool_versions_cache_path()
    if cache_path.exists():
        try:
            cached = json.loads(cache_path.read_text())
            if cached.get("prefix") == prefix and cached.get("mtime") == prefix_mtime:
                return dict(cached.get("versions") or {})
        except (OSError, ValueError):
            pass

    versions = _query_conda_versions(prefix)

    try:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps(
            {"prefix": prefix, "mtime": prefix_mtime, "versions": versions},
            indent=2,
        ))
    except OSError:
        pass  # caching is best-effort

    return versions


def _query_conda_versions(prefix: str) -> dict[str, str]:
    """Run ``conda list`` once and parse version strings for known tools."""
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

    # Per-pipeline summary above the step detail.
    pipeline_summaries = [
        (p, manifest.pipeline_status(p.value)) for p in PipelineName
    ]
    if any(status != "pending" for _, status in pipeline_summaries):
        lines.append("Pipelines:")
        for pipeline, status in pipeline_summaries:
            if status == "pending":
                continue
            lines.append(f"  {pipeline.value:<12s} {status}")
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
    return datetime.now(UTC).isoformat()
