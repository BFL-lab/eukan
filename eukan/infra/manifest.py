"""Run manifest data model and I/O.

Defines the ``RunManifest`` pydantic model serialized to ``eukan-run.json``
and the load/save/init helpers around it (``eukan status`` renders it, in
:mod:`eukan.cli.status`). Pipeline lifecycle (the
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

from pydantic import BaseModel, Field, PrivateAttr

from eukan.infra.logging import get_logger
from eukan.infra.tools_registry import load_tools

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
    # Relative to the run dir when inside it (see steps.run_dir_relative); older
    # manifests hold absolute paths, which still resolve.
    output_file: str | None = None
    output_md5: str | None = None
    # Fingerprint of the step's declared *inputs* at the time it last completed.
    # Lets resume detect "inputs changed since this step ran" and re-execute,
    # instead of reusing a stale output (see steps.fingerprint_inputs). ``None``
    # for steps that declare no inputs or were recorded before this field existed.
    input_md5: str | None = None
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
    # This copy as last loaded or saved; save_manifest writes only what changed
    # since, so copies saving to one run dir keep each other's changes.
    _synced: dict[str, Any] = PrivateAttr(default_factory=dict)


# ---------------------------------------------------------------------------
# Manifest I/O
# ---------------------------------------------------------------------------


def load_manifest(work_dir: Path) -> RunManifest | None:
    """Load an existing run manifest, or None if not found/corrupt."""
    path = work_dir / MANIFEST_FILE
    if not path.exists():
        return None
    try:
        manifest = RunManifest.model_validate_json(path.read_text())
    except (ValueError, KeyError, OSError) as exc:
        log.warning("Corrupt manifest at %s (%s), ignoring", path, exc)
        return None
    manifest._synced = manifest.model_dump(mode="json")
    return manifest


_manifest_lock = threading.Lock()


def save_manifest(work_dir: Path, manifest: RunManifest) -> None:
    """Save the run manifest to disk atomically, merged with the one there.

    Several ``eukan`` invocations, and annotation's concurrent steps, save to
    one eukan-run.json. Under the lock the file on disk is re-read, and only
    what this copy changed since it was last loaded or saved (step records,
    run-level fields) is written over it, so a writer keeps the step records
    others saved meanwhile. Writes to a temp file first, then renames --
    prevents corruption if the process crashes mid-write. Thread-safe within a
    process via ``_manifest_lock``; cross-process safe via an fcntl lock on the
    target path so two simultaneous ``eukan`` invocations against the same
    work_dir don't tear writes.

    The fcntl import is local because Windows ships without it; on
    Windows we fall back to in-process locking only.
    """
    target = work_dir / MANIFEST_FILE

    with _manifest_lock:
        try:
            import fcntl
        except ImportError:
            _write_merged(work_dir, manifest)
            return

        # Lock a sidecar file (target itself may not exist yet, and we
        # don't want to hold a fd to the file we're about to replace).
        lock_path = target.with_suffix(".lock")
        with open(lock_path, "w") as lock_fh:
            fcntl.flock(lock_fh.fileno(), fcntl.LOCK_EX)
            try:
                _write_merged(work_dir, manifest)
            finally:
                fcntl.flock(lock_fh.fileno(), fcntl.LOCK_UN)


def _write_merged(work_dir: Path, manifest: RunManifest) -> None:
    """Write *manifest*'s changes since it was last synced over the manifest in
    *work_dir* (temp file, then rename). Caller holds the locks."""
    mine = manifest.model_dump(mode="json")
    on_disk = load_manifest(work_dir)
    merged = mine
    if on_disk is not None:
        merged = on_disk.model_dump(mode="json")
        synced = manifest._synced
        for field, value in mine.items():
            if field != "steps" and value != synced.get(field):
                merged[field] = value
        for key, record in mine["steps"].items():
            if record != synced.get("steps", {}).get(key):
                merged["steps"][key] = record
    target = work_dir / MANIFEST_FILE
    tmp = target.with_suffix(".tmp")
    tmp.write_text(RunManifest.model_validate(merged).model_dump_json(indent=2) + "\n")
    tmp.replace(target)
    manifest._synced = mine


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


def _tool_versions_cache_path() -> Path:
    """Per-user cache file for the conda tool-version snapshot."""
    import os

    base = os.environ.get("XDG_CACHE_HOME") or str(Path.home() / ".cache")
    return Path(base) / "eukan" / "tool-versions.json"


def _snapshot_tool_versions() -> dict[str, str]:
    """Capture the version of every conda package in tools.toml. Best-effort, never fails.

    Cached in ``$XDG_CACHE_HOME/eukan/tool-versions.json`` (or
    ``~/.cache/eukan/tool-versions.json``) keyed by ``CONDA_PREFIX``, its
    mtime and the package list — re-running ``conda list`` on every fresh
    manifest costs hundreds of milliseconds, and the result only changes
    when packages are installed or upgraded.
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

    packages = sorted({t.conda_package for t in load_tools() if t.conda_package})
    key = {"prefix": prefix, "mtime": prefix_mtime, "packages": packages}
    cache_path = _tool_versions_cache_path()
    if cache_path.exists():
        try:
            cached = json.loads(cache_path.read_text())
            if all(cached.get(k) == v for k, v in key.items()):
                return dict(cached.get("versions") or {})
        except (OSError, ValueError):
            pass

    versions = _query_conda_versions(prefix, packages)

    try:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps({**key, "versions": versions}, indent=2))
    except OSError:
        pass  # caching is best-effort

    return versions


def _query_conda_versions(prefix: str, packages: list[str]) -> dict[str, str]:
    """Run ``conda list`` once and parse the version of each of *packages*."""
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
        if len(parts) >= 2 and parts[0] in packages:
            versions[parts[0]] = parts[1]
    return versions


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _now() -> str:
    return datetime.now(UTC).isoformat()
