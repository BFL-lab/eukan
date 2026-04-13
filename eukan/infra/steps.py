"""Pipeline step directory management and completion checking."""

from __future__ import annotations

import re
from pathlib import Path

from eukan.exceptions import ConfigurationError
from eukan.infra.logging import get_logger, validate_gff
from eukan.infra.manifest import is_step_interrupted

log = get_logger(__name__)

_SAFE_STEP_NAME = re.compile(r"^[a-zA-Z0-9_\-]+$")


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


def step_complete(work_dir: Path, step_name: str, output: str) -> Path | None:
    """Check if a pipeline step has already completed.

    For GFF3/GFF outputs, validates the file is parseable.
    For other file types, checks existence and non-emptiness.

    Returns:
        The output path if the step is complete and valid, else None.
    """
    if is_step_interrupted(work_dir, step_name):
        return None

    path = work_dir / step_name / output
    if not path.exists():
        return None

    # Only validate GFF structure for GFF files; for other types just
    # check that the file is non-empty.
    if path.suffix in (".gff3", ".gff"):
        if validate_gff(path):
            return path
        return None

    if path.stat().st_size > 0:
        return path
    return None
