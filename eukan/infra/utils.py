"""Shared utility functions for the eukan pipeline."""

from __future__ import annotations

import os
from pathlib import Path


def symlink(target: Path, link: Path) -> None:
    """Create a symlink, replacing any existing one."""
    link.unlink(missing_ok=True)
    os.symlink(target, link)
