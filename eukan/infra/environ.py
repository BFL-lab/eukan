"""Conda environment configuration for external tools.

Reads per-tool env-var and PATH fields from :mod:`tools_registry` and
separates *safe* global env-var setup from *subprocess-only* settings
(like ``LD_LIBRARY_PATH``) that would break system tools if leaked into
the parent process.

Usage::

    # Once at CLI startup (cli.py):
    from eukan.infra.environ import configure_process_env
    configure_process_env()

    # Per-subprocess (runner.py):
    from eukan.infra.environ import subprocess_env
    subprocess.run(cmd, env=subprocess_env())
"""

from __future__ import annotations

import os
from pathlib import Path

from eukan.infra.tools_registry import Tool, load_tools


def _resolve(prefix: str, rel_path: str) -> str | None:
    """Resolve a path relative to *prefix*, returning None if missing."""
    resolved = Path(prefix) / rel_path
    return str(resolved) if resolved.is_dir() else None


def _prepend(env: dict[str, str], var: str, value: str) -> None:
    """Prepend *value* to *var* unless it's already present."""
    current = env.get(var, "")
    if value not in current:
        env[var] = f"{value}{os.pathsep}{current}" if current else value


def _apply_env_vars(tool: Tool, prefix: str, env: dict[str, str]) -> None:
    """Set env vars declared by *tool* and extend PATH for add_to_path dirs."""
    for spec in tool.env_vars:
        if spec.path:
            # Always set the var when a path is declared. Some tools (e.g.
            # spaln's ALN_DBS) point at user-data directories that aren't
            # populated by the conda package — the env var still needs to
            # be set for the tool to work.
            resolved = str(Path(prefix) / spec.path)
            env.setdefault(spec.var, resolved)

        var_val = env.get(spec.var, "")
        for rel in spec.add_to_path:
            if var_val:
                d = Path(var_val) / rel
                if d.is_dir():
                    _prepend(env, "PATH", str(d))


def _apply_path_dirs(tool: Tool, prefix: str, env: dict[str, str]) -> None:
    for rel in tool.conda_path_dirs:
        resolved = _resolve(prefix, rel)
        if resolved is not None:
            _prepend(env, "PATH", resolved)


def configure_process_env() -> None:
    """Set global env vars for the current process.

    Called once at CLI startup. Skips ``conda_lib_dirs`` — those are
    subprocess-only to avoid breaking system tools (git, curl, etc.).
    """
    prefix = os.environ.get("CONDA_PREFIX")
    if not prefix:
        return

    for tool in load_tools():
        _apply_env_vars(tool, prefix, os.environ)  # type: ignore[arg-type]
        _apply_path_dirs(tool, prefix, os.environ)  # type: ignore[arg-type]


def subprocess_env(extra: dict[str, str] | None = None) -> dict[str, str] | None:
    """Build an environment dict for subprocess execution.

    Starts from ``os.environ`` (which already has global settings from
    :func:`configure_process_env`) and layers subprocess-scoped settings
    — ``LD_LIBRARY_PATH`` entries from each tool's ``conda_lib_dirs``.

    Returns ``None`` when nothing needs adjusting so subprocess.run can
    inherit the parent environment directly.
    """
    prefix = os.environ.get("CONDA_PREFIX", "")

    if not prefix and not extra:
        return None

    env = {**os.environ, **(extra or {})}

    if prefix:
        for tool in load_tools():
            for rel in tool.conda_lib_dirs:
                resolved = _resolve(prefix, rel)
                if resolved is not None:
                    _prepend(env, "LD_LIBRARY_PATH", resolved)

    return env
