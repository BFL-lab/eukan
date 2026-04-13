"""Conda environment configuration for external tools.

Reads per-tool conda fields from ``tools.toml`` and provides two
functions that separate *safe* global env-var setup from
*subprocess-only* settings (like ``LD_LIBRARY_PATH``) that would break
system tools if leaked into the parent process.

Per-tool fields in tools.toml::

    conda_env_path  = "share/snap"        # resolves env_var from $CONDA_PREFIX
    conda_env_vars  = [                   # extra env vars (beyond env_var)
        { var = "ALN_TAB", path = "share/spaln/table" },
    ]
    conda_path_dirs = ["opt/genemark"]    # prepended to PATH
    conda_lib_dirs  = ["lib"]             # prepended to LD_LIBRARY_PATH (subprocess-only)

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
from functools import lru_cache
from pathlib import Path
from typing import Any


def _load_tools_toml() -> dict[str, Any]:
    """Load and cache tools.toml."""
    from eukan.check import _find_tools_toml

    toml_path = _find_tools_toml()
    if not toml_path:
        return {}

    from eukan.infra import tomllib

    with open(toml_path, "rb") as f:
        return tomllib.load(f)


@lru_cache(maxsize=1)
def _get_tools() -> dict[str, Any]:
    return _load_tools_toml()


def _resolve(prefix: str, rel_path: str) -> str | None:
    """Resolve a relative path under CONDA_PREFIX.

    Returns the path string if the directory exists, else ``None``.
    """
    resolved = Path(prefix) / rel_path
    return str(resolved) if resolved.is_dir() else None


def _prepend(env: dict[str, str], var: str, value: str) -> None:
    """Prepend *value* to an environment variable if not already present."""
    current = env.get(var, "")
    if value not in current:
        env[var] = f"{value}{os.pathsep}{current}" if current else value


def configure_process_env() -> None:
    """Set global environment variables for the current process.

    Called once at CLI startup.  Reads ``conda_env_path``,
    ``conda_env_vars``, and ``conda_path_dirs`` from each tool in
    tools.toml and applies them.

    Skips ``conda_lib_dirs`` (``LD_LIBRARY_PATH``) — those are
    subprocess-only to avoid breaking system tools (git, curl, etc.).
    """
    prefix = os.environ.get("CONDA_PREFIX")
    if not prefix:
        return

    tools = _get_tools()

    for _name, cfg in tools.items():
        if not isinstance(cfg, dict):
            continue

        # conda_env_path: resolves the tool's env_var
        conda_path = cfg.get("conda_env_path")
        env_var = cfg.get("env_var")
        if conda_path and env_var:
            resolved = _resolve(prefix, conda_path)
            if resolved is not None:
                os.environ.setdefault(env_var, resolved)

        # conda_env_vars: additional env vars beyond env_var
        for entry in cfg.get("conda_env_vars", []):
            resolved = _resolve(prefix, entry["path"])
            if resolved is not None:
                os.environ.setdefault(entry["var"], resolved)

        # conda_path_dirs: prepend to PATH
        for rel_path in cfg.get("conda_path_dirs", []):
            resolved = _resolve(prefix, rel_path)
            if resolved is not None:
                _prepend(os.environ, "PATH", resolved)


def subprocess_env(extra: dict[str, str] | None = None) -> dict[str, str] | None:
    """Build an environment dict for subprocess execution.

    Starts from the current ``os.environ`` (which already has global
    settings applied by :func:`configure_process_env`) and adds
    subprocess-scoped settings — ``LD_LIBRARY_PATH`` entries from
    ``conda_lib_dirs`` on each tool.

    Returns ``None`` when no modifications are needed, letting
    ``subprocess.run`` inherit the parent environment directly.
    """
    prefix = os.environ.get("CONDA_PREFIX", "")

    if not prefix and not extra:
        return None

    env = {**os.environ, **(extra or {})}

    if prefix:
        tools = _get_tools()
        for _name, cfg in tools.items():
            if not isinstance(cfg, dict):
                continue
            for rel_path in cfg.get("conda_lib_dirs", []):
                resolved = _resolve(prefix, rel_path)
                if resolved is not None:
                    _prepend(env, "LD_LIBRARY_PATH", resolved)

    return env
