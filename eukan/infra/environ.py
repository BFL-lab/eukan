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

import contextlib
import os
from functools import lru_cache
from pathlib import Path

try:
    import resource
except ImportError:  # pragma: no cover - non-Unix platforms
    resource = None  # type: ignore[assignment]

from eukan.infra.tools_registry import Tool, load_tools


def _resolve(prefix: str, rel_path: str) -> str | None:
    """Resolve a path relative to *prefix*, returning None if missing."""
    resolved = Path(prefix) / rel_path
    return str(resolved) if resolved.is_dir() else None


def _prepend(env: dict[str, str], var: str, value: str) -> None:
    """Prepend *value* to *var* unless it's already present.

    Membership is tested per path component (split on ``os.pathsep``), not by
    substring, so ``/foo`` isn't wrongly treated as present when only
    ``/foobar`` is on the path (and vice-versa for a missed duplicate).
    """
    current = env.get(var, "")
    if value not in current.split(os.pathsep):
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


def _raise_open_file_limit() -> None:
    """Raise the soft open-file limit toward the hard cap.

    STAR's BAM sorter opens roughly ``runThreadN * outBAMsortingBinsN`` temp
    files at once (e.g. 40 threads * 50 bins ~= 2000); samtools and segemehl are
    similarly fd-hungry. When a login shell leaves the soft ``RLIMIT_NOFILE`` at
    the common 1024 default, those tools abort mid-run with "could not create
    output file ..." (STAR exits 109). Bump the soft limit to the hard limit
    (capped) so the child tools we spawn inherit the headroom. Best-effort: any
    failure leaves the inherited shell limit untouched.
    """
    if resource is None:
        return
    try:
        soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    except (ValueError, OSError):
        return
    ceiling = 1_048_576  # don't request an absurd value if hard is unlimited
    target = ceiling if hard == resource.RLIM_INFINITY else min(hard, ceiling)
    if soft >= target:
        return
    with contextlib.suppress(ValueError, OSError):
        resource.setrlimit(resource.RLIMIT_NOFILE, (target, hard))


def configure_process_env() -> None:
    """Set global env vars for the current process.

    Called once at CLI startup. Skips ``conda_lib_dirs`` — those are
    subprocess-only to avoid breaking system tools (git, curl, etc.).
    """
    _raise_open_file_limit()

    prefix = os.environ.get("CONDA_PREFIX")
    if not prefix:
        return

    for tool in load_tools():
        _apply_env_vars(tool, prefix, os.environ)  # type: ignore[arg-type]
        _apply_path_dirs(tool, prefix, os.environ)  # type: ignore[arg-type]


@lru_cache(maxsize=4)
def _subprocess_lib_dirs(prefix: str) -> tuple[str, ...]:
    """Resolve and cache the LD_LIBRARY_PATH additions for *prefix*.

    Cached because the tool registry and conda prefix don't change during
    a process lifetime, but ``subprocess_env`` is called for every
    external command — for pipelines that fan out to thousands of AUGUSTUS
    or alignment partitions this work would otherwise repeat constantly.
    """
    if not prefix:
        return ()
    additions: list[str] = []
    for tool in load_tools():
        for rel in tool.conda_lib_dirs:
            resolved = _resolve(prefix, rel)
            if resolved is not None:
                additions.append(resolved)
    return tuple(additions)


def subprocess_env(extra: dict[str, str] | None = None) -> dict[str, str] | None:
    """Build an environment dict for subprocess execution.

    Starts from ``os.environ`` (which already has global settings from
    :func:`configure_process_env`) and layers subprocess-scoped settings
    — ``LD_LIBRARY_PATH`` entries from each tool's ``conda_lib_dirs``.

    Returns ``None`` when nothing needs adjusting so subprocess.run can
    inherit the parent environment directly.
    """
    prefix = os.environ.get("CONDA_PREFIX", "")
    additions = _subprocess_lib_dirs(prefix)

    if not additions and not extra:
        return None

    env = {**os.environ, **(extra or {})}

    for resolved in additions:
        _prepend(env, "LD_LIBRARY_PATH", resolved)

    return env
