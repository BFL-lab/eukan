"""Shared loader for tools.toml.

Single source of truth consumed by :mod:`eukan.check`,
:mod:`eukan.infra.environ`, and ``scripts/generate-env.py``.

Schema (per-tool, all fields optional unless noted)::

    binary          = "augustus"                     # required
    version_cmd     = ["augustus", "--version"]      # default: [binary]
    required_by     = ["annotate"]                   # subcommands needing this tool
    conda_package   = "augustus"                     # omit if not on bioconda
    min_version     = "3.5"                          # conda pin; omit for unpinned
    install_hint    = "..."                          # shown when tool missing

    env_vars = [                                     # env vars the tool needs
        { var = "AUGUSTUS_CONFIG_PATH", path = "config" },
        { var = "PASAHOME", add_to_path = ["scripts"] },
    ]
    conda_path_dirs = ["opt/genemark"]               # prepended to PATH
    conda_lib_dirs  = ["lib"]                        # LD_LIBRARY_PATH (subprocess-only)

``path`` on an env_var entry is resolved under ``$CONDA_PREFIX`` and used
to auto-set the variable under conda. Entries without ``path`` are still
declared as required (checked by ``eukan check``) but must be set by the
environment (e.g., a Dockerfile ENV line).

``add_to_path`` dirs are resolved relative to the env var's value and
prepended to ``PATH`` — used for helper scripts bundled inside a tool's
install tree (e.g., ``$PASAHOME/scripts``).
"""

from __future__ import annotations

import tomllib
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

from eukan.infra.utils import package_resource


@dataclass(frozen=True)
class EnvVarSpec:
    """A single environment variable required by a tool."""

    var: str
    path: str | None = None
    add_to_path: tuple[str, ...] = ()


@dataclass(frozen=True)
class Tool:
    """An external tool dependency loaded from tools.toml."""

    name: str
    binary: str
    version_cmd: tuple[str, ...]
    required_by: tuple[str, ...] = ()
    env_vars: tuple[EnvVarSpec, ...] = ()
    conda_path_dirs: tuple[str, ...] = ()
    conda_lib_dirs: tuple[str, ...] = ()
    conda_package: str | None = None
    min_version: str | None = None
    install_hint: str | None = None

    @property
    def env_var_names(self) -> tuple[str, ...]:
        """Names of every env var this tool declares."""
        return tuple(e.var for e in self.env_vars)


def _find_tools_toml() -> Path | None:
    """Locate tools.toml shipped under ``eukan/data/``."""
    return package_resource("tools.toml")


def _parse_tool(name: str, cfg: dict[str, Any]) -> Tool:
    env_vars = tuple(
        EnvVarSpec(
            var=entry["var"],
            path=entry.get("path"),
            add_to_path=tuple(entry.get("add_to_path", ())),
        )
        for entry in cfg.get("env_vars", ())
    )
    return Tool(
        name=name,
        binary=cfg["binary"],
        version_cmd=tuple(cfg.get("version_cmd", [cfg["binary"]])),
        required_by=tuple(cfg.get("required_by", ())),
        env_vars=env_vars,
        conda_path_dirs=tuple(cfg.get("conda_path_dirs", ())),
        conda_lib_dirs=tuple(cfg.get("conda_lib_dirs", ())),
        conda_package=cfg.get("conda_package"),
        min_version=cfg.get("min_version"),
        install_hint=cfg.get("install_hint"),
    )


@lru_cache(maxsize=1)
def load_tools() -> tuple[Tool, ...]:
    """Load and cache the tool registry from tools.toml.

    Returns an empty tuple if the file cannot be found.
    """
    toml_path = _find_tools_toml()
    if toml_path is None:
        return ()

    with open(toml_path, "rb") as f:
        data = tomllib.load(f)

    return tuple(
        _parse_tool(name, cfg)
        for name, cfg in data.items()
        if isinstance(cfg, dict) and "binary" in cfg
    )
