#!/usr/bin/env python3
"""Regenerate environment.yml from tools.toml.

Usage:
    python scripts/generate-env.py [-o environment.yml]
"""

from __future__ import annotations

import argparse
from pathlib import Path

from eukan.infra.tools_registry import load_tools


def generate_environment_yml() -> str:
    """Generate a conda environment.yml from the tool registry.

    Returns the YAML content as a string.
    """
    tools = load_tools()

    # Deduplicate conda packages (two tools may name the same package)
    seen: set[str] = set()
    conda_deps: list[str] = []
    for tool in tools:
        if not tool.conda_package or tool.conda_package in seen:
            continue
        seen.add(tool.conda_package)
        if tool.conda_pin:
            pin = tool.conda_pin
        elif tool.min_version:
            pin = f">={tool.min_version}"
        else:
            pin = ""
        conda_deps.append(f"  - {tool.conda_package}{pin}")

    # Tools requiring manual install
    manual: list[str] = []
    for tool in tools:
        if not tool.conda_package and tool.install_hint:
            manual.append(f"#   - {tool.name}: {tool.install_hint}")

    # Bundled tools (no conda package, no install hint)
    bundled: list[str] = []
    for tool in tools:
        if not tool.conda_package and not tool.install_hint:
            bundled.append(f"#   - {tool.binary}: bundled with parent package")

    lines = [
        "# Auto-generated from tools.toml — do not edit manually.",
        "# Regenerate with: python scripts/generate-env.py",
        "#",
        "# Create with STRICT channel priority (bioconda requirement — without it the",
        "# solver mixes channel builds and yields inconsistent deps, e.g. r-stringi vs icu):",
        "#   CONDA_CHANNEL_PRIORITY=strict conda env create -f environment.yml",
        "#",
        "# Tools requiring manual installation:",
        *manual,
        "#",
        "# Tools bundled with parent packages (no separate install):",
        *bundled,
        "",
        "name: eukan",
        "channels:",
        # bioconda requires this exact channel order (conda-forge first), paired
        # with strict channel priority; a bioconda-first / flexible-priority solve
        # mixes builds across channels and yields inconsistent runtime deps (e.g.
        # an r-stringi linked against an icu version not in the env).
        "  - conda-forge",
        "  - bioconda",
        "  - defaults",
        "",
        "dependencies:",
        "  # Python",
        "  - python>=3.11,<3.13",
        "",
        "  # External tools",
        *sorted(conda_deps),
        "",
        "  # Build dependencies (needed to compile fitild from source)",
        "  - gsl",
        "  - liblbfgs",
        "",
        "  # Perl (needed by AUGUSTUS and GeneMark helper scripts)",
        "  - perl",
        "  - perl-yaml",
        "  - perl-hash-merge",
        "  - perl-mce",
        "  - perl-parallel-forkmanager",
        "  - perl-math-utils",
        "",
        "  # eukan + Python dependencies (installed via pip for version consistency)",
        "  - pip",
        "  - pip:",
        "    - .",
        "",
    ]

    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate environment.yml from tools.toml")
    parser.add_argument(
        "-o", "--output", type=Path, default="environment.yml",
        help="Output file path (default: environment.yml)",
    )
    args = parser.parse_args()

    content = generate_environment_yml()
    args.output.write_text(content)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
