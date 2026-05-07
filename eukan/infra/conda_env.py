"""Generate a conda environment.yml from the tool registry.

Used by ``scripts/generate-env.py`` to keep ``environment.yml`` in sync
with ``tools.toml``. Independent of the runtime health checks in
:mod:`eukan.infra.health`.
"""

from __future__ import annotations

from eukan.infra.tools_registry import load_tools


def generate_environment_yml() -> str:
    """Generate a conda environment.yml from the tool registry.

    Returns the YAML content as a string.
    """
    tools = load_tools()

    # Deduplicate conda packages (some tools share a package, e.g. hmmer)
    seen: set[str] = set()
    conda_deps: list[str] = []
    for tool in tools:
        if not tool.conda_package or tool.conda_package in seen:
            continue
        seen.add(tool.conda_package)
        pin = f">={tool.min_version}" if tool.min_version else ""
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
        "# Tools requiring manual installation:",
        *manual,
        "#",
        "# Tools bundled with parent packages (no separate install):",
        *bundled,
        "",
        "name: eukan",
        "channels:",
        "  - bioconda",
        "  - conda-forge",
        "  - defaults",
        "",
        "dependencies:",
        "  # Python",
        "  - python>=3.10,<3.13",
        "",
        "  # External tools",
        *sorted(conda_deps),
        "",
        "  # Build dependencies (needed to compile fitild from source)",
        "  - gsl",
        "  - liblbfgs",
        "",
        "  # Perl (needed by AUGUSTUS, PASA, EVM, GeneMark helper scripts)",
        "  - perl",
        "  - perl-bioperl",
        "  - perl-dbd-sqlite",
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
