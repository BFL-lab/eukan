#!/usr/bin/env python3
"""Regenerate environment.yml from tools.toml.

Usage:
    python scripts/generate-env.py [-o environment.yml]
"""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate environment.yml from tools.toml")
    parser.add_argument(
        "-o", "--output", type=Path, default="environment.yml",
        help="Output file path (default: environment.yml)",
    )
    args = parser.parse_args()

    from eukan.check import generate_environment_yml

    content = generate_environment_yml()
    args.output.write_text(content)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
