#!/usr/bin/env python3
"""Pin eukan to a published combinr release.

combinr is fetched as a pinned pre-built release binary by both the Docker image
(``docker/Dockerfile``) and the local installer (``scripts/install-extras.sh``),
each verifying a per-platform SHA-256. Those checksums only exist once the combinr
GitHub release has built its assets, so after cutting a new combinr release run::

    python scripts/refresh-combinr-shas.py 0.1.1

This downloads the per-platform release tarballs, computes their SHA-256s, and
writes both the version and the checksums into ``docker/Dockerfile`` and
``scripts/install-extras.sh``. It replaces whatever value currently sits in each
platform's SHA slot (including the ``SENTINEL_SHA256_*`` placeholders left when the
version is bumped ahead of the release), so it is idempotent and reusable for any
future bump. Review the diff and commit.
"""

from __future__ import annotations

import hashlib
import re
import sys
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DOCKERFILE = REPO / "docker" / "Dockerfile"
INSTALL = REPO / "scripts" / "install-extras.sh"
BASE = "https://github.com/BFL-lab/combinr/releases/download"

# install-extras.sh covers all four targets; the Docker image only needs linux musl x86_64.
TARGETS = [
    "x86_64-unknown-linux-musl",
    "aarch64-unknown-linux-musl",
    "x86_64-apple-darwin",
    "aarch64-apple-darwin",
]
DOCKER_TARGET = "x86_64-unknown-linux-musl"


def fetch_sha(version: str, target: str) -> str:
    url = f"{BASE}/v{version}/combinr-{target}.tar.xz"
    print(f"  {target}: {url}", file=sys.stderr)
    digest = hashlib.sha256()
    with urllib.request.urlopen(url) as resp:  # noqa: S310 (trusted GitHub release URL)
        for chunk in iter(lambda: resp.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def update_install(version: str, shas: dict[str, str]) -> int:
    """Set COMBINR_VERSION and each case arm's sha="..."; returns SHAs replaced."""
    current: str | None = None
    replaced = 0
    out: list[str] = []
    for line in INSTALL.read_text().splitlines(keepends=True):
        if re.match(r"^COMBINR_VERSION=", line):
            out.append(f'COMBINR_VERSION="{version}"\n')
            continue
        t = re.search(r'target="([^"]+)"', line)
        if t:
            current = t.group(1)
        s = re.search(r'sha="([^"]*)"', line)
        if s and current in shas:
            line = line.replace(f'sha="{s.group(1)}"', f'sha="{shas[current]}"')
            replaced += 1
            current = None
        out.append(line)
    INSTALL.write_text("".join(out))
    return replaced


def update_dockerfile(version: str, sha: str) -> None:
    text = DOCKERFILE.read_text()
    text = re.sub(r"^ARG COMBINR_VERSION=.*$", f"ARG COMBINR_VERSION={version}", text, flags=re.M)
    text = re.sub(r"^ARG COMBINR_SHA256=.*$", f"ARG COMBINR_SHA256={sha}", text, flags=re.M)
    DOCKERFILE.write_text(text)


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: refresh-combinr-shas.py <version>  (e.g. 0.1.1)", file=sys.stderr)
        return 2
    version = sys.argv[1].lstrip("v")
    print(f"Fetching combinr v{version} release checksums...", file=sys.stderr)
    shas = {t: fetch_sha(version, t) for t in TARGETS}
    replaced = update_install(version, shas)
    update_dockerfile(version, shas[DOCKER_TARGET])
    if replaced != len(TARGETS):
        print(
            f"WARNING: replaced {replaced}/{len(TARGETS)} SHA slots in install-extras.sh "
            "— inspect the file (its layout may have changed).",
            file=sys.stderr,
        )
    print(
        "Updated docker/Dockerfile and scripts/install-extras.sh. Review the diff and commit.",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
