"""Reference databases: what ``eukan db-fetch`` downloads and func-annot reads,
and the manifest that records what was fetched.

The downloaders are in :mod:`eukan.functional.dbfetch` (they need pyhmmer and
requests). This half imports neither, so settings, ``eukan check`` and the CLI
read the database names and file names from here.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Literal, get_args

from eukan.infra.logging import get_logger
from eukan.infra.utils import md5_file

log = get_logger(__name__)


@dataclass(frozen=True)
class Database:
    """A database ``eukan db-fetch`` can fetch."""

    file: str  # what func-annot reads, in the database directory
    url: str  # where db-fetch downloads it from


# ``kofam`` is the pressed HMM database; ``ko_list`` carries the per-KO score
# thresholds and definitions KOfam scoring reads. Both must be present for
# ``func-annot --homology-db kofam``.
DATABASES: dict[str, Database] = {
    "uniprot": Database(
        "uniprot_sprot.faa",
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
    ),
    "pfam": Database(
        "Pfam-A.hmm",
        "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz",
    ),
    "kofam": Database("kofam_eukaryote.hmm", "https://www.genome.jp/ftp/db/kofam/profiles.tar.gz"),
    "ko_list": Database("ko_list.tsv", "https://www.genome.jp/ftp/db/kofam/ko_list.gz"),
}

DEFAULT_DB_DIR = Path("databases")
"""Where ``eukan db-fetch`` puts the databases and func-annot and ``eukan check``
look for them by default, relative to the current directory."""

HomologyDb = Literal["uniprot", "kofam"]
"""A homology database func-annot searches: phmmer against SwissProt, or
hmmscan against KOfam. Pfam hmmscan runs with either."""

HOMOLOGY_DBS: tuple[str, ...] = get_args(HomologyDb)

HOMOLOGY_DB_TARGETS: dict[str, list[str]] = {
    "uniprot": ["uniprot"],
    "kofam": ["kofam", "ko_list"],
}
"""Per-homology-DB set of registry entries to fetch when no explicit
``--database`` list is given. Pfam is always added on top in
:func:`eukan.functional.dbfetch.fetch_databases`.
"""

MANIFEST_FILE = ".manifest.json"


# ---------------------------------------------------------------------------
# Manifest tracking
# ---------------------------------------------------------------------------


@dataclass
class DatabaseEntry:
    """Metadata for a single downloaded database."""

    name: str
    file: str
    source_url: str
    md5: str
    date: str
    size_bytes: int


def load_db_manifest(db_dir: Path) -> dict[str, DatabaseEntry]:
    """Load the database manifest from disk."""
    manifest_path = db_dir / MANIFEST_FILE
    if not manifest_path.exists():
        return {}
    try:
        data = json.loads(manifest_path.read_text())
        return {k: DatabaseEntry(**v) for k, v in data.items()}
    except (json.JSONDecodeError, TypeError, KeyError) as exc:
        log.warning("Corrupt database manifest at %s (%s), re-fetching", manifest_path, exc)
        return {}


def save_db_manifest(db_dir: Path, manifest: dict[str, DatabaseEntry]) -> None:
    """Write the database manifest to disk."""
    manifest_path = db_dir / MANIFEST_FILE
    data = {k: asdict(v) for k, v in manifest.items()}
    manifest_path.write_text(json.dumps(data, indent=2) + "\n")


def make_entry(name: str, db_dir: Path) -> DatabaseEntry:
    """Create a manifest entry for database *name*, just fetched into *db_dir*."""
    db = DATABASES[name]
    path = db_dir / db.file
    return DatabaseEntry(
        name=name,
        file=db.file,
        source_url=db.url,
        md5=md5_file(path),
        date=datetime.now(UTC).isoformat(),
        size_bytes=path.stat().st_size,
    )


def is_current(entry: DatabaseEntry, db_dir: Path) -> bool:
    """Check if a database file matches its manifest entry."""
    path = db_dir / entry.file
    if not path.exists():
        return False
    return md5_file(path) == entry.md5


def check_databases(
    db_dir: Path,
    homology_db: str | None = None,
) -> list[tuple[str, str, bool]]:
    """Check status of databases relevant to the active homology choice.

    When ``homology_db`` is given (``"uniprot"`` or ``"kofam"``), only that
    homology DB's files plus Pfam are checked — the inactive choice is
    skipped so users who never fetch it don't see false-negative warnings.
    When ``None``, checks every entry in :data:`DATABASES`.

    Returns:
        List of (name, status_message, ok) tuples.
    """
    manifest = load_db_manifest(db_dir)
    results: list[tuple[str, str, bool]] = []

    if homology_db is not None:
        wanted = set(HOMOLOGY_DB_TARGETS[homology_db]) | {"pfam"}
    else:
        wanted = set(DATABASES)

    for db_name, db in DATABASES.items():
        if db_name not in wanted:
            continue
        filename = db.file
        path = db_dir / filename

        if not path.exists():
            results.append((db_name, f"{filename} not found in {db_dir}", False))
            continue

        if db_name not in manifest:
            size = path.stat().st_size
            results.append((db_name, f"{filename} present ({size:,} bytes) but no manifest entry", False))
            continue

        entry = manifest[db_name]
        if is_current(entry, db_dir):
            results.append((db_name, f"{filename} OK (md5:{entry.md5[:12]}..., {entry.date})", True))
        else:
            results.append((db_name, f"{filename} CORRUPTED (md5 mismatch)", False))

    return results
