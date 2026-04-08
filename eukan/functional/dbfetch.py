"""Download and prepare reference databases with integrity tracking.

Downloads UniProt-SwissProt and Pfam databases, verifies integrity via MD5,
and tracks state in a manifest file so subsequent runs skip up-to-date files.
"""

from __future__ import annotations

import gzip
import json
import shutil
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

import pyhmmer
import requests

from eukan.infra.logging import get_logger, md5_file

log = get_logger(__name__)

UNIPROT_URL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
PFAM_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"

CHUNK_SIZE = 8 * 1024 * 1024  # 8MB
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


def _is_current(entry: DatabaseEntry, db_dir: Path) -> bool:
    """Check if a database file matches its manifest entry."""
    path = db_dir / entry.file
    if not path.exists():
        return False
    return md5_file(path) == entry.md5


# ---------------------------------------------------------------------------
# Download helpers
# ---------------------------------------------------------------------------


def _download_file(url: str, dest: Path) -> None:
    """Download a file with streaming and progress reporting."""
    log.info("Downloading %s...", url)
    with requests.get(url, stream=True, allow_redirects=True, timeout=600) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        downloaded = 0
        with open(dest, "wb") as f:
            for chunk in r.iter_content(chunk_size=CHUNK_SIZE):
                f.write(chunk)
                downloaded += len(chunk)
                if total:
                    pct = downloaded / total * 100
                    log.debug("  %s / %s bytes (%.1f%%)", f"{downloaded:,}", f"{total:,}", pct)


def _make_entry(name: str, filename: str, url: str, db_dir: Path) -> DatabaseEntry:
    """Create a manifest entry for a freshly downloaded/processed file."""
    path = db_dir / filename
    return DatabaseEntry(
        name=name,
        file=filename,
        source_url=url,
        md5=md5_file(path),
        date=datetime.now(timezone.utc).isoformat(),
        size_bytes=path.stat().st_size,
    )


# ---------------------------------------------------------------------------
# Database-specific download + processing
# ---------------------------------------------------------------------------


def download_uniprot(output_dir: Path) -> DatabaseEntry:
    """Download UniProt-SwissProt FASTA and decompress."""
    gz_path = output_dir / "uniprot_sprot.fasta.gz"
    faa_path = output_dir / "uniprot_sprot.faa"

    _download_file(UNIPROT_URL, gz_path)

    log.info("Decompressing UniProt FASTA...")
    with gzip.open(gz_path, "rb") as f_in, open(faa_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

    gz_path.unlink()
    return _make_entry("uniprot", "uniprot_sprot.faa", UNIPROT_URL, output_dir)


def download_pfam(output_dir: Path) -> DatabaseEntry:
    """Download Pfam-A HMM database, decompress, and press it."""
    gz_path = output_dir / "Pfam-A.hmm.gz"
    hmm_path = output_dir / "Pfam-A.hmm"

    _download_file(PFAM_URL, gz_path)

    log.info("Decompressing Pfam HMM...")
    with gzip.open(gz_path, "rb") as f_in, open(hmm_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

    gz_path.unlink()

    log.info("Pressing HMM database...")
    pyhmmer.hmmer.hmmpress(hmms=pyhmmer.plan7.HMMFile(str(hmm_path)), output=str(hmm_path))

    return _make_entry("pfam", "Pfam-A.hmm", PFAM_URL, output_dir)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

# Registry of available databases
DATABASES = {
    "uniprot": ("uniprot_sprot.faa", download_uniprot),
    "pfam": ("Pfam-A.hmm", download_pfam),
}


def fetch_databases(
    output_dir: Path,
    force: bool = False,
    databases: list[str] | None = None,
) -> None:
    """Download and prepare databases, skipping those already up to date.

    Args:
        output_dir: Directory to store databases in.
        force: If True, re-download even if manifest says current.
        databases: List of database names to fetch. None means all.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest = load_db_manifest(output_dir)
    targets = databases or list(DATABASES.keys())

    for db_name in targets:
        if db_name not in DATABASES:
            log.warning("Unknown database: %s", db_name)
            continue

        filename, download_fn = DATABASES[db_name]

        # Check if already current
        if not force and db_name in manifest and _is_current(manifest[db_name], output_dir):
            entry = manifest[db_name]
            log.info("[%s] Up to date (md5:%s..., %s). Use --force to re-download.", db_name, entry.md5[:12], entry.date)
            continue

        log.info("[%s] Fetching...", db_name)
        entry = download_fn(output_dir)
        manifest[db_name] = entry
        save_db_manifest(output_dir, manifest)
        log.info("[%s] Done. md5:%s... size:%s bytes", db_name, entry.md5[:12], f"{entry.size_bytes:,}")


def check_databases(db_dir: Path) -> list[tuple[str, str, bool]]:
    """Check status of all known databases.

    Returns:
        List of (name, status_message, ok) tuples.
    """
    manifest = load_db_manifest(db_dir)
    results: list[tuple[str, str, bool]] = []

    for db_name, (filename, _) in DATABASES.items():
        path = db_dir / filename

        if not path.exists():
            results.append((db_name, f"{filename} not found in {db_dir}", False))
            continue

        if db_name not in manifest:
            size = path.stat().st_size
            results.append((db_name, f"{filename} present ({size:,} bytes) but no manifest entry", False))
            continue

        entry = manifest[db_name]
        if _is_current(entry, db_dir):
            results.append((db_name, f"{filename} OK (md5:{entry.md5[:12]}..., {entry.date})", True))
        else:
            results.append((db_name, f"{filename} CORRUPTED (md5 mismatch)", False))

    return results
