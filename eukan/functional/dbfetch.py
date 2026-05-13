"""Download and prepare reference databases with integrity tracking.

Downloads UniProt-SwissProt, Pfam, and KOfam databases, verifies
integrity via MD5, and tracks state in a manifest file so subsequent
runs skip up-to-date files.
"""

from __future__ import annotations

import gzip
import json
import shutil
import tarfile
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path

import pyhmmer
import requests

from eukan.infra.logging import get_logger
from eukan.infra.utils import concat_files, md5_file

log = get_logger(__name__)

UNIPROT_URL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
PFAM_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
KOFAM_PROFILES_URL = "https://www.genome.jp/ftp/db/kofam/profiles.tar.gz"
KOFAM_KO_LIST_URL = "https://www.genome.jp/ftp/db/kofam/ko_list.gz"

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
    """Download a file with resume support.

    Streams to ``dest.partial`` and atomically renames to ``dest`` on
    success.  If ``dest.partial`` already exists from a previous failed
    attempt, sends a Range request to resume from where we stopped.
    Servers that don't support Range fall back to a clean restart.
    """
    log.info("Downloading %s...", url)
    partial = dest.with_suffix(dest.suffix + ".partial")
    resume_from = partial.stat().st_size if partial.exists() else 0

    headers: dict[str, str] = {}
    if resume_from > 0:
        headers["Range"] = f"bytes={resume_from}-"
        log.info("Resuming from byte %s", f"{resume_from:,}")

    with requests.get(
        url, stream=True, allow_redirects=True, timeout=600, headers=headers,
    ) as r:
        r.raise_for_status()

        # Some servers ignore Range and return 200 + the full body.
        # In that case discard any partial data and start over.
        if resume_from > 0 and r.status_code != 206:
            log.warning("Server ignored Range header; restarting download.")
            partial.unlink(missing_ok=True)
            resume_from = 0

        mode = "ab" if resume_from > 0 else "wb"
        total_remaining = int(r.headers.get("content-length", 0))
        total = resume_from + total_remaining
        downloaded = resume_from
        with open(partial, mode) as f:
            for chunk in r.iter_content(chunk_size=CHUNK_SIZE):
                f.write(chunk)
                downloaded += len(chunk)
                if total:
                    pct = downloaded / total * 100
                    log.debug("  %s / %s bytes (%.1f%%)", f"{downloaded:,}", f"{total:,}", pct)

    partial.replace(dest)


def _make_entry(name: str, filename: str, url: str, db_dir: Path) -> DatabaseEntry:
    """Create a manifest entry for a freshly downloaded/processed file."""
    path = db_dir / filename
    return DatabaseEntry(
        name=name,
        file=filename,
        source_url=url,
        md5=md5_file(path),
        date=datetime.now(UTC).isoformat(),
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


KOFAM_HMM_FILE = "kofam_eukaryote.hmm"
KOFAM_KO_LIST_FILE = "ko_list.tsv"


def _read_hal(hal_path: Path) -> list[str]:
    """Read a .hal file and return the list of HMM filenames (skipping comments)."""
    result: list[str] = []
    with open(hal_path) as f:
        for line in f:
            s = line.strip()
            if s and not s.startswith("#"):
                result.append(s)
    return result


def download_kofam(output_dir: Path) -> DatabaseEntry:
    """Download KOfam profiles, filter to eukaryote subset, concatenate, press.

    Uses ``eukaryote.hal`` shipped inside ``profiles.tar.gz`` to pick the
    eukaryote-only KO subset (~16k of ~27k total) — keeps the pressed DB
    roughly proportional to the subset size. The per-KO ``.hmm`` files
    are extracted to a temporary directory and removed after pressing.
    """
    tarball = output_dir / "profiles.tar.gz"
    extract_dir = output_dir / "_kofam_profiles"
    final_hmm = output_dir / KOFAM_HMM_FILE

    _download_file(KOFAM_PROFILES_URL, tarball)

    log.info("Extracting KOfam profiles to %s...", extract_dir)
    if extract_dir.exists():
        shutil.rmtree(extract_dir)
    extract_dir.mkdir(parents=True)
    with tarfile.open(tarball, "r:gz") as tf:
        tf.extractall(extract_dir)
    tarball.unlink()

    # Tarball layout: profiles/eukaryote.hal + profiles/K*.hmm
    profiles_root = extract_dir / "profiles"
    if not profiles_root.is_dir():
        # Some mirrors flatten the tarball.
        profiles_root = extract_dir

    hal = profiles_root / "eukaryote.hal"
    if not hal.is_file():
        raise FileNotFoundError(
            f"eukaryote.hal not found inside KOfam tarball under {profiles_root}"
        )

    hmm_files = _read_hal(hal)
    log.info("Concatenating %d eukaryotic KOfam HMM files...", len(hmm_files))
    sources = [profiles_root / name for name in hmm_files]
    missing = [str(p) for p in sources if not p.exists()]
    if missing:
        raise FileNotFoundError(
            f"{len(missing)} HMM files listed in eukaryote.hal are missing "
            f"(first: {missing[0]})"
        )
    concat_files(sources, final_hmm)

    log.info("Removing extracted per-KO HMM files...")
    shutil.rmtree(extract_dir)

    log.info("Pressing KOfam HMM database (%s)...", final_hmm.name)
    pyhmmer.hmmer.hmmpress(
        hmms=pyhmmer.plan7.HMMFile(str(final_hmm)), output=str(final_hmm),
    )

    return _make_entry("kofam", KOFAM_HMM_FILE, KOFAM_PROFILES_URL, output_dir)


def download_ko_list(output_dir: Path) -> DatabaseEntry:
    """Download the KOfam ``ko_list`` TSV (per-KO thresholds + definitions)."""
    gz_path = output_dir / "ko_list.gz"
    tsv_path = output_dir / KOFAM_KO_LIST_FILE

    _download_file(KOFAM_KO_LIST_URL, gz_path)

    log.info("Decompressing ko_list...")
    with gzip.open(gz_path, "rb") as f_in, open(tsv_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

    gz_path.unlink()
    return _make_entry("ko_list", KOFAM_KO_LIST_FILE, KOFAM_KO_LIST_URL, output_dir)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

# Registry of available databases. ``kofam`` produces a pressed HMM;
# ``ko_list`` carries the per-KO score thresholds + definitions used by
# KOfam scoring. Both must be present for ``func-annot --homology-db kofam``.
DATABASES = {
    "uniprot": ("uniprot_sprot.faa", download_uniprot),
    "pfam": ("Pfam-A.hmm", download_pfam),
    "kofam": (KOFAM_HMM_FILE, download_kofam),
    "ko_list": (KOFAM_KO_LIST_FILE, download_ko_list),
}


HOMOLOGY_DB_TARGETS: dict[str, list[str]] = {
    "uniprot": ["uniprot"],
    "kofam": ["kofam", "ko_list"],
}
"""Per-homology-DB set of registry entries to fetch when no explicit
``--database`` list is given. Pfam is always added on top in
:func:`fetch_databases`.
"""


def fetch_databases(
    output_dir: Path,
    force: bool = False,
    databases: list[str] | None = None,
    homology_db: str = "uniprot",
) -> None:
    """Download and prepare databases, skipping those already up to date.

    Args:
        output_dir: Directory to store databases in.
        force: If True, re-download even if manifest says current.
        databases: Explicit list of registry entries to fetch. When given,
            overrides ``homology_db`` (the caller is asking for an exact
            subset, e.g. ``["pfam"]`` to refresh only Pfam).
        homology_db: Choice of homology DB to fetch alongside Pfam when
            ``databases`` is not given. ``"uniprot"`` (default) keeps the
            current SwissProt-based pipeline; ``"kofam"`` triggers the
            KOfam profiles + ko_list pair.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest = load_db_manifest(output_dir)

    if databases is None:
        if homology_db not in HOMOLOGY_DB_TARGETS:
            raise ValueError(
                f"Unknown homology_db {homology_db!r}; expected one of "
                f"{sorted(HOMOLOGY_DB_TARGETS)}"
            )
        targets = HOMOLOGY_DB_TARGETS[homology_db] + ["pfam"]
    else:
        targets = databases

    for db_name in targets:
        if db_name not in DATABASES:
            log.warning("Unknown database: %s", db_name)
            continue

        _filename, download_fn = DATABASES[db_name]

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

    for db_name, (filename, _) in DATABASES.items():
        if db_name not in wanted:
            continue
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
