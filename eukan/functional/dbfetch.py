"""Download and prepare reference databases with integrity tracking.

Downloads UniProt-SwissProt, Pfam, and KOfam databases, verifies
integrity via MD5, and tracks state in a manifest file so subsequent
runs skip up-to-date files. The registry (names, URLs, file names) and the
manifest helpers are in :mod:`eukan.infra.databases`.
"""

from __future__ import annotations

import gzip
import shutil
import tarfile
from collections.abc import Callable
from pathlib import Path

import pyhmmer
import requests

from eukan.infra.databases import (
    DATABASES,
    HOMOLOGY_DB_TARGETS,
    DatabaseEntry,
    is_current,
    load_db_manifest,
    make_entry,
    save_db_manifest,
)
from eukan.infra.logging import get_logger
from eukan.infra.utils import concat_files

log = get_logger(__name__)

CHUNK_SIZE = 8 * 1024 * 1024  # 8MB


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


def _fetch_gunzip(url: str, gz_path: Path, out_path: Path) -> None:
    """Download *url* to *gz_path*, decompress it to *out_path*, delete the .gz."""
    _download_file(url, gz_path)
    log.info("Decompressing %s...", gz_path.name)
    with gzip.open(gz_path, "rb") as f_in, open(out_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    gz_path.unlink()


# ---------------------------------------------------------------------------
# Database-specific download + processing
# ---------------------------------------------------------------------------


def download_uniprot(output_dir: Path) -> DatabaseEntry:
    """Download UniProt-SwissProt FASTA and decompress."""
    db = DATABASES["uniprot"]
    _fetch_gunzip(db.url, output_dir / "uniprot_sprot.fasta.gz", output_dir / db.file)
    return make_entry("uniprot", output_dir)


def download_pfam(output_dir: Path) -> DatabaseEntry:
    """Download Pfam-A HMM database, decompress, and press it."""
    db = DATABASES["pfam"]
    hmm_path = output_dir / db.file
    _fetch_gunzip(db.url, output_dir / "Pfam-A.hmm.gz", hmm_path)

    log.info("Pressing HMM database...")
    # HMMFile would read the .h3m an earlier fetch pressed while hmmpress truncates it.
    for stale in output_dir.glob(f"{db.file}.h3?"):
        stale.unlink()
    with pyhmmer.plan7.HMMFile(str(hmm_path)) as hf:
        pyhmmer.hmmer.hmmpress(hmms=hf, output=str(hmm_path))

    return make_entry("pfam", output_dir)


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
    db = DATABASES["kofam"]
    tarball = output_dir / "profiles.tar.gz"
    extract_dir = output_dir / "_kofam_profiles"
    final_hmm = output_dir / db.file

    _download_file(db.url, tarball)

    log.info("Extracting KOfam profiles to %s...", extract_dir)
    if extract_dir.exists():
        shutil.rmtree(extract_dir)
    extract_dir.mkdir(parents=True)
    with tarfile.open(tarball, "r:gz") as tf:
        # The "data" filter refuses members that would land outside
        # extract_dir; it exists from Python 3.11.4 (pyproject allows 3.11.0).
        if hasattr(tarfile, "data_filter"):
            tf.extractall(extract_dir, filter="data")
        else:
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
    # As in download_pfam: first remove the pressed files of an earlier fetch.
    for stale in output_dir.glob(f"{db.file}.h3?"):
        stale.unlink()
    with pyhmmer.plan7.HMMFile(str(final_hmm)) as hf:
        pyhmmer.hmmer.hmmpress(hmms=hf, output=str(final_hmm))

    return make_entry("kofam", output_dir)


def download_ko_list(output_dir: Path) -> DatabaseEntry:
    """Download the KOfam ``ko_list`` TSV (per-KO thresholds + definitions)."""
    db = DATABASES["ko_list"]
    _fetch_gunzip(db.url, output_dir / "ko_list.gz", output_dir / db.file)
    return make_entry("ko_list", output_dir)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

# The downloader of each database in the registry (infra.databases.DATABASES).
_DOWNLOADERS: dict[str, Callable[[Path], DatabaseEntry]] = {
    "uniprot": download_uniprot,
    "pfam": download_pfam,
    "kofam": download_kofam,
    "ko_list": download_ko_list,
}


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
        if db_name not in _DOWNLOADERS:
            log.warning("Unknown database: %s", db_name)
            continue

        download_fn = _DOWNLOADERS[db_name]

        # Check if already current
        if not force and db_name in manifest and is_current(manifest[db_name], output_dir):
            entry = manifest[db_name]
            log.info("[%s] Up to date (md5:%s..., %s). Use --force to re-download.", db_name, entry.md5[:12], entry.date)
            continue

        log.info("[%s] Fetching...", db_name)
        entry = download_fn(output_dir)
        manifest[db_name] = entry
        save_db_manifest(output_dir, manifest)
        log.info("[%s] Done. md5:%s... size:%s bytes", db_name, entry.md5[:12], f"{entry.size_bytes:,}")

