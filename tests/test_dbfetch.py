"""Tests for eukan.functional.dbfetch — resumable downloads and database preparation."""

import gzip
import io
import tarfile
from pathlib import Path
from unittest.mock import patch

import pyhmmer
import pytest

from eukan.functional import dbfetch
from eukan.functional.dbfetch import _download_file
from eukan.infra.databases import DATABASES
from eukan.infra.utils import md5_file


class _FakeResponse:
    """Stand-in for requests.Response usable as a context manager."""

    def __init__(self, body: bytes, status_code: int = 200, headers: dict[str, str] | None = None):
        self.body = body
        self.status_code = status_code
        self.headers = headers or {"content-length": str(len(body))}

    def __enter__(self):
        return self

    def __exit__(self, *_):
        return False

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")

    def iter_content(self, chunk_size: int):
        for i in range(0, len(self.body), chunk_size):
            yield self.body[i : i + chunk_size]


class TestMd5:
    def test_consistent_hash(self, tmp_path):
        """Same content should produce same hash."""
        f = tmp_path / "test.txt"
        f.write_text("hello world")
        h1 = md5_file(f)
        h2 = md5_file(f)
        assert h1 == h2
        assert len(h1) == 32

    def test_different_content(self, tmp_path):
        """Different content should produce different hash."""
        f1 = tmp_path / "a.txt"
        f2 = tmp_path / "b.txt"
        f1.write_text("hello")
        f2.write_text("world")
        assert md5_file(f1) != md5_file(f2)


class TestDownloadResume:
    """Cover the Range-header resume logic in _download_file."""

    def test_fresh_download_writes_full_body(self, tmp_path):
        dest = tmp_path / "file.bin"
        body = b"x" * (dbfetch.CHUNK_SIZE * 2 + 17)

        captured: dict = {}

        def fake_get(url, stream, allow_redirects, timeout, headers):
            captured["headers"] = headers
            return _FakeResponse(body, status_code=200)

        with patch.object(dbfetch.requests, "get", side_effect=fake_get):
            _download_file("https://example.com/file.bin", dest)

        assert dest.read_bytes() == body
        assert not dest.with_suffix(dest.suffix + ".partial").exists()
        assert captured["headers"] == {}

    def test_resume_sends_range_header_and_appends(self, tmp_path):
        dest = tmp_path / "file.bin"
        partial = dest.with_suffix(dest.suffix + ".partial")
        already = b"head" * 10
        partial.write_bytes(already)
        rest = b"tail" * 50

        captured: dict = {}

        def fake_get(url, stream, allow_redirects, timeout, headers):
            captured["headers"] = headers
            return _FakeResponse(rest, status_code=206)

        with patch.object(dbfetch.requests, "get", side_effect=fake_get):
            _download_file("https://example.com/file.bin", dest)

        assert captured["headers"]["Range"] == f"bytes={len(already)}-"
        assert dest.read_bytes() == already + rest
        assert not partial.exists()

    def test_server_ignoring_range_restarts_cleanly(self, tmp_path):
        """If the server returns 200 instead of 206, drop partial and restart."""
        dest = tmp_path / "file.bin"
        partial = dest.with_suffix(dest.suffix + ".partial")
        partial.write_bytes(b"stale-bytes-from-prior-attempt")
        body = b"complete-body"

        def fake_get(url, stream, allow_redirects, timeout, headers):
            # Server signals it sent the full body, not a partial range.
            return _FakeResponse(body, status_code=200)

        with patch.object(dbfetch.requests, "get", side_effect=fake_get):
            _download_file("https://example.com/file.bin", dest)

        assert dest.read_bytes() == body  # no "stale-bytes" prefix
        assert not partial.exists()

    def test_atomic_rename_only_after_success(self, tmp_path):
        dest = tmp_path / "file.bin"
        partial = dest.with_suffix(dest.suffix + ".partial")

        class BoomResponse(_FakeResponse):
            def iter_content(self, chunk_size: int):
                yield b"first-chunk"
                raise RuntimeError("network died mid-stream")

        with patch.object(dbfetch.requests, "get", return_value=BoomResponse(b"")):
            try:
                _download_file("https://example.com/file.bin", dest)
            except RuntimeError:
                pass

        assert not dest.exists()
        assert partial.exists() and partial.read_bytes() == b"first-chunk"


# pyhmmer's bundled test HMMs (10 RRE families), as tests/test_functional_search.py uses them.
RREFAM = Path(pyhmmer.__file__).parent / "tests" / "data" / "hmms" / "txt" / "RREFam.hmm"
needs_rrefam = pytest.mark.skipif(
    not RREFAM.exists(), reason="pyhmmer installed without its bundled test HMMs",
)


def _serve_gzipped_rrefam(url, dest):
    with gzip.open(dest, "wb") as fh:
        fh.write(RREFAM.read_bytes())


def _serve_rrefam_profiles_tarball(url, dest):
    """A profiles.tar.gz with one .hmm per RRE family, all listed in eukaryote.hal."""
    with pyhmmer.plan7.HMMFile(str(RREFAM)) as hf:
        hmms = list(hf)
    with tarfile.open(dest, "w:gz") as tf:
        def add(name: str, data: bytes) -> None:
            info = tarfile.TarInfo(f"profiles/{name}")
            info.size = len(data)
            tf.addfile(info, io.BytesIO(data))

        for hmm in hmms:
            buf = io.BytesIO()
            hmm.write(buf)
            add(f"{hmm.name}.hmm", buf.getvalue())
        add("eukaryote.hal", "".join(f"{hmm.name}.hmm\n" for hmm in hmms).encode())


def _pressed_profile_count(hmm: Path) -> int:
    with pyhmmer.plan7.HMMFile(str(hmm)) as hf:
        return len(hf.optimized_profiles())


@needs_rrefam
class TestDownloadPfam:
    def test_decompresses_presses_and_records(self, tmp_path):
        """The fetched .gz is decompressed to Pfam-A.hmm and removed, the HMMs
        are pressed, and the manifest entry hashes the decompressed file."""
        with patch.object(dbfetch, "_download_file", side_effect=_serve_gzipped_rrefam):
            entry = dbfetch.download_pfam(tmp_path)

        hmm = tmp_path / "Pfam-A.hmm"
        assert hmm.read_bytes() == RREFAM.read_bytes()
        assert not (tmp_path / "Pfam-A.hmm.gz").exists()
        assert entry.md5 == md5_file(hmm)
        assert _pressed_profile_count(hmm) == 10

    def test_refetch_presses_over_earlier_pressed_files(self, tmp_path):
        """A second fetch (db-fetch --force) presses again over the .h3* files of
        the first; HMMFile used to read the old .h3m while hmmpress truncated it
        ("failed to read mat[...]"), leaving the pressed files empty."""
        with patch.object(dbfetch, "_download_file", side_effect=_serve_gzipped_rrefam):
            dbfetch.download_pfam(tmp_path)
            dbfetch.download_pfam(tmp_path)

        assert _pressed_profile_count(tmp_path / "Pfam-A.hmm") == 10


@needs_rrefam
class TestDownloadKofam:
    def test_refetch_presses_over_earlier_pressed_files(self, tmp_path):
        """As for Pfam: a second fetch presses again over the first one's .h3* files."""
        with patch.object(dbfetch, "_download_file", side_effect=_serve_rrefam_profiles_tarball):
            dbfetch.download_kofam(tmp_path)
            dbfetch.download_kofam(tmp_path)

        assert _pressed_profile_count(tmp_path / DATABASES["kofam"].file) == 10


class TestKofamExtraction:
    @pytest.mark.skipif(
        not hasattr(tarfile, "data_filter"), reason="extraction filters need Python >= 3.11.4",
    )
    def test_member_escaping_extract_dir_is_refused(self, tmp_path):
        """A profiles.tar.gz member named ../escaped.txt must not be written
        outside the extraction directory."""

        def fake_download(url, dest):
            data = b"escaped\n"
            info = tarfile.TarInfo("../escaped.txt")
            info.size = len(data)
            with tarfile.open(dest, "w:gz") as tf:
                tf.addfile(info, io.BytesIO(data))

        with patch.object(dbfetch, "_download_file", side_effect=fake_download):
            with pytest.raises(tarfile.OutsideDestinationError):
                dbfetch.download_kofam(tmp_path)

        assert not (tmp_path / "escaped.txt").exists()
