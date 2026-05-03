"""Tests for eukan.dbfetch — manifest tracking, integrity, and resumable downloads."""

from unittest.mock import patch

from eukan.functional import dbfetch
from eukan.functional.dbfetch import (
    DatabaseEntry,
    _download_file,
    _is_current,
    check_databases,
)
from eukan.functional.dbfetch import (
    load_db_manifest as load_manifest,
)
from eukan.functional.dbfetch import (
    save_db_manifest as save_manifest,
)
from eukan.infra.logging import md5_file


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


class TestManifest:
    def test_round_trip(self, tmp_path):
        """Save and load should preserve entries."""
        entry = DatabaseEntry(
            name="test", file="test.faa", source_url="https://example.com",
            md5="abc123", date="2026-01-01", size_bytes=1000,
        )
        save_manifest(tmp_path, {"test": entry})
        loaded = load_manifest(tmp_path)
        assert "test" in loaded
        assert loaded["test"].md5 == "abc123"
        assert loaded["test"].size_bytes == 1000

    def test_load_missing(self, tmp_path):
        """Missing manifest should return empty dict."""
        assert load_manifest(tmp_path) == {}

    def test_load_corrupt(self, tmp_path):
        """Corrupt manifest should return empty dict."""
        (tmp_path / ".manifest.json").write_text("not json{{{")
        assert load_manifest(tmp_path) == {}


class TestIsCurrent:
    def test_matching_file(self, tmp_path):
        """File with correct md5 should be current."""
        f = tmp_path / "test.faa"
        f.write_text("data")
        md5 = md5_file(f)
        entry = DatabaseEntry("test", "test.faa", "", md5, "", f.stat().st_size)
        assert _is_current(entry, tmp_path)

    def test_modified_file(self, tmp_path):
        """File with wrong md5 should not be current."""
        f = tmp_path / "test.faa"
        f.write_text("data")
        entry = DatabaseEntry("test", "test.faa", "", "wrong_md5", "", 0)
        assert not _is_current(entry, tmp_path)

    def test_missing_file(self, tmp_path):
        """Missing file should not be current."""
        entry = DatabaseEntry("test", "missing.faa", "", "abc", "", 0)
        assert not _is_current(entry, tmp_path)


class TestCheckDatabases:
    def test_missing_databases(self, tmp_path):
        """Empty directory should report all missing."""
        results = check_databases(tmp_path)
        assert len(results) == 2
        assert all(not ok for _, _, ok in results)

    def test_present_without_manifest(self, tmp_path):
        """File present but no manifest should warn."""
        (tmp_path / "uniprot_sprot.faa").write_text(">seq1\nMKK\n")
        results = check_databases(tmp_path)
        uniprot_result = [r for r in results if r[0] == "uniprot"][0]
        assert not uniprot_result[2]  # not OK — no manifest
        assert "no manifest" in uniprot_result[1]

    def test_present_with_valid_manifest(self, tmp_path):
        """File matching manifest should be OK."""
        f = tmp_path / "uniprot_sprot.faa"
        f.write_text(">seq1\nMKK\n")
        entry = DatabaseEntry(
            "uniprot", "uniprot_sprot.faa", "https://example.com",
            md5_file(f), "2026-01-01", f.stat().st_size,
        )
        save_manifest(tmp_path, {"uniprot": entry})
        results = check_databases(tmp_path)
        uniprot_result = [r for r in results if r[0] == "uniprot"][0]
        assert uniprot_result[2]  # OK


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
