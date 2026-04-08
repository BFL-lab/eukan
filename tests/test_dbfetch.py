"""Tests for eukan.dbfetch — manifest tracking and integrity checks."""

from pathlib import Path

from eukan.functional.dbfetch import (
    DatabaseEntry,
    _is_current,
    load_db_manifest as load_manifest,
    save_db_manifest as save_manifest,
    check_databases,
)
from eukan.infra.logging import md5_file


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
