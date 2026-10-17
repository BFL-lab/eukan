"""tests/testdata.py, the test-data downloader behind run_pipeline.py setup-test-data."""

from tests import testdata


def test_accession_files_resolve_outside_the_repo_root(tmp_path, monkeypatch):
    """The tracked accession lists are found whatever the working directory."""
    monkeypatch.chdir(tmp_path)
    for path in (
        testdata.GENOME_ACCESSION_FILE,
        testdata.GENOME_ASSEMBLY_ACCESSION_FILE,
        testdata.NEIGHBOR_ACCESSION_FILE,
        testdata.SRA_ACCESSION_FILE,
    ):
        assert path.is_file(), path
