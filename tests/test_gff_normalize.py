"""Tests for eukan.gff.normalize.normalize_to_gff3."""

from __future__ import annotations

from pathlib import Path

from eukan.gff.normalize import normalize_to_gff3


def _write_gff(path: Path, content: str) -> None:
    path.write_text(content)


class TestNormalizeToGff3:
    def test_fix_contig_names_strips_fasta_descriptions(self, tmp_path):
        """``fix_contig_names=True`` must invoke the transform, not pass the bool."""
        src = tmp_path / "in.gff"
        # GeneMark-style: seqid is the entire FASTA header up to the first whitespace,
        # but tools occasionally copy the description in too.
        _write_gff(
            src,
            "##gff-version 3\n"
            "chr1 some description\tGeneMark.hmm\tgene\t100\t200\t.\t+\t.\tID=g1\n"
            "chr1 some description\tGeneMark.hmm\tmRNA\t100\t200\t.\t+\t.\tID=g1.t1;Parent=g1\n"
            "chr1 some description\tGeneMark.hmm\tCDS\t100\t200\t.\t+\t0\tID=cds1;Parent=g1.t1\n",
        )
        out = tmp_path / "out.gff3"

        result = normalize_to_gff3(src, out, fix_contig_names=True)

        assert result == out
        text = out.read_text()
        # Description must be stripped from the seqid column (tab-separated col 1).
        for line in text.splitlines():
            if line.startswith("#") or not line.strip():
                continue
            seqid = line.split("\t", 1)[0]
            assert " " not in seqid, f"seqid still has whitespace: {seqid!r}"

    def test_no_transform_passes_through(self, tmp_path):
        """``fix_contig_names=False`` (default) must not mutate seqid."""
        src = tmp_path / "in.gff"
        _write_gff(
            src,
            "##gff-version 3\n"
            "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=g1\n"
            "chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=g1.t1;Parent=g1\n"
            "chr1\tsrc\tCDS\t1\t100\t.\t+\t0\tID=cds1;Parent=g1.t1\n",
        )
        out = tmp_path / "out.gff3"

        normalize_to_gff3(src, out)

        for line in out.read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            assert line.split("\t", 1)[0] == "chr1"
