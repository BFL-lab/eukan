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


_SPACED = (
    "##gff-version 3\n"
    "chr1 desc\tGeneMark.hmm3\tgene\t100\t200\t.\t+\t.\tID=g1\n"
    "chr1 desc\tGeneMark.hmm3\tmRNA\t100\t200\t.\t+\t.\tID=g1.t1;Parent=g1\n"
    "chr1 desc\tGeneMark.hmm3\tCDS\t100\t200\t.\t+\t0\tID=cds1;Parent=g1.t1\n"
)


class TestNormalizePasses:
    """normalize_to_gff3 builds no DB it does not need."""

    def test_featuredb_input_without_parse_transform_is_not_copied(self, tmp_path, monkeypatch):
        from eukan.gff import create_gff_db, normalize

        src = tmp_path / "in.gff"
        _write_gff(src, _SPACED)
        db = create_gff_db(src)
        built = []
        real_create = normalize.create_gff_db
        monkeypatch.setattr(
            normalize, "create_gff_db",
            lambda *a, **k: built.append(a) or real_create(*a, **k),
        )
        out = tmp_path / "out.gff3"

        normalize_to_gff3(db, out, fix_contig_names=True)

        assert built == []
        assert [ln.split("\t")[2] for ln in out.read_text().splitlines()] == [
            "gene", "mRNA", "exon", "CDS",
        ]

    def test_post_transform_and_fix_contig_names_share_one_pass(self, tmp_path, monkeypatch):
        from eukan.gff import normalize

        passes = []
        real_transform_db = normalize.transform_db

        def counting(db, transform, **kwargs):
            passes.append(transform)
            return real_transform_db(db, transform, **kwargs)

        monkeypatch.setattr(normalize, "transform_db", counting)
        src = tmp_path / "in.gff"
        _write_gff(src, _SPACED)
        out = tmp_path / "out.gff3"

        def set_source(f):
            f.source = "genemark"
            return f

        normalize_to_gff3(src, out, post_transform=set_source, fix_contig_names=True)

        assert len(passes) == 1
        for line in out.read_text().splitlines():
            cols = line.split("\t")
            assert (cols[0], cols[1]) == ("chr1", "genemark")

    def test_filtering_post_transform_with_fix_contig_names(self, tmp_path):
        """A post_transform that drops features still composes with the seqid fix."""
        src = tmp_path / "in.gff"
        _write_gff(src, _SPACED)
        out = tmp_path / "out.gff3"

        normalize_to_gff3(
            src, out,
            post_transform=lambda f: None if f.featuretype == "CDS" else f,
            fix_contig_names=True,
        )

        rows = [ln.split("\t") for ln in out.read_text().splitlines()]
        assert [r[2] for r in rows] == ["gene", "mRNA", "exon"]
        assert {r[0] for r in rows} == {"chr1"}
