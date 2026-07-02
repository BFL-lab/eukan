"""Tests for GeneMark source-column homogenization.

GeneMark stamps column 2 of its GFF3 as ``GeneMark.hmm3`` (version-dependent),
while the consensus weights map uses the token ``genemark`` — the consensus
engine matches weights by source token, so without normalization the predictions
would be silently weighted zero. ``_genemark_homogenize_source`` rewrites the
source regardless of GeneMark's value.
"""

from __future__ import annotations

import gffutils

from eukan.annotation.genemark import _genemark_homogenize_source
from eukan.gff.normalize import normalize_to_gff3


class TestGenemarkHomogenizeSource:
    """The transform must rewrite source regardless of GeneMark's value."""

    def _feature(self, source: str) -> gffutils.Feature:
        return gffutils.Feature(
            seqid="chr1", source=source, featuretype="gene",
            start=1, end=300, strand="+", attributes={"ID": ["g1"]},
        )

    def test_overrides_genemark_hmm3(self):
        f = self._feature("GeneMark.hmm3")
        assert _genemark_homogenize_source(f).source == "genemark"

    def test_overrides_genemark_hmm(self):
        f = self._feature("GeneMark.hmm")
        assert _genemark_homogenize_source(f).source == "genemark"

    def test_normalize_pipeline_applies_transform(self, tmp_path):
        """Integration: normalize_to_gff3 must wire the post_transform in.

        Mirror how run_genemark calls normalize_to_gff3 — with both the
        source-homogenize post_transform and fix_contig_names — and
        verify column 2 in the output is ``genemark`` end-to-end.
        """
        src = tmp_path / "in.gff"
        src.write_text(
            "##gff-version 3\n"
            "chr1\tGeneMark.hmm3\tgene\t100\t200\t.\t+\t.\tID=g1\n"
            "chr1\tGeneMark.hmm3\tmRNA\t100\t200\t.\t+\t.\tID=g1.t1;Parent=g1\n"
            "chr1\tGeneMark.hmm3\tCDS\t100\t200\t.\t+\t0\tID=cds1;Parent=g1.t1\n"
        )
        out = tmp_path / "out.gff3"

        normalize_to_gff3(
            src, out,
            post_transform=_genemark_homogenize_source,
            fix_contig_names=True,
        )

        for line in out.read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split("\t")
            assert cols[1] == "genemark", f"source not normalized: {line!r}"
