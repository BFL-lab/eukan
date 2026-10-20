"""Tests for eukan.annotation.consensus — transcript-ORF reintroduction.

``build_consensus_models`` itself is exercised end-to-end (combinr through the
prettify tail) in test_combinr_consensus.py. This file covers the consensus-only
``_patch_in_unmatched_orfs`` step, which folds transcript ORFs that don't overlap
any consensus gene back into the model set, and the ``--orf-patch`` modes that
decide whether it runs and which ORFs it adds.
"""

from __future__ import annotations

from typing import ClassVar

from eukan.annotation import consensus
from eukan.annotation.consensus import _patch_in_unmatched_orfs, _write_prettified_gff3
from eukan.gff import create_gff_db
from eukan.settings import PipelineConfig


def _gene_gff(
    chrom: str, start: int, end: int, gid: str, *, cds: bool = True, strand: str = "+",
) -> str:
    rows = [
        f"{chrom}\tx\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gid}",
        f"{chrom}\tx\tmRNA\t{start}\t{end}\t.\t{strand}\t.\tID={gid}.t1;Parent={gid}",
    ]
    if cds:
        rows.append(
            f"{chrom}\tx\tCDS\t{start}\t{end}\t.\t{strand}\t0\tID={gid}.cds;Parent={gid}.t1"
        )
    return "##gff-version 3\n" + "\n".join(rows) + "\n"


def _gene_ids(db) -> set[str]:
    return {f.id for f in db.features_of_type("gene")}


# One consensus gene as combinr writes it.
_CONSENSUS_GFF = (
    "##gff-version 3\n"
    "chr1\tcombinr\tgene\t1\t300\t.\t+\t.\tID=g1\n"
    "chr1\tcombinr\tmRNA\t1\t300\t.\t+\t.\tID=g1.t1;Parent=g1;support=consensus\n"
    "chr1\tcombinr\texon\t1\t300\t.\t+\t.\tID=g1.exon1;Parent=g1.t1\n"
    "chr1\tcombinr\tCDS\t1\t300\t.\t+\t0\tID=g1.cds1;Parent=g1.t1\n"
)


def _orf_gff(start: int, end: int, gid: str, *, cds: tuple[int, int] | None = None) -> str:
    """One transcript ORF gene on chr1 as find_orfs writes it: gene and mRNA rows
    derived by gffutils, exon and CDS rows from the assembly. The CDS spans the
    whole transcript unless *cds* gives its span."""
    cds_start, cds_end = cds or (start, end)
    rows = [
        f"chr1\tgffutils_derived\tgene\t{start}\t{end}\t.\t+\t.\tID={gid}",
        f"chr1\tgffutils_derived\tmRNA\t{start}\t{end}\t.\t+\t.\tParent={gid};ID={gid}_mRNA",
        f"chr1\tcombinr-assembly\texon\t{start}\t{end}\t.\t+\t.\tID={gid}:exon:1;Parent={gid}_mRNA",
        f"chr1\tcombinr-assembly\tCDS\t{cds_start}\t{cds_end}\t.\t+\t0"
        f"\tID={gid}:CDS:1;Parent={gid}_mRNA",
    ]
    return "\n".join(rows) + "\n"


class TestPatchInUnmatchedOrfs:
    def _consdb(self, tmp_path):
        path = tmp_path / "consensus_models.gff3"
        path.write_text(_gene_gff("chr1", 1, 300, "g1"))
        return create_gff_db(path)

    def test_no_orf_file_returns_unchanged(self, tmp_path):
        consdb = self._consdb(tmp_path)
        result = _patch_in_unmatched_orfs(consdb, tmp_path / "missing.gff3")
        assert result is consdb

    def test_nonoverlapping_orf_is_added(self, tmp_path):
        consdb = self._consdb(tmp_path)
        orf = tmp_path / "transcript_orfs.gff3"
        orf.write_text(_gene_gff("chr1", 1000, 1300, "orf1"))

        result = _patch_in_unmatched_orfs(consdb, orf)

        assert _gene_ids(result) == {"g1", "orf1"}

    def test_orfs_are_added_in_place(self, tmp_path):
        """The ORFs go into the consensus DB itself, not a rebuilt copy."""
        consdb = self._consdb(tmp_path)
        orf = tmp_path / "transcript_orfs.gff3"
        orf.write_text(_gene_gff("chr1", 1000, 1300, "orf1"))

        result = _patch_in_unmatched_orfs(consdb, orf)

        assert result is consdb
        assert _gene_ids(consdb) == {"g1", "orf1"}
        assert [f.id for f in consdb.children("orf1.t1")] == ["orf1.cds"]

    def test_overlapping_orf_not_added(self, tmp_path):
        consdb = self._consdb(tmp_path)
        orf = tmp_path / "transcript_orfs.gff3"
        orf.write_text(_gene_gff("chr1", 50, 250, "orf1"))  # overlaps g1

        result = _patch_in_unmatched_orfs(consdb, orf)

        assert result is consdb

    def test_antisense_orf_not_added(self, tmp_path):
        """An ORF on the opposite strand of a consensus gene overlaps it too."""
        consdb = self._consdb(tmp_path)
        orf = tmp_path / "transcript_orfs.gff3"
        orf.write_text(_gene_gff("chr1", 50, 250, "orf1", strand="-"))  # antisense to g1

        result = _patch_in_unmatched_orfs(consdb, orf)

        assert result is consdb
        assert _gene_ids(consdb) == {"g1"}

    def test_log_counts_genes_not_features(self, tmp_path, caplog):
        """One ORF gene (gene + mRNA + CDS) is one reintroduced ORF, not three."""
        consdb = self._consdb(tmp_path)
        orf = tmp_path / "transcript_orfs.gff3"
        orf.write_text(_gene_gff("chr1", 1000, 1300, "orf1"))

        with caplog.at_level("INFO", logger="eukan.annotation.consensus"):
            _patch_in_unmatched_orfs(consdb, orf)

        assert "Reintroduced 1 transcript ORFs" in caplog.text

    def test_orf_without_cds_not_added(self, tmp_path):
        consdb = self._consdb(tmp_path)
        orf = tmp_path / "transcript_orfs.gff3"
        orf.write_text(_gene_gff("chr1", 1000, 1300, "orf1", cds=False))

        result = _patch_in_unmatched_orfs(consdb, orf)

        assert result is consdb


class TestPatchLabels:
    """The patched genes stay recognisable: every row the patch adds says
    eukan-orf-patch in column 2 and every mRNA it adds carries
    support=transcript_orf, beside combinr's support=consensus."""

    def _patched(self, tmp_path):
        cons = tmp_path / "consensus_models.gff3"
        cons.write_text(_CONSENSUS_GFF)
        orf = tmp_path / "transcript_orfs.gff3"
        orf.write_text(_orf_gff(1000, 1299, "orf1"))
        return _patch_in_unmatched_orfs(create_gff_db(cons), orf)

    def test_added_rows_are_labelled_and_consensus_rows_are_not(self, tmp_path):
        db = self._patched(tmp_path)

        assert {f.id: f.source for f in db.all_features()} == {
            "g1": "combinr", "g1.t1": "combinr", "g1.exon1": "combinr", "g1.cds1": "combinr",
            "orf1": "eukan-orf-patch", "orf1_mRNA": "eukan-orf-patch",
            "orf1:exon:1": "eukan-orf-patch", "orf1:CDS:1": "eukan-orf-patch",
        }

    def test_support_attribute_on_added_mrnas_only(self, tmp_path):
        db = self._patched(tmp_path)

        assert {f.id: f.attributes.get("support") for f in db.all_features()} == {
            "g1": None, "g1.t1": ["consensus"], "g1.exon1": None, "g1.cds1": None,
            "orf1": None, "orf1_mRNA": ["transcript_orf"], "orf1:exon:1": None,
            "orf1:CDS:1": None,
        }

    def test_labels_reach_final_gff3(self, tmp_path):
        """Every final.gff3 row is eukan or, for a patched gene, eukan-orf-patch."""
        _write_prettified_gff3(self._patched(tmp_path), "TEST", tmp_path / "final.gff3")

        rows = [ln.split("\t") for ln in (tmp_path / "final.gff3").read_text().splitlines()]
        assert [(c[2], c[1]) for c in rows] == [
            ("gene", "eukan"), ("mRNA", "eukan"), ("exon", "eukan"), ("CDS", "eukan"),
            ("gene", "eukan-orf-patch"), ("mRNA", "eukan-orf-patch"),
            ("exon", "eukan-orf-patch"), ("CDS", "eukan-orf-patch"),
        ]
        consensus_mrna, patched_mrna = (c[8] for c in rows if c[2] == "mRNA")
        assert "support=consensus" in consensus_mrna
        assert "support=transcript_orf" in patched_mrna


# The protein alignment (prot.gff3) over the ORFs of TestOrfPatchModes.
_PROT_GFF = (
    "##gff-version 3\n"
    # A CDS block on orf_supported's strand overlaps its CDS.
    "chr1\tprot_align\tgene\t1100\t1200\t.\t+\t.\tID=pg1\n"
    "chr1\tprot_align\tmRNA\t1100\t1200\t.\t+\t.\tID=pm1;Parent=pg1\n"
    "chr1\tprot_align\tCDS\t1100\t1200\t.\t+\t0\tID=pc1;Parent=pm1\n"
    # The only CDS block over orf_antisense is on the other strand.
    "chr1\tprot_align\tgene\t2100\t2200\t.\t-\t.\tID=pg2\n"
    "chr1\tprot_align\tmRNA\t2100\t2200\t.\t-\t.\tID=pm2;Parent=pg2\n"
    "chr1\tprot_align\tCDS\t2100\t2200\t.\t-\t0\tID=pc2;Parent=pm2\n"
    # orf_in_gap lies inside a spliced alignment, between two of its CDS blocks.
    "chr1\tprot_align\tgene\t3000\t4000\t.\t+\t.\tID=pg3\n"
    "chr1\tprot_align\tmRNA\t3000\t4000\t.\t+\t.\tID=pm3;Parent=pg3\n"
    "chr1\tprot_align\tCDS\t3000\t3100\t.\t+\t0\tID=pc3a;Parent=pm3\n"
    "chr1\tprot_align\tCDS\t3900\t4000\t.\t+\t0\tID=pc3b;Parent=pm3\n"
    # A CDS block over orf_utr's 5' UTR, short of its CDS.
    "chr1\tprot_align\tgene\t4900\t5050\t.\t+\t.\tID=pg4\n"
    "chr1\tprot_align\tmRNA\t4900\t5050\t.\t+\t.\tID=pm4;Parent=pg4\n"
    "chr1\tprot_align\tCDS\t4900\t5050\t.\t+\t0\tID=pc4;Parent=pm4\n"
)


class TestOrfPatchModes:
    """build_consensus_models patches transcript ORFs in only under --orf-patch:
    'all' adds every ORF gene clear of the consensus, 'supported' only those whose
    CDS overlaps a CDS block of the protein alignment on the same strand."""

    ORFS: ClassVar[dict[str, tuple[int, int, tuple[int, int] | None]]] = {
        "orf_supported": (1000, 1299, None),
        "orf_antisense": (2000, 2299, None),
        "orf_in_gap": (3400, 3599, None),
        "orf_utr": (5000, 5299, (5100, 5299)),
        "orf_bare": (6000, 6299, None),
    }

    def _patched(self, tmp_path, monkeypatch, orf_patch) -> set[str]:
        """The ORFs final.gff3 holds after build_consensus_models (combinr stubbed)."""
        def fake_combinr(config, sdir, evidence, *, transcripts=None):
            (sdir / "consensus_models.gff3").write_text(_CONSENSUS_GFF)
            return sdir / "consensus_models.gff3"

        monkeypatch.setattr(consensus, "run_combinr_consensus", fake_combinr)
        (tmp_path / "orf_finder").mkdir()
        (tmp_path / "orf_finder" / "transcript_orfs.gff3").write_text("".join(
            _orf_gff(start, end, gid, cds=cds) for gid, (start, end, cds) in self.ORFS.items()
        ))
        prot = tmp_path / "prot.gff3"
        prot.write_text(_PROT_GFF)
        genome = tmp_path / "genome.fa"
        genome.write_text(">chr1\nACGT\n")
        config = PipelineConfig(
            genome=genome, proteins=[genome], work_dir=tmp_path, orf_patch=orf_patch,
        )

        final = consensus.build_consensus_models(config, prot, tmp_path / "augustus.gff3")

        starts = {
            int(c[3]) for c in (ln.split("\t") for ln in final.read_text().splitlines())
            if c[2] == "gene" and c[1] == "eukan-orf-patch"
        }
        return {gid for gid, (start, _end, _cds) in self.ORFS.items() if start in starts}

    def test_off_adds_nothing_and_says_how_to_turn_it_on(self, tmp_path, monkeypatch, caplog):
        with caplog.at_level("INFO", logger="eukan.annotation.consensus"):
            assert self._patched(tmp_path, monkeypatch, None) == set()
        assert (
            "ORF patch off; pass --orf-patch supported|all to reintroduce transcript ORFs"
            in caplog.text
        )

    def test_all_adds_every_orf_clear_of_the_consensus(self, tmp_path, monkeypatch):
        assert self._patched(tmp_path, monkeypatch, "all") == set(self.ORFS)

    def test_supported_adds_only_orfs_on_a_same_strand_protein_cds(self, tmp_path, monkeypatch):
        assert self._patched(tmp_path, monkeypatch, "supported") == {"orf_supported"}


class TestWritePrettifiedGff3:
    def test_phases_are_fixed_without_rebuilding_the_db(self, tmp_path, monkeypatch):
        path = tmp_path / "consensus_models.gff3"
        path.write_text(
            "chr1\tx\tgene\t1\t600\t.\t+\t.\tID=g1\n"
            "chr1\tx\tmRNA\t1\t600\t.\t+\t.\tID=g1.t1;Parent=g1\n"
            "chr1\tx\tCDS\t1\t100\t.\t+\t0\tID=g1.cds1;Parent=g1.t1\n"
            "chr1\tx\tCDS\t301\t600\t.\t+\t0\tID=g1.cds2;Parent=g1.t1\n"
        )
        consdb = create_gff_db(path)
        builds = []
        real_create = consensus.create_gff_db
        monkeypatch.setattr(
            consensus, "create_gff_db",
            lambda *a, **k: builds.append(a) or real_create(*a, **k),
        )

        _write_prettified_gff3(consdb, "TEST", tmp_path / "final.gff3")

        assert builds == []
        cds = [ln.split("\t") for ln in (tmp_path / "final.gff3").read_text().splitlines()
               if ln.split("\t")[2] == "CDS"]
        # 100 nt at phase 0 leaves phase 2 for the second CDS.
        assert [(c[3], c[7]) for c in cds] == [("1", "0"), ("301", "2")]
