"""Tests for non-canonical splice site handling and genetic code integration."""

from __future__ import annotations

import json
from pathlib import Path

from eukan.gencode import GeneticCode

# ---------------------------------------------------------------------------
# GeneticCode.genemark_flag
# ---------------------------------------------------------------------------


class TestGeneticCodeGenemarkFlag:
    def test_code_1_no_flag(self):
        """Code 1 (standard) should not produce a --gcode flag."""
        gc = GeneticCode(1)
        assert gc.genemark_flag == []

    def test_code_6_produces_flag(self):
        """Code 6 (Ciliate) is supported by GeneMark."""
        gc = GeneticCode(6)
        assert gc.genemark_flag == ["--gcode=6"]

    def test_code_26_produces_flag(self):
        """Code 26 is supported (remapped to 6 internally by GeneMark)."""
        gc = GeneticCode(26)
        assert gc.genemark_flag == ["--gcode=26"]

    def test_unsupported_code_no_flag(self):
        """Unsupported codes should return empty list, not error."""
        gc = GeneticCode(10)  # Euplotes — not in GeneMark's set
        assert gc.genemark_flag == []
        assert gc.is_genemark_supported is False

    def test_is_genemark_supported(self):
        assert GeneticCode(1).is_genemark_supported is True
        assert GeneticCode(6).is_genemark_supported is True
        assert GeneticCode(12).is_genemark_supported is False


# ---------------------------------------------------------------------------
# _splice_type_to_augustus
# ---------------------------------------------------------------------------


class TestSpliceTypeToAugustus:
    def test_standard_types(self):
        from eukan.annotation.augustus import _splice_type_to_augustus

        assert _splice_type_to_augustus("GT-AG") == "gtag"
        assert _splice_type_to_augustus("GC-AG") == "gcag"
        assert _splice_type_to_augustus("AT-AC") == "atac"
        assert _splice_type_to_augustus("TT-CA") == "ttca"

    def test_unknown_returns_none(self):
        from eukan.annotation.augustus import _splice_type_to_augustus

        assert _splice_type_to_augustus("unknown") is None

    def test_malformed_returns_none(self):
        from eukan.annotation.augustus import _splice_type_to_augustus

        assert _splice_type_to_augustus("GTT-AG") is None  # 3-char donor
        assert _splice_type_to_augustus("G-A") is None     # 1-char


# ---------------------------------------------------------------------------
# _get_splice_sites_flag
# ---------------------------------------------------------------------------


class TestGetSpliceSitesFlag:
    def _make_config(self, tmp_path, allow_noncanonical=False):
        from eukan.settings import PipelineConfig

        genome = tmp_path / "genome.fa"
        genome.write_text(">chr1\nACGT\n")
        return PipelineConfig(
            genome=genome,
            proteins=[genome],
            work_dir=tmp_path,
            allow_noncanonical_splice=allow_noncanonical,
        )

    def test_no_summary_no_flag(self, tmp_path):
        """Without splice_site_summary.json and no --splice-permissive, no flag."""
        from eukan.annotation.augustus import _get_splice_sites_flag

        config = self._make_config(tmp_path)
        assert _get_splice_sites_flag(config) == []

    def test_no_summary_permissive_blanket(self, tmp_path):
        """With --splice-permissive but no summary, fall back to atac."""
        from eukan.annotation.augustus import _get_splice_sites_flag

        config = self._make_config(tmp_path, allow_noncanonical=True)
        result = _get_splice_sites_flag(config)
        assert len(result) == 1
        assert "atac" in result[0]

    def test_summary_with_sufficient_evidence(self, tmp_path):
        """Types passing both absolute and proportional thresholds should be included."""
        from eukan.annotation.augustus import _get_splice_sites_flag

        # 5000 total junctions; AT-AC at 100 = 2%, TT-CA at 60 = 1.2%
        # Both above _MIN_JUNCTIONS=10, _MIN_UNIQUE_READS=50, _MIN_FRACTION=1%
        summary = {
            "GT-AG": {"count": 4840, "unique_reads": 500000},
            "AT-AC": {"count": 100, "unique_reads": 200},
            "TT-CA": {"count": 60, "unique_reads": 80},
        }
        (tmp_path / "splice_site_summary.json").write_text(json.dumps(summary))
        config = self._make_config(tmp_path)
        result = _get_splice_sites_flag(config)

        assert len(result) == 1
        flag_value = result[0].split("=")[1]
        types = set(flag_value.split(","))
        assert "gtag" not in types  # builtin
        assert "atac" in types
        assert "ttca" in types

    def test_summary_below_absolute_threshold(self, tmp_path):
        """Types below absolute count threshold should be excluded."""
        from eukan.annotation.augustus import _get_splice_sites_flag

        summary = {
            "GT-AG": {"count": 500, "unique_reads": 10000},
            "AT-AC": {"count": 3, "unique_reads": 10},  # below _MIN_JUNCTIONS
        }
        (tmp_path / "splice_site_summary.json").write_text(json.dumps(summary))
        config = self._make_config(tmp_path)
        assert _get_splice_sites_flag(config) == []

    def test_summary_below_proportional_threshold(self, tmp_path):
        """Types passing absolute but failing proportional threshold should be excluded."""
        from eukan.annotation.augustus import _get_splice_sites_flag

        # AT-AC: 15 junctions (passes _MIN_JUNCTIONS=10) but
        # 15/10015 = 0.15% (fails _MIN_FRACTION=1%)
        summary = {
            "GT-AG": {"count": 10000, "unique_reads": 5000000},
            "AT-AC": {"count": 15, "unique_reads": 80},
        }
        (tmp_path / "splice_site_summary.json").write_text(json.dumps(summary))
        config = self._make_config(tmp_path)
        assert _get_splice_sites_flag(config) == []

    def test_summary_permissive_lowers_threshold(self, tmp_path):
        """With --splice-permissive, even 1 junction should pass."""
        from eukan.annotation.augustus import _get_splice_sites_flag

        summary = {
            "GT-AG": {"count": 50000, "unique_reads": 1000000},
            "AT-AC": {"count": 1, "unique_reads": 1},
        }
        (tmp_path / "splice_site_summary.json").write_text(json.dumps(summary))
        config = self._make_config(tmp_path, allow_noncanonical=True)
        result = _get_splice_sites_flag(config)
        assert len(result) == 1
        assert "atac" in result[0]

    def test_builtin_types_excluded(self, tmp_path):
        """gcag is already built-in to AUGUSTUS — should not appear in flag."""
        from eukan.annotation.augustus import _get_splice_sites_flag

        summary = {
            "GT-AG": {"count": 50000, "unique_reads": 1000000},
            "GC-AG": {"count": 100, "unique_reads": 500},
        }
        (tmp_path / "splice_site_summary.json").write_text(json.dumps(summary))
        config = self._make_config(tmp_path)
        # Both GT-AG (gtag) and GC-AG (gcag) are builtin
        assert _get_splice_sites_flag(config) == []


# ---------------------------------------------------------------------------
# _analyze_splice_sites
# ---------------------------------------------------------------------------


class TestAnalyzeSpliceSites:
    def _write_genome(self, path: Path) -> None:
        """Write a small genome FASTA with known splice site dinucleotides.

        chr1 is 100bp.  We'll place junctions so that:
          - junction at 11-50: donor = genome[10:12], acceptor = genome[49:51]
          - junction at 61-90: donor = genome[60:62], acceptor = genome[89:91]
        """
        # Build a 100bp sequence with specific dinucleotides at known positions
        seq = list("A" * 100)
        # Junction 1 (pos 11-50): GT-AG (canonical)
        seq[10] = "G"; seq[11] = "T"  # donor at pos 11 (1-based)
        seq[49] = "A"; seq[50] = "G"  # acceptor ending at pos 50 (1-based) — but wait...
        # Actually: acceptor = last 2 bases of intron = genome[end-2:end] (0-based)
        # For intron ending at 50 (1-based): genome[48:50] = seq[48], seq[49]
        seq[48] = "A"; seq[49] = "G"

        # Junction 2 (pos 61-90): AT-AC (non-canonical, motif=0 in STAR)
        seq[60] = "A"; seq[61] = "T"  # donor at pos 61
        seq[88] = "A"; seq[89] = "C"  # acceptor at pos 90: genome[88:90]

        fasta_str = ">chr1\n" + "".join(seq) + "\n"
        path.write_text(fasta_str)

    def _write_sj_tab(self, path: Path) -> None:
        """Write a minimal SJ.out.tab with known junctions."""
        # Columns: chrom, start, end, strand, motif, annotated, uniq, multi, overhang
        lines = [
            # Canonical GT-AG junction (motif=1)
            "chr1\t11\t50\t1\t1\t0\t100\t5\t30\n",
            # Non-canonical junction (motif=0) — dinucleotides extracted from genome
            "chr1\t61\t90\t1\t0\t0\t50\t2\t25\n",
        ]
        path.write_text("".join(lines))

    def test_basic_analysis(self, tmp_path):
        from eukan.assembly.star import _analyze_splice_sites

        genome = tmp_path / "genome.fa"
        self._write_genome(genome)
        sj_file = tmp_path / "STAR_SJ.out.tab"
        self._write_sj_tab(sj_file)

        _analyze_splice_sites(sj_file, genome, tmp_path)

        summary_path = tmp_path / "splice_site_summary.json"
        assert summary_path.exists()

        with open(summary_path) as f:
            summary = json.load(f)

        # Canonical junction uses STAR's motif classification
        assert "GT-AG" in summary
        assert summary["GT-AG"]["count"] == 1
        assert summary["GT-AG"]["unique_reads"] == 100

        # Non-canonical junction has dinucleotides extracted from genome
        assert "AT-AC" in summary
        assert summary["AT-AC"]["count"] == 1
        assert summary["AT-AC"]["unique_reads"] == 50

    def test_empty_sj_file(self, tmp_path):
        from eukan.assembly.star import _analyze_splice_sites

        genome = tmp_path / "genome.fa"
        genome.write_text(">chr1\nACGTACGT\n")
        sj_file = tmp_path / "STAR_SJ.out.tab"
        sj_file.write_text("")

        _analyze_splice_sites(sj_file, genome, tmp_path)

        with open(tmp_path / "splice_site_summary.json") as f:
            summary = json.load(f)
        assert summary == {}
