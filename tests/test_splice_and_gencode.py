"""Tests for non-canonical splice site handling and genetic code integration."""

from __future__ import annotations

import json
import logging
import re
from pathlib import Path

import pytest

from eukan.infra.genetic_code import GeneticCode

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

    def test_non_acgt_returns_none(self):
        """N (assembly gaps) or other non-ACGT bases can't be passed to AUGUSTUS."""
        from eukan.annotation.augustus import _splice_type_to_augustus

        for name in ("NN-AC", "nn-ac", "CT-NC", "AG-NG", "NN-NN", "RY-AG"):
            assert _splice_type_to_augustus(name) is None, name

    def test_lowercase_input(self):
        from eukan.annotation.augustus import _splice_type_to_augustus

        assert _splice_type_to_augustus("gc-ag") == "gcag"


# ---------------------------------------------------------------------------
# _get_splice_sites_flag
# ---------------------------------------------------------------------------


class TestGetSpliceSitesFlag:
    # The N-containing types that the BFL-lab/eukan#8 --splice-permissive run
    # passed to AUGUSTUS, which then rejected the whole option.
    _REPORTED_N_TYPES = (
        "AG-NG", "CT-NC", "CT-NN", "CT-NT", "GT-NN", "NN-AC",
        "NN-AG", "NN-CG", "NN-GC", "NN-NC", "TT-NN",
    )

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

    def test_summary_permissive_keeps_absolute_floors(self, tmp_path):
        """--splice-permissive still rejects types below the junction/read floors."""
        from eukan.annotation.augustus import _get_splice_sites_flag

        summary = {
            "GT-AG": {"count": 50000, "unique_reads": 1000000},
            "AT-AC": {"count": 1, "unique_reads": 1},    # below _MIN_JUNCTIONS
            "TT-CA": {"count": 20, "unique_reads": 10},  # below _MIN_UNIQUE_READS
        }
        (tmp_path / "splice_site_summary.json").write_text(json.dumps(summary))
        config = self._make_config(tmp_path, allow_noncanonical=True)
        assert _get_splice_sites_flag(config) == []

    def test_summary_permissive_drops_proportional_threshold(self, tmp_path):
        """--splice-permissive admits a type the default rejects only for its share."""
        from eukan.annotation.augustus import _get_splice_sites_flag

        # AT-AC: 15 junctions, 200 reads (passes the floors) but
        # 15/10015 = 0.15% (fails _MIN_FRACTION=1%, which permissive drops)
        summary = {
            "GT-AG": {"count": 10000, "unique_reads": 5000000},
            "AT-AC": {"count": 15, "unique_reads": 200},
        }
        (tmp_path / "splice_site_summary.json").write_text(json.dumps(summary))
        config = self._make_config(tmp_path, allow_noncanonical=True)
        assert _get_splice_sites_flag(config) == ["--allow_hinted_splicesites=atac"]

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

    @pytest.mark.parametrize("permissive", [False, True])
    def test_non_acgt_types_excluded(self, tmp_path, permissive):
        """Types with N never reach AUGUSTUS, in either mode, however well supported."""
        from eukan.annotation.augustus import _get_splice_sites_flag

        summary = {
            "GT-AG": {"count": 20000, "unique_reads": 1000000},
            "AT-AC": {"count": 2000, "unique_reads": 10000},
            "NN-AC": {"count": 2000, "unique_reads": 10000},
            "nn-ac": {"count": 2000, "unique_reads": 10000},
            "CT-NC": {"count": 2000, "unique_reads": 10000},
            "AG-NG": {"count": 2000, "unique_reads": 10000},
        }
        (tmp_path / "splice_site_summary.json").write_text(json.dumps(summary))
        config = self._make_config(tmp_path, allow_noncanonical=permissive)
        assert _get_splice_sites_flag(config) == ["--allow_hinted_splicesites=atac"]

    def test_mixed_summary_keeps_only_acgt_types(self, tmp_path, caplog):
        """Only ACGT types are allowed; the dropped N types are named in one warning."""
        from eukan.annotation.augustus import _get_splice_sites_flag

        # 57700 total junctions; every type passes all default thresholds.
        summary = {
            "GT-AG": {"count": 50000, "unique_reads": 1000000},
            "GC-AG": {"count": 2000, "unique_reads": 9000},  # builtin
            "CT-GC": {"count": 1500, "unique_reads": 8000},
            "AT-AC": {"count": 1200, "unique_reads": 4000},
            "CT-NC": {"count": 900, "unique_reads": 4000},
            "NN-AC": {"count": 800, "unique_reads": 4000},
            "TT-CA": {"count": 700, "unique_reads": 3000},
            "unknown": {"count": 600, "unique_reads": 3000},
        }
        (tmp_path / "splice_site_summary.json").write_text(json.dumps(summary))
        config = self._make_config(tmp_path)
        with caplog.at_level(logging.WARNING, logger="eukan.annotation.augustus"):
            result = _get_splice_sites_flag(config)

        assert result == ["--allow_hinted_splicesites=atac,ctgc,ttca"]
        warnings = [
            r.getMessage() for r in caplog.records
            if r.name == "eukan.annotation.augustus" and r.levelno == logging.WARNING
        ]
        assert len(warnings) == 1
        assert "Ignoring 2 splice-site types with non-ACGT bases" in warnings[0]
        assert "CT-NC (900 junctions)" in warnings[0]
        assert "NN-AC (800 junctions)" in warnings[0]
        assert "unknown" not in warnings[0]

    def test_emitted_tokens_are_valid_augustus_patterns(self, tmp_path):
        """Regression for BFL-lab/eukan#8: every emitted pattern is four of a/c/g/t."""
        from eukan.annotation.augustus import _get_splice_sites_flag

        summary = {"GT-AG": {"count": 50000, "unique_reads": 1000000}}
        for splice_type in (*self._REPORTED_N_TYPES, "AT-AC", "CT-GC", "CT-AC", "GT-AT"):
            summary[splice_type] = {"count": 100, "unique_reads": 1000}
        (tmp_path / "splice_site_summary.json").write_text(json.dumps(summary))
        config = self._make_config(tmp_path, allow_noncanonical=True)
        result = _get_splice_sites_flag(config)

        assert len(result) == 1
        flag, value = result[0].split("=")
        assert flag == "--allow_hinted_splicesites"
        tokens = value.split(",")
        assert all(re.fullmatch(r"[acgt]{4}", t) for t in tokens), tokens
        assert sorted(tokens) == ["atac", "ctac", "ctgc", "gtat"]

    @pytest.mark.parametrize("permissive", [False, True])
    def test_only_non_acgt_and_builtin_types_no_flag(self, tmp_path, permissive):
        """Nothing left to allow → no flag, never an empty --allow_hinted_splicesites=."""
        from eukan.annotation.augustus import _get_splice_sites_flag

        summary = {
            "GT-AG": {"count": 20000, "unique_reads": 1000000},
            "GC-AG": {"count": 2000, "unique_reads": 10000},
            "NN-AC": {"count": 2000, "unique_reads": 10000},
            "CT-NC": {"count": 2000, "unique_reads": 10000},
        }
        (tmp_path / "splice_site_summary.json").write_text(json.dumps(summary))
        config = self._make_config(tmp_path, allow_noncanonical=permissive)
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
        from eukan.assembly.align_hints import analyze_splice_sites

        genome = tmp_path / "genome.fa"
        self._write_genome(genome)
        sj_file = tmp_path / "STAR_SJ.out.tab"
        self._write_sj_tab(sj_file)

        analyze_splice_sites(sj_file, genome, tmp_path)

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
        from eukan.assembly.align_hints import analyze_splice_sites

        genome = tmp_path / "genome.fa"
        genome.write_text(">chr1\nACGTACGT\n")
        sj_file = tmp_path / "STAR_SJ.out.tab"
        sj_file.write_text("")

        analyze_splice_sites(sj_file, genome, tmp_path)

        with open(tmp_path / "splice_site_summary.json") as f:
            summary = json.load(f)
        assert summary == {}


# ---------------------------------------------------------------------------
# sanitize_genome_fasta
# ---------------------------------------------------------------------------


class TestSanitizeGenomeFasta:
    def test_genomes_in_one_work_dir_do_not_share_the_copy(self, tmp_path):
        """A second genome with descriptions (e.g. the masked one after the raw one)
        gets its own sanitized copy, not the one built from the first genome."""
        from eukan.validation import sanitize_genome_fasta

        raw = tmp_path / "raw.fasta"
        raw.write_text(">chr1 assembly v1\nACGTACGT\n")
        masked = tmp_path / "raw.masked.fasta"
        masked.write_text(">chr1 assembly v1\nacgtACGT\n")
        work_dir = tmp_path / "annotate"
        work_dir.mkdir()

        first = sanitize_genome_fasta(raw, work_dir)
        second = sanitize_genome_fasta(masked, work_dir)

        assert first != second
        assert first.read_text() == ">chr1\nACGTACGT\n"
        assert second.read_text() == ">chr1\nacgtACGT\n"

    def test_same_name_genome_with_older_mtime_is_resanitized(self, tmp_path):
        """A genome replaced by a same-named file that kept an old mtime (wget, tar,
        cp -p) is sanitized afresh, not served from the earlier copy."""
        import os

        from eukan.validation import sanitize_genome_fasta

        genome = tmp_path / "genome.fasta"
        genome.write_text(">chr1 v1\nAAAA\n")
        sanitize_genome_fasta(genome, tmp_path)
        genome.write_text(">chr1 v2\nCCCC\n")
        os.utime(genome, (0, 0))

        assert sanitize_genome_fasta(genome, tmp_path).read_text() == ">chr1\nCCCC\n"
