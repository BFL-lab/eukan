"""Tests for the prep-submission table2asn wrapper."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from pydantic import ValidationError

from eukan.exceptions import ConfigurationError
from eukan.settings import SubmissionConfig
from eukan.submission import build_command, shell_repr
from eukan.submission.pipeline import _parse_val_report

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def submission_inputs(tmp_path: Path) -> dict[str, Path]:
    """Create minimal valid input files (genome, gff3, template) under tmp_path."""
    genome = tmp_path / "genome.fasta"
    genome.write_text(">chr1\nACGT\n")
    gff3 = tmp_path / "final.mod.gff3"
    gff3.write_text("##gff-version 3\nchr1\teukan\tgene\t1\t4\t.\t+\t.\tID=g1\n")
    template = tmp_path / "submission-template.sbt"
    template.write_text("Submit-block ::= { }\n")
    return {"genome": genome, "gff3": gff3, "template": template, "work_dir": tmp_path}


def _base_config(inputs: dict[str, Path], **overrides) -> SubmissionConfig:
    """Build a SubmissionConfig with sensible defaults for tests."""
    kwargs = dict(
        work_dir=inputs["work_dir"],
        genome=inputs["genome"],
        gff3=inputs["gff3"],
        template=inputs["template"],
        organism="Foo bar",
        outdir=inputs["work_dir"] / "submission",
    )
    kwargs.update(overrides)
    return SubmissionConfig(**kwargs)


# ---------------------------------------------------------------------------
# build_command
# ---------------------------------------------------------------------------


class TestBuildCommand:
    def test_includes_fixed_flags(self, submission_inputs):
        cmd = build_command(_base_config(submission_inputs))
        for flag in ("-split-logs", "-W", "-J", "-Z", "-euk", "-T"):
            assert flag in cmd
        assert cmd[:1] == ["table2asn"]
        # -V b is a flag/value pair
        assert "-V" in cmd
        assert cmd[cmd.index("-V") + 1] == "b"

    def test_organism_only_composes_source_info(self, submission_inputs):
        cmd = build_command(_base_config(submission_inputs, organism="Homo sapiens"))
        j_value = cmd[cmd.index("-j") + 1]
        assert j_value == "[organism=Homo sapiens]"

    def test_organism_and_isolate_compose_source_info(self, submission_inputs):
        cmd = build_command(_base_config(
            submission_inputs, organism="Homo sapiens", isolate="ABC123",
        ))
        assert cmd[cmd.index("-j") + 1] == "[organism=Homo sapiens] [isolate=ABC123]"

    def test_source_info_overrides_structured_qualifiers(self, submission_inputs):
        cmd = build_command(_base_config(
            submission_inputs,
            organism="ignored",
            isolate="ignored",
            source_info="[organism=Foo] [country=Canada]",
        ))
        assert cmd[cmd.index("-j") + 1] == "[organism=Foo] [country=Canada]"

    def test_missing_source_qualifier_raises(self, submission_inputs):
        config = _base_config(submission_inputs, organism=None, isolate=None)
        # source_info is also None; build_command should raise
        with pytest.raises(ConfigurationError, match="source qualifier"):
            build_command(config)

    def test_locus_tag_prefix_appended(self, submission_inputs):
        cmd = build_command(_base_config(submission_inputs, locus_tag_prefix="ABC"))
        assert "-locus-tag-prefix" in cmd
        assert cmd[cmd.index("-locus-tag-prefix") + 1] == "ABC"

    def test_extra_args_round_trip(self, submission_inputs):
        cmd = build_command(_base_config(
            submission_inputs, extra_args=["-split-dr", "-huge"],
        ))
        # extra args should appear at the end, in order
        assert cmd[-2:] == ["-split-dr", "-huge"]

    def test_paths_serialized_as_strings(self, submission_inputs):
        cmd = build_command(_base_config(submission_inputs))
        # -i, -f, -t, -o, -outdir all take string paths
        for flag in ("-i", "-f", "-t", "-o", "-outdir"):
            value = cmd[cmd.index(flag) + 1]
            assert isinstance(value, str)

    def test_default_cleanup_mode_assembly(self, submission_inputs):
        cmd = build_command(_base_config(submission_inputs))
        assert cmd[cmd.index("-c") + 1] == "befw"
        assert cmd[cmd.index("-M") + 1] == "n"
        assert cmd[cmd.index("-a") + 1] == "r10k"

    def test_overridden_cleanup_mode_assembly(self, submission_inputs):
        cmd = build_command(_base_config(
            submission_inputs, cleanup="bef", mode="b", assembly_type="r5k",
        ))
        assert cmd[cmd.index("-c") + 1] == "bef"
        assert cmd[cmd.index("-M") + 1] == "b"
        assert cmd[cmd.index("-a") + 1] == "r5k"


# ---------------------------------------------------------------------------
# Auto-discovery
# ---------------------------------------------------------------------------


class TestAutoDiscovery:
    def test_default_output_under_outdir(self, submission_inputs):
        config = _base_config(submission_inputs)
        assert config.output == config.outdir / f"{config.genome.stem}.sqn"

    def test_explicit_output_preserved(self, submission_inputs):
        explicit = submission_inputs["work_dir"] / "custom.sqn"
        config = _base_config(submission_inputs, output=explicit)
        assert config.output == explicit

    def test_discovers_genome_from_manifest(self, submission_inputs):
        # Construct a minimal eukan-run.json with a genome path, omit explicit
        # genome arg, and check the validator picks it up.
        work_dir = submission_inputs["work_dir"]
        manifest_path = work_dir / "eukan-run.json"
        manifest_path.write_text(json.dumps({
            "version": "1",
            "status": "completed",
            "started_at": "2026-05-06T00:00:00",
            "genome": str(submission_inputs["genome"]),
            "proteins": [],
            "kingdom": None,
            "genetic_code": "1",
            "num_cpu": 1,
            "has_transcripts": False,
            "tool_versions": {},
            "steps": {},
        }))
        config = SubmissionConfig(
            work_dir=work_dir,
            template=submission_inputs["template"],
            organism="Foo bar",
            gff3=submission_inputs["gff3"],
        )
        assert config.genome == submission_inputs["genome"]

    def test_discovers_gff3_prefers_mod(self, submission_inputs, tmp_path: Path):
        work_dir = tmp_path / "wd"
        work_dir.mkdir()
        # Both files present; .mod.gff3 should win.
        plain = work_dir / "final.gff3"
        plain.write_text("##gff-version 3\n")
        mod = work_dir / "final.mod.gff3"
        mod.write_text("##gff-version 3\n")
        config = SubmissionConfig(
            work_dir=work_dir,
            genome=submission_inputs["genome"],
            template=submission_inputs["template"],
            organism="Foo bar",
        )
        assert config.gff3 == mod

    def test_discovers_gff3_falls_back_to_plain(self, submission_inputs, tmp_path: Path):
        work_dir = tmp_path / "wd"
        work_dir.mkdir()
        plain = work_dir / "final.gff3"
        plain.write_text("##gff-version 3\n")
        config = SubmissionConfig(
            work_dir=work_dir,
            genome=submission_inputs["genome"],
            template=submission_inputs["template"],
            organism="Foo bar",
        )
        assert config.gff3 == plain

    def test_missing_genome_raises(self, submission_inputs, tmp_path: Path):
        # work_dir has no eukan-run.json; no --genome → ValidationError.
        empty = tmp_path / "empty"
        empty.mkdir()
        with pytest.raises(ValidationError, match="genome"):
            SubmissionConfig(
                work_dir=empty,
                template=submission_inputs["template"],
                gff3=submission_inputs["gff3"],
                organism="Foo bar",
            )

    def test_missing_gff3_raises(self, submission_inputs, tmp_path: Path):
        empty = tmp_path / "empty"
        empty.mkdir()
        with pytest.raises(ValidationError, match="gff3"):
            SubmissionConfig(
                work_dir=empty,
                template=submission_inputs["template"],
                genome=submission_inputs["genome"],
                organism="Foo bar",
            )


# ---------------------------------------------------------------------------
# Validation report parser
# ---------------------------------------------------------------------------


class TestParseValReport:
    def test_counts_severities(self, tmp_path: Path):
        val = tmp_path / "out.val"
        val.write_text(
            "ERROR: valid [SEQ_FEAT.MissingCDSproduct] No CDS product\n"
            "WARNING: valid [SEQ_INST.LowQualityScore] Low quality\n"
            "ERROR: valid [SEQ_FEAT.PartialProblem] Partial issue\n"
            "FATAL: valid [SEQ_INST.SeqDataLenWrong] Length mismatch\n"
            "INFO: valid [SEQ_INST.Note] Just FYI\n",
        )
        counts = _parse_val_report(val)
        assert counts == {"FATAL": 1, "ERROR": 2, "WARNING": 1, "INFO": 1}

    def test_missing_file_returns_zero_counts(self, tmp_path: Path):
        counts = _parse_val_report(tmp_path / "does-not-exist.val")
        assert counts == {"FATAL": 0, "ERROR": 0, "WARNING": 0, "INFO": 0}

    def test_empty_file_returns_zero_counts(self, tmp_path: Path):
        empty = tmp_path / "empty.val"
        empty.touch()
        assert _parse_val_report(empty) == {
            "FATAL": 0, "ERROR": 0, "WARNING": 0, "INFO": 0,
        }


# ---------------------------------------------------------------------------
# shell_repr
# ---------------------------------------------------------------------------


def test_shell_repr_quotes_arguments_with_spaces(submission_inputs):
    cmd = build_command(_base_config(submission_inputs, organism="Homo sapiens"))
    rendered = shell_repr(cmd)
    # The -j value contains a space and brackets; both must be quoted.
    assert "'[organism=Homo sapiens]'" in rendered
