"""CLI-boundary error formatting tests for the EukanGroup wrapper."""

from __future__ import annotations

import errno

import click
import pytest
from click.testing import CliRunner
from pydantic import BaseModel, ValidationError, model_validator

from eukan.cli._framework import CONTEXT_SETTINGS, EukanGroup, _format_disk_full
from eukan.exceptions import ExternalToolError

# What a mapper prints before its real error: well over the 500 stderr chars kept.
_PROGRESS = "[M::worker_pipeline::12.345*3.21] mapped 500000 sequences\n" * 40


@pytest.fixture
def app_factory():
    """Build a tiny Click app whose subcommand raises a configured exception."""

    def _build(exc: BaseException):
        @click.group(cls=EukanGroup, context_settings=CONTEXT_SETTINGS)
        def app() -> None:
            pass

        @app.command()
        def boom() -> None:
            raise exc

        return app

    return _build


class TestFormatDiskFullUnit:
    def test_oserror_enospc_matches(self):
        exc = OSError(errno.ENOSPC, "No space left on device", "/work/out.json")
        result = _format_disk_full(exc)
        assert result is not None
        title, details = result
        assert title == "Error: no space left on device"
        assert "Path: /work/out.json" in details

    def test_oserror_edquot_matches(self):
        exc = OSError(errno.EDQUOT, "Disk quota exceeded", "/work/cache.bin")
        result = _format_disk_full(exc)
        assert result is not None

    def test_oserror_other_errno_does_not_match(self):
        exc = OSError(errno.EACCES, "Permission denied", "/etc/shadow")
        assert _format_disk_full(exc) is None


class TestExternalToolErrorDiskFull:
    """format_for_cli reports a tool that ran out of disk space as such."""

    def test_enospc_stderr_gets_the_disk_full_title(self):
        exc = ExternalToolError(
            "STAR failed",
            tool="STAR",
            returncode=137,
            cmd=["STAR", "--runMode", "alignReads"],
            stderr_snippet="EXITING because of FATAL ERROR: No space left on device\n",
            step="annotation/star",
        )
        title, details = exc.format_for_cli()
        assert title == "Error: no space left on device"
        assert "Tool: STAR (exit 137)" in details
        assert "Step: annotation/star" in details
        assert details[-1] == (
            "Free up disk space on the device hosting the work directory and re-run."
        )

    def test_quota_stderr_gets_the_disk_full_title(self):
        exc = ExternalToolError(
            "rnaspades failed",
            tool="rnaSPAdes",
            returncode=1,
            stderr_snippet="write error: Disk quota exceeded",
        )
        assert exc.format_for_cli()[0] == "Error: no space left on device"

    def test_enospc_after_progress_output_gets_the_disk_full_title(self):
        """The error line comes last, after the progress output; the snippet must
        keep it for the detector to see it."""
        exc = ExternalToolError(
            "samtools failed (exit 1)", tool="samtools", returncode=1,
            stderr_snippet=_PROGRESS + "samtools sort: No space left on device\n",
        )
        assert exc.format_for_cli()[0] == "Error: no space left on device"

    def test_samtools_illegal_seek_gets_the_disk_full_title(self):
        """htslib reports a write to a full disk as 'Illegal seek'; the CLI used to
        miss it, although assembly/bam_utils.py already treated it as disk-full."""
        exc = ExternalToolError(
            "samtools failed (exit 1)", tool="samtools", returncode=1,
            stderr_snippet='samtools sort: failed writing to "o.bam": Illegal seek\n',
        )
        assert exc.format_for_cli()[0] == "Error: no space left on device"

    def test_unrelated_stderr_keeps_the_generic_title(self):
        exc = ExternalToolError(
            "STAR failed",
            tool="STAR",
            returncode=1,
            stderr_snippet="segfault at 0xdeadbeef",
        )
        assert exc.format_for_cli()[0] == "Error: STAR failed (exit 1)"


class TestExternalToolErrorFormat:
    def test_shows_message_and_stderr_tail(self):
        """A wrapper's explanation (bam_utils' sort re-wrap) and the error at the end
        of a chatty stderr both reach the user; a generic message is not repeated."""
        exc = ExternalToolError(
            "samtools sort failed writing the coordinate-sorted BAM",
            tool="samtools", returncode=1,
            stderr_snippet=_PROGRESS + 'samtools sort: failed writing to "o.bam": Illegal seek\n',
        )
        details = "\n".join(exc.format_for_cli()[1])
        assert "samtools sort failed writing the coordinate-sorted BAM" in details
        assert "Illegal seek" in details

        generic = ExternalToolError("samtools failed (exit 1)", tool="samtools", returncode=1)
        assert "samtools failed (exit 1)" not in generic.format_for_cli()[1]


class TestCliBoundary:
    def test_enospc_oserror_prints_clean_message_no_traceback(self, app_factory):
        exc = OSError(errno.ENOSPC, "No space left on device", "/work/out.json")
        runner = CliRunner()
        result = runner.invoke(app_factory(exc), ["boom"])
        assert result.exit_code == 1
        assert "Error: no space left on device" in result.output
        assert "/work/out.json" in result.output
        # The clean handler must not surface a Python traceback.
        assert "Traceback" not in result.output

    def test_enospc_external_tool_prints_clean_message(self, app_factory):
        exc = ExternalToolError(
            "STAR failed",
            tool="STAR",
            returncode=137,
            cmd=["STAR"],
            stderr_snippet="No space left on device",
            step="annotation/star",
        )
        runner = CliRunner()
        result = runner.invoke(app_factory(exc), ["boom"])
        assert result.exit_code == 1
        assert "Error: no space left on device" in result.output
        assert "STAR" in result.output
        # A disk-full report points at freeing space, not at -v.
        assert "Run with -v" not in result.output

    @pytest.mark.parametrize("stderr", [
        "samtools sort: No space left on device",
        "[E::bgzf_flush] File write failed (wrong size)",
    ])
    def test_disk_full_tool_error_keeps_message_and_hint(self, app_factory, stderr):
        """The disk-full report used to drop the raiser's explanation and hint
        (bam_utils' samtools-sort re-wrap carries both), and missed htslib's
        write-failure wording altogether."""
        exc = ExternalToolError(
            "samtools sort failed writing the coordinate-sorted BAM",
            tool="samtools", returncode=1, cmd=["samtools", "sort"],
            stderr_snippet=f"{stderr}\n", step="assembly/minimap2",
            hint="Free space on the run-dir filesystem, then re-run.",
        )
        result = CliRunner().invoke(app_factory(exc), ["boom"])
        assert result.exit_code == 1
        assert result.output == (
            "Error: no space left on device\n"
            "  Tool: samtools (exit 1)\n"
            "  samtools sort failed writing the coordinate-sorted BAM\n"
            "  Step: assembly/minimap2\n"
            f"  stderr: {stderr}\n\n"
            "  Hint: Free space on the run-dir filesystem, then re-run.\n"
        )

    def test_unrelated_oserror_still_propagates(self, app_factory):
        exc = OSError(errno.EACCES, "Permission denied", "/etc/shadow")
        runner = CliRunner()
        result = runner.invoke(app_factory(exc), ["boom"])
        # EACCES is not one of our domain errors, so it propagates as-is.
        assert result.exit_code != 0
        assert "no space left on device" not in result.output.lower()


class _Reads(BaseModel):
    left: int = 0

    @model_validator(mode="after")
    def _paired(self) -> _Reads:
        if self.left < 0:
            raise ValueError("Paired-end mode requires both --left and --right.")
        return self


def _validation_error(**kwargs) -> ValidationError:
    with pytest.raises(ValidationError) as info:
        _Reads(**kwargs)
    return info.value


class TestConfigErrors:
    def test_model_level_error_prints_without_an_empty_location(self, app_factory):
        result = CliRunner().invoke(app_factory(_validation_error(left=-1)), ["boom"])
        assert result.exit_code == 1
        assert "Error: invalid configuration" in result.output
        assert "  Value error, Paired-end mode requires both" in result.output
        assert "  : " not in result.output

    def test_field_error_keeps_its_location(self, app_factory):
        result = CliRunner().invoke(app_factory(_validation_error(left="x")), ["boom"])
        assert result.exit_code == 1
        assert "  left: Input should be a valid integer" in result.output


class TestEukanErrors:
    def test_title_details_and_hint(self, app_factory):
        exc = ExternalToolError(
            "spaln crashed", tool="spaln", returncode=2, cmd=["spaln"],
            stderr_snippet="bad\n", step="annotation/prot_align",
        )
        exc.hint = "Check the protein FASTA."
        result = CliRunner().invoke(app_factory(exc), ["boom"])
        assert result.exit_code == 1
        assert result.output == (
            "Error: spaln failed (exit 2)\n"
            "  spaln crashed\n"
            "  Step: annotation/prot_align\n"
            "  stderr: bad\n\n"
            "  Run with -v for full command and stderr output.\n"
            "  Hint: Check the protein FASTA.\n"
        )

    @pytest.mark.parametrize("exc", [SystemExit(3), click.exceptions.Exit(3)])
    def test_exits_pass_through(self, app_factory, exc):
        assert CliRunner().invoke(app_factory(exc), ["boom"]).exit_code == 3
