"""CLI-boundary error formatting tests for the EukanGroup wrapper."""

from __future__ import annotations

import errno

import click
import pytest
from click.testing import CliRunner

from eukan.cli._framework import CONTEXT_SETTINGS, EukanGroup, _format_disk_full
from eukan.exceptions import ExternalToolError


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

    def test_external_tool_with_enospc_stderr_matches(self):
        exc = ExternalToolError(
            "STAR failed",
            tool="STAR",
            returncode=137,
            cmd=["STAR", "--runMode", "alignReads"],
            stderr_snippet="EXITING because of FATAL ERROR: No space left on device\n",
            step="annotation/star",
        )
        result = _format_disk_full(exc)
        assert result is not None
        title, details = result
        assert title == "Error: no space left on device"
        assert any("STAR" in d for d in details)
        assert any("annotation/star" in d for d in details)

    def test_external_tool_with_quota_stderr_matches(self):
        exc = ExternalToolError(
            "rnaspades failed",
            tool="rnaSPAdes",
            returncode=1,
            stderr_snippet="write error: Disk quota exceeded",
        )
        assert _format_disk_full(exc) is not None

    def test_external_tool_with_unrelated_stderr_does_not_match(self):
        exc = ExternalToolError(
            "STAR failed",
            tool="STAR",
            returncode=1,
            stderr_snippet="segfault at 0xdeadbeef",
        )
        assert _format_disk_full(exc) is None


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
        # The standard ExternalToolError details should NOT be shown — the
        # disk-full handler short-circuits first.
        assert "Run with -v" not in result.output

    def test_unrelated_oserror_still_propagates(self, app_factory):
        exc = OSError(errno.EACCES, "Permission denied", "/etc/shadow")
        runner = CliRunner()
        result = runner.invoke(app_factory(exc), ["boom"])
        # EACCES is not one of our domain errors, so it propagates as-is.
        assert result.exit_code != 0
        assert "no space left on device" not in result.output.lower()
