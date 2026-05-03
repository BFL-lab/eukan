"""Tests for eukan.infra.runner — subprocess execution + Ctrl-C handling."""

from __future__ import annotations

import subprocess

import pytest

from eukan.exceptions import ExternalToolError
from eukan.infra.runner import (
    _RUNNING,
    _track,
    _untrack,
    run_cmd,
    terminate_all_children,
)


class TestTerminateAllChildren:
    def test_no_running_processes_is_noop(self):
        assert set() == _RUNNING
        terminate_all_children()  # must not raise

    def test_terminates_tracked_process(self):
        """A long-running tracked child should be terminated within the grace period."""
        proc = subprocess.Popen(
            ["sleep", "60"],
            start_new_session=True,
        )
        _track(proc)
        try:
            assert proc.returncode is None
            terminate_all_children(grace_period=2.0)
            assert proc.returncode is not None
            assert proc.poll() is not None
        finally:
            _untrack(proc)
            if proc.returncode is None:
                proc.kill()
                proc.wait()

    def test_kills_process_that_ignores_sigterm(self):
        """A child that traps SIGTERM is SIGKILLed once the grace period elapses."""
        proc = subprocess.Popen(
            ["python", "-c", "import signal, time; signal.signal(signal.SIGTERM, signal.SIG_IGN); time.sleep(60)"],
            start_new_session=True,
        )
        _track(proc)
        try:
            terminate_all_children(grace_period=1.0)
            assert proc.returncode is not None
        finally:
            _untrack(proc)
            if proc.returncode is None:
                proc.kill()
                proc.wait()


class TestRunCmd:
    def test_in_file_routes_stdin(self, tmp_path):
        """run_cmd should pipe a file into the child's stdin via in_file=."""
        (tmp_path / "input.txt").write_text("line1\nline2\nline3\n")
        run_cmd(
            ["wc", "-l"],
            cwd=tmp_path,
            in_file="input.txt",
            out_file="count.txt",
        )
        assert (tmp_path / "count.txt").read_text().strip().split()[0] == "3"

    def test_out_file_streams_directly(self, tmp_path):
        """out_file should land bytes from child stdout in the named file."""
        run_cmd(
            ["printf", "hello"],
            cwd=tmp_path,
            out_file="out.txt",
        )
        assert (tmp_path / "out.txt").read_bytes() == b"hello"

    def test_err_file_streams_directly(self, tmp_path):
        """err_file should capture stderr without it appearing in the exception."""
        run_cmd(
            ["sh", "-c", "echo oops 1>&2"],
            cwd=tmp_path,
            err_file="err.txt",
        )
        assert (tmp_path / "err.txt").read_bytes() == b"oops\n"

    def test_nonzero_exit_raises(self, tmp_path):
        with pytest.raises(ExternalToolError) as exc_info:
            run_cmd(["false"], cwd=tmp_path)
        assert exc_info.value.returncode != 0

    def test_timeout_raises(self, tmp_path):
        with pytest.raises(ExternalToolError):
            run_cmd(["sleep", "10"], cwd=tmp_path, timeout=1)

    def test_process_is_untracked_after_completion(self, tmp_path):
        before = len(_RUNNING)
        run_cmd(["true"], cwd=tmp_path)
        assert len(_RUNNING) == before
