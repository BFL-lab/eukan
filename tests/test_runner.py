"""Tests for eukan.infra.runner — subprocess execution + Ctrl-C handling."""

from __future__ import annotations

import os
import subprocess

import pytest

from eukan.exceptions import ExternalToolError, MissingToolError
from eukan.infra import runner
from eukan.infra.runner import (
    _RUNNING,
    _track,
    _untrack,
    run_cmd,
    run_piped,
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

    def test_missing_binary_raises_missing_tool_error(self, tmp_path):
        """Tools not on PATH should raise MissingToolError, not FileNotFoundError."""
        with pytest.raises(MissingToolError) as exc_info:
            run_cmd(["definitely-not-a-real-binary-xyz"], cwd=tmp_path)
        assert exc_info.value.tool == "definitely-not-a-real-binary-xyz"
        assert "not found on PATH" in str(exc_info.value)
        assert exc_info.value.hint and "eukan check" in exc_info.value.hint

    def test_missing_binary_with_path_uses_basename(self, tmp_path):
        """Path-style commands surface the binary basename, not the full path."""
        with pytest.raises(MissingToolError) as exc_info:
            run_cmd(["/opt/nope/missing-tool"], cwd=tmp_path)
        assert exc_info.value.tool == "missing-tool"

    def test_missing_cwd_does_not_become_missing_tool_error(self, tmp_path):
        """A bogus cwd should still surface as FileNotFoundError, not be misclassified."""
        bogus_cwd = tmp_path / "does_not_exist"
        with pytest.raises(FileNotFoundError):
            run_cmd(["true"], cwd=bogus_cwd)

    def test_non_utf8_stderr_does_not_crash_communicate(self, tmp_path):
        """External tools emitting raw bytes on stderr must not crash run_cmd."""
        # 0xf9 is invalid as the start of a UTF-8 sequence (the actual byte
        # AUGUSTUS surfaced in the wild).
        run_cmd(
            ["python", "-c", "import sys; sys.stderr.buffer.write(b'pre \\xf9 post\\n')"],
            cwd=tmp_path,
        )

    def test_non_utf8_stderr_on_failure_surfaces_replaced_text(self, tmp_path):
        """A failing tool with non-UTF-8 stderr should still produce a usable snippet."""
        with pytest.raises(ExternalToolError) as exc_info:
            run_cmd(
                ["python", "-c",
                 "import sys; sys.stderr.buffer.write(b'BOOM \\xf9 BOOM\\n'); sys.exit(2)"],
                cwd=tmp_path,
            )
        # The undecodable byte is replaced rather than crashing the runner.
        assert exc_info.value.returncode == 2
        assert "BOOM" in exc_info.value.stderr_snippet


class TestRunCmdResourceCleanup:
    """Redirect file handles must not leak if setup raises before Popen."""

    @staticmethod
    def _fd_count() -> int:
        return len(os.listdir(f"/proc/{os.getpid()}/fd"))

    @pytest.mark.skipif(
        not os.path.isdir(f"/proc/{os.getpid()}/fd"),
        reason="needs /proc fd listing",
    )
    def test_setup_failure_does_not_leak_redirect_handles(self, tmp_path, monkeypatch):
        """A raise between opening the out/err files and starting the child
        must still close them. Pre-ExitStack this leaked 2 fds per call."""
        def boom(_cmd):
            raise RuntimeError("tool-name boom")

        monkeypatch.setattr(runner, "_tool_name", boom)

        baseline = self._fd_count()
        for _ in range(40):
            with pytest.raises(RuntimeError, match="tool-name boom"):
                run_cmd(["true"], cwd=tmp_path, out_file="o.txt", err_file="e.txt")
        # No monotonic growth: handles were closed each time.
        assert self._fd_count() <= baseline + 1


class TestRunPiped:
    def test_missing_first_command_raises_missing_tool_error(self, tmp_path):
        with pytest.raises(MissingToolError) as exc_info:
            run_piped(["definitely-not-a-real-binary-xyz"], ["cat"], cwd=tmp_path)
        assert exc_info.value.tool == "definitely-not-a-real-binary-xyz"

    def test_missing_second_command_raises_missing_tool_error(self, tmp_path):
        with pytest.raises(MissingToolError) as exc_info:
            run_piped(["true"], ["definitely-not-a-real-binary-xyz"], cwd=tmp_path)
        assert exc_info.value.tool == "definitely-not-a-real-binary-xyz"

    def test_producer_failure_surfaces_even_when_consumer_succeeds(self, tmp_path):
        # cat drains the input and exits 0, but the producer exits non-zero —
        # this must NOT be silently swallowed (the segemehl-OOM-masking bug).
        with pytest.raises(ExternalToolError) as exc_info:
            run_piped(["bash", "-c", "echo hi; exit 7"], ["cat"], cwd=tmp_path)
        assert exc_info.value.returncode == 7
        assert exc_info.value.tool == "bash"

    def test_consumer_failure_reports_producer_exit(self, tmp_path):
        # Consumer fails on an empty stream; the producer's non-zero exit is
        # folded into the error message so the real cause is visible.
        with pytest.raises(ExternalToolError) as exc_info:
            run_piped(
                ["bash", "-c", "exit 9"],
                ["bash", "-c", "head -c1 >/dev/null; exit 1"],
                cwd=tmp_path,
            )
        assert "exited 9" in exc_info.value.stderr_snippet

    def test_chatty_producer_stderr_does_not_deadlock(self, tmp_path):
        # >64 KB on the producer's stderr would deadlock an unread PIPE; the
        # temp-file drain must let the pipeline complete normally.
        out = run_piped(
            ["python3", "-c", "import sys; sys.stderr.write('X' * 200000); print('data')"],
            ["cat"],
            cwd=tmp_path,
        )
        assert out.strip() == "data"

    def test_producer_sigpipe_is_not_an_error(self, tmp_path):
        # `yes` produces forever; the consumer reads one byte and exits 0, so
        # the producer dies of SIGPIPE. That early-close is legitimate, not a
        # failure, and must not raise.
        out = run_piped(["yes"], ["head", "-c", "1"], cwd=tmp_path)
        assert out == "y"
