"""Safe subprocess execution for external bioinformatics tools."""

from __future__ import annotations

import contextlib
import os
import signal
import subprocess
import tempfile
import threading
import time
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

from eukan.exceptions import ExternalToolError, MissingToolError
from eukan.infra.environ import subprocess_env as _subprocess_env
from eukan.infra.logging import get_logger

log = get_logger(__name__)


# Registry of currently-running child processes so a top-level SIGINT
# handler (installed by the CLI) can terminate them on shutdown.
_REGISTRY_LOCK = threading.Lock()
_RUNNING: set[subprocess.Popen] = set()


def _track(proc: subprocess.Popen) -> None:
    with _REGISTRY_LOCK:
        _RUNNING.add(proc)


def _untrack(proc: subprocess.Popen) -> None:
    with _REGISTRY_LOCK:
        _RUNNING.discard(proc)


def _signal_group(proc: subprocess.Popen, sig: int) -> None:
    """Send *sig* to *proc*'s whole process group, not only to *proc*.

    Every child starts in its own session (see :func:`_tracked_popen`), so its
    process-group id is its PID and the group also holds the worker processes
    the tool forks, which signalling the PID alone would leave running. A
    group that is already gone, or that may not be signalled, is ignored.
    """
    with contextlib.suppress(ProcessLookupError, PermissionError):
        os.killpg(proc.pid, sig)


@contextmanager
def _tracked_popen(
    cmd: list[str],
    *,
    cwd: Path,
    env: dict[str, str] | None,
    missing_tool: str | None = None,
    **popen_kwargs,
) -> Iterator[subprocess.Popen]:
    """Spawn a Popen registered for SIGINT cleanup, with isolated session.

    On ``KeyboardInterrupt`` the child's process group is terminated (5s
    grace) then SIGKILLed; in either case the process is unregistered before the
    context exits. ``FileNotFoundError`` is translated to
    :class:`MissingToolError` when *missing_tool* is provided and the
    cwd actually exists (so a missing cwd surfaces as the original
    FileNotFoundError, not a misleading "tool not found").
    """
    try:
        proc = subprocess.Popen(
            cmd, cwd=cwd, env=env, start_new_session=True, **popen_kwargs,
        )
    except FileNotFoundError as exc:
        if missing_tool is not None and cwd.is_dir():
            raise MissingToolError(missing_tool) from exc
        raise
    _track(proc)
    try:
        yield proc
    except KeyboardInterrupt:
        _signal_group(proc, signal.SIGTERM)
        with contextlib.suppress(subprocess.TimeoutExpired):
            proc.communicate(timeout=5)
        if proc.returncode is None:
            _signal_group(proc, signal.SIGKILL)
        raise
    finally:
        _untrack(proc)


def terminate_all_children(grace_period: float = 5.0) -> None:
    """Terminate every tracked child's process group, then SIGKILL stragglers.

    Called from the CLI's SIGINT handler so Ctrl-C doesn't orphan running
    bioinformatics tools or the worker processes they fork.
    """
    with _REGISTRY_LOCK:
        procs = list(_RUNNING)
    if not procs:
        return
    log.warning("Interrupted; terminating %d running child process(es)...", len(procs))
    for p in procs:
        _signal_group(p, signal.SIGTERM)
    deadline = time.monotonic() + grace_period
    for p in procs:
        remaining = max(0.0, deadline - time.monotonic())
        try:
            p.wait(timeout=remaining)
        except subprocess.TimeoutExpired:
            _signal_group(p, signal.SIGKILL)


def _tool_name(cmd: list[str]) -> str:
    """Extract the tool binary name from a command list.

    ``["STAR", ...]`` returns ``"STAR"``.
    """
    for token in cmd:
        if token.startswith("-"):
            continue
        return Path(token).name
    return cmd[0] if cmd else "unknown"


def run_cmd(
    cmd: list[str],
    *,
    cwd: Path,
    in_file: str | None = None,
    out_file: str | None = None,
    binary: bool = False,
    extra_env: dict[str, str] | None = None,
) -> None:
    """Run an external command safely.

    Args:
        cmd: Command as a list of strings (never a shell string).
        cwd: Working directory for the subprocess.
        in_file: If set, the file at this path within *cwd* is opened
            for reading and connected to the child's stdin (no Python
            buffering — the child reads from the fd directly).
        out_file: If set, stdout is streamed straight to this filename
            within *cwd* via the child's fd (no in-process capture);
            otherwise stdout goes to ``/dev/null``. Stdout is never
            captured, so nothing is returned.
        binary: Read the child's stderr as bytes and decode it here,
            instead of opening the pipe in text mode.
        extra_env: Additional environment variables merged on top of
            the current environment for this subprocess only.

    Raises:
        ExternalToolError: If the command exits non-zero.
    """
    cmd_str = " ".join(cmd)
    log.debug("Running: %s (cwd=%s)", cmd_str, cwd)

    # An ExitStack owns the redirect file handles so they're closed
    # unconditionally — even if setup between opening them and starting the
    # subprocess (fileno(), env build, Popen) raises before _tracked_popen's
    # own cleanup could run.
    with contextlib.ExitStack() as stack:
        in_handle = stack.enter_context(open(cwd / in_file, "rb")) if in_file else None
        out_handle = stack.enter_context(open(cwd / out_file, "wb")) if out_file else None

        stdin_target: int | None = in_handle.fileno() if in_handle else None
        stdout_target: int = out_handle.fileno() if out_handle else subprocess.DEVNULL

        tool = _tool_name(cmd)
        # External tools occasionally emit non-UTF-8 bytes on stderr (e.g.
        # AUGUSTUS surfaces raw locale-encoded bytes from its C++ layer); a
        # strict decoder would crash before we ever see the message. errors
        # only applies in text mode — pass it only when text=True; binary
        # mode reads raw bytes, decoded the same way below.
        text_kwargs: dict[str, object] = (
            {"text": True, "errors": "replace"} if not binary else {"text": False}
        )
        with _tracked_popen(
            cmd,
            cwd=cwd,
            env=_subprocess_env(extra_env),
            missing_tool=tool,
            stdin=stdin_target,
            stdout=stdout_target,
            stderr=subprocess.PIPE,
            **text_kwargs,
        ) as proc:
            _, stderr = proc.communicate()

    stderr_text = stderr if isinstance(stderr, str) else stderr.decode(errors="replace")
    if stderr_text:
        log.debug("stderr from %s:\n%s", tool, stderr_text)

    if proc.returncode != 0:
        log.error("Command failed (exit %d): %s", proc.returncode, cmd_str)
        raise ExternalToolError(
            f"{tool} failed (exit {proc.returncode})",
            tool=tool, returncode=proc.returncode, cmd=cmd,
            stderr_snippet=stderr_text,
        )


def run_piped(
    cmd1: list[str],
    cmd2: list[str],
    *,
    cwd: Path,
    out_file: str | None = None,
) -> str:
    """Run two commands connected by a pipe (cmd1 | cmd2).

    Returns:
        The stdout of cmd2 as a string.
    """
    cmd_str = f"{' '.join(cmd1)} | {' '.join(cmd2)}"
    log.debug("Running pipe: %s", cmd_str)

    env = _subprocess_env()
    # Drain the producer's stderr to a temp file rather than an unread PIPE: an
    # unread stderr PIPE deadlocks any producer chatty enough to fill the
    # ~64 KB buffer before the consumer drains its stdout. Capturing it also
    # lets us report *why* the producer died (e.g. an OOM kill) instead of
    # only the consumer's downstream symptom.
    with tempfile.TemporaryFile() as p1_err, _tracked_popen(
        cmd1, cwd=cwd, env=env, missing_tool=_tool_name(cmd1),
        stdout=subprocess.PIPE, stderr=p1_err,
    ) as p1:
        assert p1.stdout is not None  # PIPE was requested above
        try:
            with _tracked_popen(
                cmd2, cwd=cwd, env=env, missing_tool=_tool_name(cmd2),
                stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            ) as p2:
                p1.stdout.close()
                stdout, stderr = p2.communicate()
        except BaseException:
            # Tear down the producer if the consumer raises (or KeyboardInterrupt).
            _signal_group(p1, signal.SIGKILL)
            p1.wait()
            raise
        p1.wait()
        p1_err.seek(0)
        producer_err = p1_err.read().decode(errors="replace")

    output = stdout.decode(errors="replace")

    if out_file:
        with open(cwd / out_file, "w") as f:
            f.write(output)

    if p2.returncode != 0:
        stderr_text = stderr.decode(errors="replace")
        # Fold in the producer's fate so a crashed/OOM-killed producer feeding
        # an empty stream isn't masked by the consumer's "no input" error. Its
        # exit line goes after its stderr, so the kept tail of the snippet still
        # shows it. A producer killed by SIGPIPE only saw the consumer stop
        # reading: its stderr would just push the consumer's error out.
        if p1.returncode not in (0, None, -13, 141):
            stderr_text = (
                f"{producer_err.strip()}\n"
                f"{_tool_name(cmd1)} (producer) exited {p1.returncode}\n{stderr_text}"
            ).lstrip()
        raise ExternalToolError(
            f"{_tool_name(cmd2)} failed (exit {p2.returncode})",
            tool=_tool_name(cmd2),
            returncode=p2.returncode,
            cmd=cmd2,
            stderr_snippet=stderr_text,
        )

    # Consumer succeeded but the producer failed: its output was truncated, so
    # the pipeline result is unreliable. SIGPIPE (-13 / 141) is excluded — that
    # just means the consumer stopped reading early, which is legitimate.
    if p1.returncode not in (0, None, -13, 141):
        raise ExternalToolError(
            f"{_tool_name(cmd1)} failed (exit {p1.returncode})",
            tool=_tool_name(cmd1),
            returncode=p1.returncode,
            cmd=cmd1,
            stderr_snippet=producer_err,
        )

    return output

