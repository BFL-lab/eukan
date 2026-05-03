"""Safe subprocess execution for external bioinformatics tools."""

from __future__ import annotations

import contextlib
import subprocess
import threading
import time
from pathlib import Path

from eukan.exceptions import ExternalToolError
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


def terminate_all_children(grace_period: float = 5.0) -> None:
    """Terminate every tracked child process, then SIGKILL stragglers.

    Called from the CLI's SIGINT handler so Ctrl-C doesn't orphan running
    bioinformatics tools.
    """
    with _REGISTRY_LOCK:
        procs = list(_RUNNING)
    if not procs:
        return
    log.warning("Interrupted; terminating %d running child process(es)...", len(procs))
    for p in procs:
        with contextlib.suppress(OSError):
            p.terminate()
    deadline = time.monotonic() + grace_period
    for p in procs:
        remaining = max(0.0, deadline - time.monotonic())
        try:
            p.wait(timeout=remaining)
        except subprocess.TimeoutExpired:
            with contextlib.suppress(OSError):
                p.kill()


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
    err_file: str | None = None,
    binary: bool = False,
    timeout: int | None = None,
    extra_env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess:
    """Run an external command safely.

    Args:
        cmd: Command as a list of strings (never a shell string).
        cwd: Working directory for the subprocess.
        in_file: If set, the file at this path within *cwd* is opened
            for reading and connected to the child's stdin (no Python
            buffering — the child reads from the fd directly).
        out_file: If set, stdout is streamed straight to this filename
            within *cwd* via the child's fd (no in-process capture).
        err_file: If set, stderr is streamed straight to this filename
            within *cwd* via the child's fd.
        binary: Only meaningful when *out_file* is unset; controls
            whether captured stdout is returned as bytes or str.
        timeout: Optional timeout in seconds.
        extra_env: Additional environment variables merged on top of
            the current environment for this subprocess only.

    Returns:
        The CompletedProcess result.

    Raises:
        ExternalToolError: If the command exits non-zero or times out.
    """
    cmd_str = " ".join(cmd)
    log.debug("Running: %s (cwd=%s)", cmd_str, cwd)

    # File handles' lifetimes span the subprocess call; closed in finally.
    in_handle = open(cwd / in_file, "rb") if in_file else None  # noqa: SIM115
    out_handle = open(cwd / out_file, "wb") if out_file else None  # noqa: SIM115
    err_handle = open(cwd / err_file, "wb") if err_file else None  # noqa: SIM115

    stdin_target: int | None = in_handle.fileno() if in_handle else None
    stdout_target: int = out_handle.fileno() if out_handle else subprocess.DEVNULL
    stderr_target: int = err_handle.fileno() if err_handle else subprocess.PIPE

    # Run via Popen so we can register the child for SIGINT cleanup.
    # start_new_session=True isolates the child from terminal SIGINT, so
    # only our handler decides when to terminate it.
    proc = subprocess.Popen(
        cmd,
        cwd=cwd,
        env=_subprocess_env(extra_env),
        stdin=stdin_target,
        stdout=stdout_target,
        stderr=stderr_target,
        text=not binary,
        start_new_session=True,
    )
    _track(proc)
    try:
        try:
            _, stderr = proc.communicate(timeout=timeout)
        except subprocess.TimeoutExpired as exc:
            proc.kill()
            with contextlib.suppress(subprocess.TimeoutExpired):
                proc.communicate(timeout=5)
            raise ExternalToolError(
                f"{_tool_name(cmd)} timed out after {timeout}s",
                tool=_tool_name(cmd),
                returncode=-1,
                cmd=cmd,
                stderr_snippet=str(exc),
            ) from exc
        except KeyboardInterrupt:
            proc.terminate()
            try:
                proc.communicate(timeout=5)
            except subprocess.TimeoutExpired:
                proc.kill()
            raise
    finally:
        _untrack(proc)
        for h in (in_handle, out_handle, err_handle):
            if h is not None:
                h.close()

    if stderr is None:
        stderr_text = ""
    else:
        stderr_text = stderr if isinstance(stderr, str) else stderr.decode(errors="replace")
    if stderr_text:
        log.debug("stderr from %s:\n%s", _tool_name(cmd), stderr_text)

    if proc.returncode != 0:
        log.error("Command failed (exit %d): %s", proc.returncode, cmd_str)
        raise ExternalToolError(
            f"{_tool_name(cmd)} failed (exit {proc.returncode})",
            tool=_tool_name(cmd),
            returncode=proc.returncode,
            cmd=cmd,
            stderr_snippet=stderr_text,
        )

    return subprocess.CompletedProcess(
        args=cmd, returncode=proc.returncode, stdout=None, stderr=stderr,
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
    p1 = subprocess.Popen(
        cmd1, cwd=cwd, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        start_new_session=True,
    )
    assert p1.stdout is not None  # PIPE was requested above
    _track(p1)
    p2: subprocess.Popen | None = None
    try:
        p2 = subprocess.Popen(
            cmd2, cwd=cwd, env=env, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            start_new_session=True,
        )
        _track(p2)
        p1.stdout.close()

        stdout, stderr = p2.communicate()
    except (Exception, KeyboardInterrupt):
        p1.kill()
        if p2 is not None:
            p2.kill()
        p1.wait()
        if p2 is not None:
            p2.wait()
        raise
    finally:
        p1.wait()
        _untrack(p1)
        if p2 is not None:
            _untrack(p2)

    output = stdout.decode()

    if out_file:
        with open(cwd / out_file, "w") as f:
            f.write(output)

    if p2.returncode != 0:
        stderr_text = stderr.decode(errors="replace")
        raise ExternalToolError(
            f"{_tool_name(cmd2)} failed (exit {p2.returncode})",
            tool=_tool_name(cmd2),
            returncode=p2.returncode,
            cmd=cmd2,
            stderr_snippet=stderr_text,
        )

    return output


def run_shell(
    cmd_str: str,
    *,
    cwd: Path,
) -> subprocess.CompletedProcess:
    """Run a shell command string with proper error handling.

    Use this only when the command requires shell features (redirections,
    pipes built by external tools, etc.) that cannot be expressed as a
    list.  Prefer :func:`run_cmd` wherever possible.

    Raises:
        ExternalToolError: If the command exits non-zero.
    """
    log.debug("Running (shell): %s (cwd=%s)", cmd_str, cwd)

    proc = subprocess.Popen(
        cmd_str, shell=True, cwd=cwd, env=_subprocess_env(),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        start_new_session=True,
    )
    _track(proc)
    try:
        try:
            stdout, stderr = proc.communicate()
        except KeyboardInterrupt:
            proc.terminate()
            try:
                proc.communicate(timeout=5)
            except subprocess.TimeoutExpired:
                proc.kill()
            raise
    finally:
        _untrack(proc)

    if stderr:
        log.debug("stderr (shell):\n%s", stderr)

    if proc.returncode != 0:
        tool = cmd_str.split()[0] if cmd_str.strip() else "shell"
        log.error("Shell command failed (exit %d): %s", proc.returncode, cmd_str[:200])
        raise ExternalToolError(
            f"Shell command failed (exit {proc.returncode})",
            tool=tool,
            returncode=proc.returncode,
            cmd=[cmd_str],
            stderr_snippet=stderr,
        )

    return subprocess.CompletedProcess(
        args=cmd_str, returncode=proc.returncode, stdout=stdout, stderr=stderr,
    )
