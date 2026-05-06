"""Safe subprocess execution for external bioinformatics tools."""

from __future__ import annotations

import contextlib
import subprocess
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


@contextmanager
def _tracked_popen(
    cmd: list[str] | str,
    *,
    cwd: Path,
    env: dict[str, str] | None,
    missing_tool: str | None = None,
    **popen_kwargs,
) -> Iterator[subprocess.Popen]:
    """Spawn a Popen registered for SIGINT cleanup, with isolated session.

    On ``KeyboardInterrupt`` the child is terminated (5s grace) then
    SIGKILLed; in either case the process is unregistered before the
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
        proc.terminate()
        with contextlib.suppress(subprocess.TimeoutExpired):
            proc.communicate(timeout=5)
        if proc.returncode is None:
            proc.kill()
        raise
    finally:
        _untrack(proc)


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

    tool = _tool_name(cmd)
    try:
        with _tracked_popen(
            cmd,
            cwd=cwd,
            env=_subprocess_env(extra_env),
            missing_tool=tool,
            stdin=stdin_target,
            stdout=stdout_target,
            stderr=stderr_target,
            text=not binary,
        ) as proc:
            try:
                _, stderr = proc.communicate(timeout=timeout)
            except subprocess.TimeoutExpired as exc:
                proc.kill()
                with contextlib.suppress(subprocess.TimeoutExpired):
                    proc.communicate(timeout=5)
                raise ExternalToolError(
                    f"{tool} timed out after {timeout}s",
                    tool=tool, returncode=-1, cmd=cmd,
                    stderr_snippet=str(exc),
                ) from exc
    finally:
        for h in (in_handle, out_handle, err_handle):
            if h is not None:
                h.close()

    if stderr is None:
        stderr_text = ""
    else:
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
    with _tracked_popen(
        cmd1, cwd=cwd, env=env, missing_tool=_tool_name(cmd1),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
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
            p1.kill()
            p1.wait()
            raise
        p1.wait()

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

    with _tracked_popen(
        cmd_str, cwd=cwd, env=_subprocess_env(),
        shell=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    ) as proc:
        stdout, stderr = proc.communicate()

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
