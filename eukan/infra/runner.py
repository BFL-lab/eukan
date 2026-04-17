"""Safe subprocess execution for external bioinformatics tools."""

from __future__ import annotations

import subprocess
from pathlib import Path

from eukan.exceptions import ExternalToolError
from eukan.infra.environ import subprocess_env as _subprocess_env
from eukan.infra.logging import get_logger

log = get_logger(__name__)


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
    out_file: str | None = None,
    binary: bool = False,
    timeout: int | None = None,
    extra_env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess:
    """Run an external command safely.

    Args:
        cmd: Command as a list of strings (never a shell string).
        cwd: Working directory for the subprocess.
        out_file: If set, stdout is written to this filename within *cwd*.
        binary: If True, capture stdout as bytes (for BAM, gzip, etc.).
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

    stdout_path = cwd / out_file if out_file else None

    try:
        result = subprocess.run(
            cmd,
            cwd=cwd,
            env=_subprocess_env(extra_env),
            stdout=subprocess.PIPE if stdout_path else subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=not binary,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired as exc:
        raise ExternalToolError(
            f"{_tool_name(cmd)} timed out after {timeout}s",
            tool=_tool_name(cmd),
            returncode=-1,
            cmd=cmd,
            stderr_snippet=str(exc),
        ) from exc

    if stdout_path and result.stdout:
        mode = "wb" if binary else "w"
        with open(stdout_path, mode) as f:
            f.write(result.stdout)

    stderr_text = (
        result.stderr
        if isinstance(result.stderr, str)
        else result.stderr.decode(errors="replace")
    )
    if stderr_text:
        log.debug("stderr from %s:\n%s", _tool_name(cmd), stderr_text)

    if result.returncode != 0:
        log.error("Command failed (exit %d): %s", result.returncode, cmd_str)
        raise ExternalToolError(
            f"{_tool_name(cmd)} failed (exit {result.returncode})",
            tool=_tool_name(cmd),
            returncode=result.returncode,
            cmd=cmd,
            stderr_snippet=stderr_text,
        )

    return result


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
        cmd1, cwd=cwd, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    try:
        p2 = subprocess.Popen(
            cmd2, cwd=cwd, env=env, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        p1.stdout.close()

        stdout, stderr = p2.communicate()
    except Exception:
        p1.kill()
        p1.wait()
        raise
    finally:
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

    result = subprocess.run(
        cmd_str, shell=True, cwd=cwd, env=_subprocess_env(), capture_output=True, text=True,
    )

    if result.stderr:
        log.debug("stderr (shell):\n%s", result.stderr)

    if result.returncode != 0:
        # Try to extract tool name from the command string
        tool = cmd_str.split()[0] if cmd_str.strip() else "shell"
        log.error("Shell command failed (exit %d): %s", result.returncode, cmd_str[:200])
        raise ExternalToolError(
            f"Shell command failed (exit {result.returncode})",
            tool=tool,
            returncode=result.returncode,
            cmd=[cmd_str],
            stderr_snippet=result.stderr,
        )

    return result
