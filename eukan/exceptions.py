"""Eukan exception hierarchy.

Structured, domain-specific exceptions carrying actionable context for
CLI error formatting. The CLI's top-level handler in ``cli.py``
dispatches on the category classes below, so the hierarchy stays
minimal and isinstance-friendly::

    EukanError
    ├── ConfigurationError
    │   └── InvalidOptionError
    ├── ValidationError
    │   ├── FastaValidationError
    │   └── GFFValidationError
    ├── DependencyError
    │   └── ToolEnvError
    └── ExternalToolError
"""

from __future__ import annotations

from pathlib import Path


class EukanError(Exception):
    """Base exception for all eukan errors.

    Args:
        message: Human-readable error description.
        hint: Optional remediation advice shown to the user.
    """

    def __init__(self, message: str, *, hint: str | None = None) -> None:
        self.hint = hint
        super().__init__(message)


class ConfigurationError(EukanError):
    """Bad settings, missing options, or invalid combinations."""


class InvalidOptionError(ConfigurationError):
    """An option value is out of range or an incompatible combination was given."""


class ValidationError(EukanError):
    """Invalid input files."""

    _kind: str = "input"

    def __init__(self, path: str | Path, message: str = "", **kwargs) -> None:
        self.path = Path(path)
        msg = f"Invalid {self._kind}: {self.path}"
        if message:
            msg = f"{msg} — {message}"
        super().__init__(msg, **kwargs)


class FastaValidationError(ValidationError):
    """Unparseable, empty, or malformed FASTA file."""

    _kind = "FASTA"


class GFFValidationError(ValidationError):
    """Structural problems in a GFF3 file."""

    _kind = "GFF3"


class DependencyError(EukanError):
    """Missing or broken external tools or databases."""


class ToolEnvError(DependencyError):
    """Required environment variable for a tool is not set."""

    def __init__(self, tool: str, *, env_var: str) -> None:
        self.tool = tool
        self.env_var = env_var
        super().__init__(
            f"{tool} requires environment variable {env_var} to be set",
            hint=f"Set {env_var} to the {tool} configuration directory.",
        )


class ExternalToolError(EukanError):
    """An external tool exited non-zero.

    Attributes:
        tool: Binary name (e.g. ``"gmes_petap.pl"``).
        returncode: Process exit code.
        cmd: Full command as a list of strings.
        stderr_snippet: First 500 characters of stderr output.
        step: Pipeline step name, if known.
    """

    def __init__(
        self,
        message: str,
        *,
        tool: str,
        returncode: int = 1,
        cmd: list[str] | None = None,
        stderr_snippet: str = "",
        step: str = "",
        **kwargs,
    ) -> None:
        self.tool = tool
        self.returncode = returncode
        self.cmd = cmd or []
        self.stderr_snippet = stderr_snippet[:500]
        self.step = step
        super().__init__(message, **kwargs)

    def __str__(self) -> str:
        msg = f"{self.tool} failed (exit {self.returncode})"
        if self.step:
            msg += f" during '{self.step}'"
        return msg
