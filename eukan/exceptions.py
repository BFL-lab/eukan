"""Eukan exception hierarchy.

Provides structured, domain-specific exceptions that carry actionable context
for CLI error formatting and programmatic handling.

Hierarchy::

    EukanError
    ├── ConfigurationError
    │   ├── MissingOptionError
    │   └── InvalidOptionError
    ├── ValidationError
    │   ├── FastaValidationError
    │   ├── GFFValidationError
    │   └── ReadsValidationError
    ├── DependencyError
    │   ├── ToolNotFoundError
    │   ├── ToolVersionError
    │   ├── ToolEnvError
    │   └── DatabaseError
    ├── PipelineStepError
    │   ├── AnnotationError
    │   │   ├── GeneMarkError
    │   │   ├── AugustusError
    │   │   ├── SnapError
    │   │   ├── CodingQuarryError
    │   │   ├── AlignmentError
    │   │   └── EVMError
    │   ├── AssemblyError
    │   │   ├── STARError
    │   │   ├── TrinityError
    │   │   └── PASAError
    │   ├── FunctionalAnnotationError
    │   └── ExternalToolError
    └── GFFError
        ├── MalformedGFFError
        └── GFFTransformError
"""

from __future__ import annotations

from pathlib import Path


# ---------------------------------------------------------------------------
# Base
# ---------------------------------------------------------------------------


class EukanError(Exception):
    """Base exception for all eukan errors.

    Args:
        message: Human-readable error description.
        hint: Optional remediation advice shown to the user.
    """

    def __init__(self, message: str, *, hint: str | None = None) -> None:
        self.hint = hint
        super().__init__(message)


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------


class ConfigurationError(EukanError):
    """Bad settings, missing options, or invalid combinations."""


class MissingOptionError(ConfigurationError):
    """A required CLI option or config field was not provided."""


class InvalidOptionError(ConfigurationError):
    """An option value is out of range or an incompatible combination was given."""


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


class ValidationError(EukanError):
    """Invalid input files."""


class FastaValidationError(ValidationError):
    """Unparseable, empty, or malformed FASTA file."""

    def __init__(self, path: str | Path, message: str = "", **kwargs) -> None:
        self.path = Path(path)
        msg = f"Invalid FASTA: {self.path}"
        if message:
            msg = f"{msg} — {message}"
        super().__init__(msg, **kwargs)


class GFFValidationError(ValidationError):
    """Structural problems in a GFF3 file."""

    def __init__(self, path: str | Path, message: str = "", **kwargs) -> None:
        self.path = Path(path)
        msg = f"Invalid GFF3: {self.path}"
        if message:
            msg = f"{msg} — {message}"
        super().__init__(msg, **kwargs)


class ReadsValidationError(ValidationError):
    """Missing, empty, or unreadable sequencing read files."""


# ---------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------


class DependencyError(EukanError):
    """Missing or broken external tools or databases."""


class ToolNotFoundError(DependencyError):
    """External tool not on PATH or not executable."""

    def __init__(self, tool: str, *, install_hint: str | None = None) -> None:
        self.tool = tool
        self.install_hint = install_hint
        msg = f"Tool not found: {tool}"
        super().__init__(msg, hint=install_hint)


class ToolVersionError(DependencyError):
    """Tool found but below the minimum required version."""

    def __init__(self, tool: str, *, found: str, required: str) -> None:
        self.tool = tool
        self.found_version = found
        self.required_version = required
        super().__init__(
            f"{tool}: found version {found}, requires >= {required}",
            hint=f"Upgrade {tool} to at least version {required}.",
        )


class ToolEnvError(DependencyError):
    """Required environment variable for a tool is not set."""

    def __init__(self, tool: str, *, env_var: str) -> None:
        self.tool = tool
        self.env_var = env_var
        super().__init__(
            f"{tool} requires environment variable {env_var} to be set",
            hint=f"Set {env_var} to the {tool} configuration directory.",
        )


class DatabaseError(DependencyError):
    """Missing, corrupt, or stale database files."""


# ---------------------------------------------------------------------------
# Pipeline steps
# ---------------------------------------------------------------------------


class PipelineStepError(EukanError):
    """A named pipeline step failed."""

    def __init__(self, message: str, *, step: str = "", **kwargs) -> None:
        self.step = step
        super().__init__(message, **kwargs)


class ExternalToolError(PipelineStepError):
    """An external tool exited non-zero.

    Attributes:
        tool: Binary name (e.g. ``"gmes_petap.pl"``).
        returncode: Process exit code.
        cmd: Full command as a list of strings.
        stderr_snippet: First 500 characters of stderr output.
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
        super().__init__(message, step=step, **kwargs)

    def __str__(self) -> str:
        msg = f"{self.tool} failed (exit {self.returncode})"
        if self.step:
            msg += f" during '{self.step}'"
        return msg


# ---------------------------------------------------------------------------
# Annotation pipeline errors
# ---------------------------------------------------------------------------


class AnnotationError(PipelineStepError):
    """Failure during the annotation pipeline."""


class GeneMarkError(AnnotationError, ExternalToolError):
    """GeneMark gene prediction failed."""


class AugustusError(AnnotationError, ExternalToolError):
    """AUGUSTUS training or prediction failed."""


class SnapError(AnnotationError, ExternalToolError):
    """SNAP training or prediction failed."""


class CodingQuarryError(AnnotationError, ExternalToolError):
    """CodingQuarry prediction failed."""


class AlignmentError(AnnotationError, ExternalToolError):
    """Protein alignment (spaln/gth) failed."""


class EVMError(AnnotationError, ExternalToolError):
    """EvidenceModeler consensus building failed."""


# ---------------------------------------------------------------------------
# Assembly pipeline errors
# ---------------------------------------------------------------------------


class AssemblyError(PipelineStepError):
    """Failure during the assembly pipeline."""


class STARError(AssemblyError, ExternalToolError):
    """STAR read mapping failed."""


class TrinityError(AssemblyError, ExternalToolError):
    """Trinity assembly failed."""


class PASAError(AssemblyError, ExternalToolError):
    """PASA spliced alignment failed."""


# ---------------------------------------------------------------------------
# Functional annotation pipeline errors
# ---------------------------------------------------------------------------


class FunctionalAnnotationError(PipelineStepError):
    """Failure during functional annotation."""


# ---------------------------------------------------------------------------
# GFF errors
# ---------------------------------------------------------------------------


class GFFError(EukanError):
    """GFF3 parsing or transformation failure."""


class MalformedGFFError(GFFError):
    """Structural problems in GFF3 input."""


class GFFTransformError(GFFError):
    """Failure during a gffutils transform pass."""
