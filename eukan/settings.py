"""Centralized configuration via pydantic-settings.

Settings are resolved in this order (last wins):
  1. Defaults defined here
  2. [tool.eukan] section in pyproject.toml
  3. Environment variables prefixed with EUKAN_
  4. CLI flags (applied by cli.py when constructing models)
"""

from __future__ import annotations

import os
import random
import string
from enum import Enum
from functools import cached_property
from pathlib import Path
from typing import Any

from pydantic import Field, computed_field, model_validator
from pydantic_settings import BaseSettings, SettingsConfigDict

from eukan.infra.logging import get_logger


def _rand_string(length: int = 5) -> str:
    return "".join(random.choice(string.ascii_uppercase) for _ in range(length))


def _pyproject_toml_settings(settings: BaseSettings) -> dict[str, Any]:
    """Load [tool.eukan] from pyproject.toml if it exists."""
    import tomllib
    pyproject = Path("pyproject.toml")
    if not pyproject.exists():
        return {}
    try:
        with open(pyproject, "rb") as f:
            data = tomllib.load(f)
        return data.get("tool", {}).get("eukan", {})
    except (OSError, ValueError, KeyError):
        return {}


def _pyproject_settings_sources(toml_key: str | None = None):
    """Create a settings_customise_sources classmethod for pyproject.toml.

    Args:
        toml_key: Sub-key under [tool.eukan] (e.g., "assemble").
            If None, reads directly from [tool.eukan].
    """
    def settings_customise_sources(cls, settings_cls, **kwargs):
        from pydantic_settings import PydanticBaseSettingsSource

        class PyprojectSource(PydanticBaseSettingsSource):
            def get_field_value(self, field, field_name):
                data = _pyproject_toml_settings(self.settings_cls)
                section = data.get(toml_key, {}) if toml_key else data
                val = section.get(field_name)
                return val, field_name, val is not None

            def __call__(self):
                data = _pyproject_toml_settings(self.settings_cls)
                return data.get(toml_key, {}) if toml_key else data

        return (
            kwargs.get("init_settings"),
            kwargs.get("env_settings"),
            PyprojectSource(settings_cls),
        )

    return classmethod(settings_customise_sources)


class Kingdom(str, Enum):
    fungus = "fungus"
    protist = "protist"
    animal = "animal"
    plant = "plant"


# ---------------------------------------------------------------------------
# Shared base for configs that drive multi-step jobs out of a work_dir
# ---------------------------------------------------------------------------


class _StepRunSettings(BaseSettings):
    """Common fields, validators, and accessors shared by step-driven configs.

    Subclasses (PipelineConfig, AssemblyConfig) provide their own
    ``model_config`` (env_prefix, settings sources) and may override
    field defaults (e.g. ``genetic_code``).
    """

    work_dir: Path = Field(default_factory=Path.cwd)
    # Filled in from work_dir by _default_manifest_dir if not provided;
    # always a Path after construction.
    manifest_dir: Path = Field(default_factory=Path.cwd)
    num_cpu: int = Field(default_factory=lambda: os.cpu_count() or 1)
    genetic_code: str = "1"

    @model_validator(mode="before")
    @classmethod
    def _default_manifest_dir(cls, data: Any) -> Any:
        # Run before field validation so manifest_dir can be a non-Optional
        # Path -- mypy then sees it as guaranteed throughout the codebase.
        if isinstance(data, dict) and not data.get("manifest_dir"):
            data["manifest_dir"] = data.get("work_dir") or Path.cwd()
        return data

    @cached_property
    def genetic_code_obj(self):
        """Return a :class:`~eukan.gencode.GeneticCode` for this config's code."""
        from eukan.gencode import GeneticCode
        return GeneticCode(self.genetic_code)


# ---------------------------------------------------------------------------
# Pipeline settings (eukan annotate)
# ---------------------------------------------------------------------------


class PipelineConfig(_StepRunSettings):
    """Configuration for the annotation pipeline.

    Fields can be set via:
      - [tool.eukan] in pyproject.toml
      - EUKAN_ prefixed env vars (e.g., EUKAN_NUM_CPU=8)
      - CLI flags (override at construction time)

    Field organization below:
      1. Required raw inputs
      2. Defaulted raw inputs
      3. Optional assembly-evidence paths (auto-discovered if not set)
      4. Validators (fill in derived defaults)
      5. Computed properties (``is_fungus``, ``has_transcripts``,
         ``genetic_code_obj``) — derived on access, never stored.
    """

    model_config = SettingsConfigDict(
        env_prefix="EUKAN_",
        env_nested_delimiter="__",
        extra="ignore",
    )

    # --- Required (must come from CLI) ---
    genome: Path
    proteins: list[Path]

    # --- Defaulted (overridable via config/env/CLI) ---
    name: str = ""  # derived from genome stem if not set
    shortname: str = Field(default_factory=_rand_string)
    kingdom: Kingdom | None = None
    genetic_code: str = "11"  # override base default
    weights: list[int] = Field(default_factory=lambda: [2, 1, 3])
    spaln_ssp: bool = False
    allow_noncanonical_splice: bool = False

    # --- Optional transcript evidence (auto-discovered from work_dir if not set) ---
    transcripts_fasta: Path | None = None
    transcripts_gff: Path | None = None
    rnaseq_hints: Path | None = None
    strand_specific: bool = False
    utrs_db: Path | None = None

    # Well-known assembly output filenames
    _ASSEMBLY_FILES: dict[str, str] = {
        "transcripts_fasta": "nr_transcripts.fasta",
        "transcripts_gff": "nr_transcripts.gff3",
        "rnaseq_hints": "hints_rnaseq.gff",
    }

    # --- Validators ------------------------------------------------------

    @model_validator(mode="after")
    def _derive_name(self) -> PipelineConfig:
        if not self.name:
            object.__setattr__(self, "name", self.genome.stem)
        return self

    @model_validator(mode="after")
    def _discover_assembly_outputs(self) -> PipelineConfig:
        """Auto-discover assembly outputs in work_dir when not explicitly provided."""
        log = get_logger(__name__)

        # If the user already set all three explicitly, nothing to do
        if all([self.transcripts_fasta, self.transcripts_gff, self.rnaseq_hints]):
            return self

        # If the user set some but not all explicitly, don't override their intent
        explicitly_set = {
            field: getattr(self, field)
            for field in self._ASSEMBLY_FILES
            if getattr(self, field) is not None
        }
        if explicitly_set:
            return self

        # Scan work_dir for assembly outputs
        found: dict[str, Path] = {}
        missing: list[str] = []
        for field, filename in self._ASSEMBLY_FILES.items():
            path = self.work_dir / filename
            if path.exists():
                found[field] = path
            else:
                missing.append(filename)

        if found and missing:
            log.warning(
                "Partial assembly outputs in %s: found %s but missing %s. "
                "Run `eukan assemble` to completion or remove partial files.",
                self.work_dir,
                ", ".join(p.name for p in found.values()),
                ", ".join(missing),
            )
        elif found:
            log.info(
                "Auto-discovered assembly outputs in %s", self.work_dir,
            )
            for field, path in found.items():
                object.__setattr__(self, field, path)

        return self

    # --- Computed properties --------------------------------------------

    @computed_field  # type: ignore[prop-decorator]
    @property
    def is_fungus(self) -> bool:
        return self.kingdom == Kingdom.fungus

    @computed_field  # type: ignore[prop-decorator]
    @property
    def is_protist(self) -> bool:
        return self.kingdom == Kingdom.protist

    @computed_field  # type: ignore[prop-decorator]
    @property
    def has_transcripts(self) -> bool:
        return all([self.transcripts_fasta, self.transcripts_gff, self.rnaseq_hints])

    settings_customise_sources = _pyproject_settings_sources()


# ---------------------------------------------------------------------------
# Assembly settings (eukan assemble)
# ---------------------------------------------------------------------------


def _default_trinity_memory_gb(meminfo_path: str = "/proc/meminfo") -> int:
    """Safe default for Trinity ``--max_memory``, in GiB.

    Trinity (Jellyfish in genome-guided mode, Inchworm in de novo) reliably
    overshoots its ``--max_memory`` soft cap during k-mer counting. We size
    the cap from ``MemAvailable`` (the kernel's estimate of memory free for
    new processes) rather than ``MemTotal``, so the cap reflects what the
    machine can actually spare. Falls back to half of ``MemTotal``, then to
    4 GiB. Always at least 4 GiB — Trinity needs that much to run at all.
    """
    try:
        avail_kb = total_kb = 0
        with open(meminfo_path) as f:
            for line in f:
                if line.startswith("MemAvailable:"):
                    avail_kb = int(line.split()[1])
                elif line.startswith("MemTotal:"):
                    total_kb = int(line.split()[1])
        if avail_kb > 0:
            return max(4, int(avail_kb * 0.6 / (1024 * 1024)))
        if total_kb > 0:
            return max(4, total_kb // (2 * 1024 * 1024))
    except (OSError, ValueError):
        pass
    return 4


class AssemblyConfig(_StepRunSettings):
    """Configuration for transcriptome assembly."""

    model_config = SettingsConfigDict(
        env_prefix="EUKAN_ASSEMBLE_",
        extra="ignore",
    )

    genome: Path
    left_reads: Path | None = None
    right_reads: Path | None = None
    single_reads: Path | None = None
    min_intron_len: int = 20
    max_intron_len: int = 5000
    phred_quality: int = 33
    strand_specific: str | None = None
    align_mode: str = "Local"
    jaccard_clip: bool = False
    splice_permissive: bool = False

    @computed_field  # type: ignore[prop-decorator]
    @property
    def name(self) -> str:
        return self.genome.stem

    @computed_field  # type: ignore[prop-decorator]
    @property
    def reads_args_star(self) -> list[str]:
        if self.left_reads and self.right_reads:
            return [str(self.left_reads), str(self.right_reads)]
        elif self.single_reads:
            return [str(self.single_reads)]
        raise ValueError("No read files provided")

    @computed_field  # type: ignore[prop-decorator]
    @property
    def reads_args_trinity(self) -> list[str]:
        if self.left_reads and self.right_reads:
            return ["--left", str(self.left_reads), "--right", str(self.right_reads)]
        elif self.single_reads:
            return ["--single", str(self.single_reads)]
        raise ValueError("No read files provided")

    memory_gb: int = Field(
        default_factory=lambda: _default_trinity_memory_gb(),
        description="Trinity --max_memory cap in GiB.",
    )

    settings_customise_sources = _pyproject_settings_sources("assemble")


# ---------------------------------------------------------------------------
# Functional annotation settings (eukan func-annot)
# ---------------------------------------------------------------------------


class FunctionalConfig(BaseSettings):
    """Configuration for functional annotation."""

    model_config = SettingsConfigDict(
        env_prefix="EUKAN_FUNC_",
        extra="ignore",
    )

    proteins: Path
    uniprot_db: Path = Path("databases/uniprot_sprot.faa")
    pfam_db: Path = Path("databases/Pfam-A.hmm")
    gff3_path: Path | None = None
    num_cpu: int = Field(default_factory=lambda: os.cpu_count() or 1)
    evalue: str = "1e-1"

    settings_customise_sources = _pyproject_settings_sources("func-annot")
