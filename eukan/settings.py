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

from eukan.infra.logging import get_logger

from pydantic import Field, computed_field, model_validator
from pydantic_settings import BaseSettings, SettingsConfigDict


def _rand_string(length: int = 5) -> str:
    return "".join(random.choice(string.ascii_uppercase) for _ in range(length))


def _pyproject_toml_settings(settings: BaseSettings) -> dict[str, Any]:
    """Load [tool.eukan] from pyproject.toml if it exists."""
    pyproject = Path("pyproject.toml")
    if not pyproject.exists():
        return {}
    try:
        from eukan.infra import tomllib
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
    @classmethod
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

    return settings_customise_sources


class Kingdom(str, Enum):
    fungus = "fungus"
    protist = "protist"
    animal = "animal"
    plant = "plant"


# ---------------------------------------------------------------------------
# Pipeline settings (eukan annotate)
# ---------------------------------------------------------------------------


class PipelineConfig(BaseSettings):
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
    work_dir: Path = Field(default_factory=Path.cwd)
    manifest_dir: Path | None = None  # defaults to work_dir if not set
    name: str = ""  # derived from genome stem if not set
    shortname: str = Field(default_factory=_rand_string)
    kingdom: Kingdom | None = None
    num_cpu: int = Field(default_factory=lambda: os.cpu_count() or 1)
    genetic_code: str = "11"
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
    def _default_manifest_dir(self) -> "PipelineConfig":
        if not self.manifest_dir:
            object.__setattr__(self, "manifest_dir", self.work_dir)
        return self

    @model_validator(mode="after")
    def _derive_name(self) -> "PipelineConfig":
        if not self.name:
            object.__setattr__(self, "name", self.genome.stem)
        return self

    @model_validator(mode="after")
    def _discover_assembly_outputs(self) -> "PipelineConfig":
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

    @cached_property
    def genetic_code_obj(self):
        """Return a :class:`~eukan.gencode.GeneticCode` for this config's code."""
        from eukan.gencode import GeneticCode
        return GeneticCode(self.genetic_code)

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


class AssemblyConfig(BaseSettings):
    """Configuration for transcriptome assembly."""

    model_config = SettingsConfigDict(
        env_prefix="EUKAN_ASSEMBLE_",
        extra="ignore",
    )

    genome: Path
    work_dir: Path = Field(default_factory=Path.cwd)
    manifest_dir: Path | None = None  # defaults to work_dir if not set
    left_reads: Path | None = None
    right_reads: Path | None = None
    single_reads: Path | None = None
    min_intron_len: int = 20
    max_intron_len: int = 5000
    phred_quality: int = 33
    num_cpu: int = Field(default_factory=lambda: os.cpu_count() or 1)
    strand_specific: str | None = None
    align_mode: str = "Local"
    jaccard_clip: bool = False
    genetic_code: str = "1"
    splice_permissive: bool = False

    @model_validator(mode="after")
    def _default_manifest_dir(self) -> "AssemblyConfig":
        if not self.manifest_dir:
            object.__setattr__(self, "manifest_dir", self.work_dir)
        return self

    @cached_property
    def genetic_code_obj(self):
        """Return a :class:`~eukan.gencode.GeneticCode` for this config's code."""
        from eukan.gencode import GeneticCode
        return GeneticCode(self.genetic_code)

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

    @computed_field  # type: ignore[prop-decorator]
    @property
    def memory_gb(self) -> str:
        try:
            with open("/proc/meminfo") as f:
                for line in f:
                    if line.startswith("MemTotal:"):
                        kb = int(line.split()[1])
                        return f"{kb // (2 * 1024 * 1024)}G"
        except (FileNotFoundError, ValueError):
            pass
        return "4G"

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
