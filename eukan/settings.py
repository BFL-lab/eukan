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
from typing import Any, Literal

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

    @model_validator(mode="after")
    def _ensure_dirs_exist(self) -> _StepRunSettings:
        # Each step's CLI now points work_dir at a subdir under cwd
        # (e.g. ./annotate/, ./repeats/) which may not exist yet.
        self.work_dir.mkdir(parents=True, exist_ok=True)
        self.manifest_dir.mkdir(parents=True, exist_ok=True)
        return self

    @cached_property
    def genetic_code_obj(self):
        """Return a :class:`~eukan.infra.genetic_code.GeneticCode` for this config's code."""
        from eukan.infra.genetic_code import GeneticCode
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
    combinr_path: Path | None = None
    """Explicit path to the combinr binary; resolved from PATH when unset."""
    combinr_stringent_overlap: float = 0.0
    """combinr ``--stringent-overlap`` for the ``consensus --alt-splice`` isoform
    grouping: two transcript isoforms attach to one gene only when their genomic spans
    overlap >= this percent of the shorter span. 0 = off (any overlap groups them).
    Raise it (e.g. 30) to stop short or tip-overlapping transcripts from welding
    collinear neighbours into one multi-isoform gene in dense / trans-spliced genomes.
    Does not affect the consensus region partitioner."""

    # --- Optional transcript evidence (auto-discovered from work_dir if not set) ---
    transcripts_fasta: Path | None = None
    transcripts_gff: Path | None = None
    rnaseq_hints: Path | None = None
    strand_specific: bool = False

    # Well-known assembly output filenames — sourced from the cross-pipeline
    # artifact registry so renaming an artifact only edits one file.
    @staticmethod
    def _assembly_files() -> dict[str, str]:
        from eukan.infra.artifacts import Artifact
        return {
            "transcripts_fasta": Artifact.NR_TRANSCRIPTS_FASTA.value,
            "transcripts_gff": Artifact.NR_TRANSCRIPTS_GFF.value,
            "rnaseq_hints": Artifact.RNASEQ_HINTS.value,
        }

    # --- Validators ------------------------------------------------------

    @model_validator(mode="after")
    def _derive_name(self) -> PipelineConfig:
        if not self.name:
            object.__setattr__(self, "name", self.genome.stem)
        return self

    @model_validator(mode="after")
    def _discover_assembly_outputs(self) -> PipelineConfig:
        """Auto-discover assembly outputs (own work_dir or sibling assemble/)."""
        log = get_logger(__name__)
        from eukan.infra import artifacts
        from eukan.infra.artifacts import Artifact

        # If the user already set all three explicitly, nothing to do
        if all([self.transcripts_fasta, self.transcripts_gff, self.rnaseq_hints]):
            return self

        # If the user set some but not all explicitly, don't override their intent
        field_to_artifact = {
            "transcripts_fasta": Artifact.NR_TRANSCRIPTS_FASTA,
            "transcripts_gff":   Artifact.NR_TRANSCRIPTS_GFF,
            "rnaseq_hints":      Artifact.RNASEQ_HINTS,
        }
        explicitly_set = {
            field: getattr(self, field)
            for field in field_to_artifact
            if getattr(self, field) is not None
        }
        if explicitly_set:
            return self

        # Search work_dir and (when in step layout) sibling assemble/
        found: dict[str, Path] = {}
        missing: list[str] = []
        for field, art in field_to_artifact.items():
            path = artifacts.find(self.work_dir, art)
            if path is not None:
                found[field] = path
            else:
                missing.append(art.value)

        if found and missing:
            log.warning(
                "Partial assembly outputs: found %s but missing %s. "
                "Run `eukan assemble` to completion or remove partial files.",
                ", ".join(str(p) for p in found.values()),
                ", ".join(missing),
            )
        elif found:
            log.info(
                "Auto-discovered assembly outputs: %s",
                ", ".join(str(p) for p in found.values()),
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


def _default_assembly_memory_gb(meminfo_path: str = "/proc/meminfo") -> int:
    """Safe default for the de novo assembler's memory cap, in GiB.

    rnaSPAdes sizes its ``-m`` budget from this. We derive it from
    ``MemAvailable`` (the kernel's estimate of memory free for new processes)
    rather than ``MemTotal``, so the cap reflects what the machine can actually
    spare. Falls back to half of ``MemTotal``, then to a 4 GiB floor.
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
    non_canonical: Literal["auto", "force", "off"] = "auto"
    """Non-canonical splice mapping control. ``auto`` (default) layers the
    non-canonical minimap2 flags (``-J 0 -C 3 --splice-flank=no``) on top of the
    splice preset only when the soft-clip diagnostic calls non-canonical splicing
    ``EXTENSIVE``; ``force`` always applies them; ``off`` never does."""
    jaccard_clip: bool = True
    """Run the in-house jaccard clipping step over the Trinity transcript FASTAs
    (both de novo and genome-guided), splitting two adjacent loci fused into one
    contig. On by default; a no-op on single-end reads (read-pair bridging is the
    clip signal). This replaces Trinity's own ``--jaccard_clip`` (the standalone
    step is faster — STAR, not bowtie2 — and tunable via the ``jaccard_*`` knobs)."""
    jaccard_greediness: float = 1.5
    """Coverage-adaptive slack for jaccard fusion-trough detection. At low read-pair
    depth the pseudocount keeps a real fusion junction's jaccard above the fixed
    trough floor, so the contig is never split; this widens the trough gate toward
    the deepest jaccard physically reachable at the local depth (times this factor),
    making low-coverage fusions splittable while leaving high-coverage behaviour
    unchanged. 0 disables the adaptation (Trinity-faithful fixed floor)."""
    jaccard_max_trough: float = 0.05
    """Jaccard trough-depth gate (Trinity's 0.05): a candidate fusion junction's
    per-position jaccard must dip to <= this for the contig to be cut there. LOWER
    it for more stringent clipping (a deeper, cleaner bridging trough required, fewer
    splits); raise it to clip on shallower dips. At low read-pair depth this floor is
    adapted upward by ``jaccard_greediness`` (bounded by ``jaccard_max_adaptive_trough``)."""
    jaccard_min_delta: float = 0.35
    """Jaccard flanking-hill requirement (Trinity's 0.35): the jaccard must rise at
    least this far above the trough on BOTH sides within the scan window for the dip
    to be called a fusion junction. RAISE it for more stringent clipping (a sharper
    hill-trough-hill shape demanded); lower it to accept fainter junctions."""
    jaccard_max_adaptive_trough: float = 0.30
    """Ceiling on the coverage-adaptive trough gate: ``jaccard_greediness`` only relaxes
    ``jaccard_max_trough`` up to this value, so even at very low read-pair depth a dip
    shallower than this is never treated as a fusion. LOWER it to keep low-coverage
    clipping stringent; raise it to allow splitting on fainter low-coverage troughs."""
    diagnose_softclips: bool = True

    # --- Genome-guided assembly (StringTie) stringency ---
    # DORMANT: StringTie is no longer in the active pipeline (Trinity covers both
    # de novo and genome-guided). These fields configure the kept-but-unwired
    # stringtie module and have no effect unless it is re-wired into _steps_for.
    stringtie_min_coverage: float = 1.5
    """StringTie ``-c``: minimum per-bp read coverage for a transcript to be
    assembled. Raised above StringTie's default of 1 to suppress low-coverage
    spurious models."""
    stringtie_min_isoform_fraction: float = 0.1
    """StringTie ``-f``: minimum isoform abundance as a fraction of a locus's
    dominant isoform. Raised above StringTie's default of 0.01 to drop minor
    noise isoforms that inflate the genome-guided set."""
    stringtie_min_junction_coverage: float = 1.0
    """StringTie ``-j``: minimum number of spliced reads spanning a junction for it
    to be kept. Left at StringTie's default of 1, but exposed so it can be raised:
    a single spurious junction read (e.g. from noisy mapping in dense or
    trans-spliced genomes) otherwise becomes a splice-graph edge and inflates
    isoforms / fuses neighbouring loci."""

    # --- de novo + combinr consolidation routine ---
    rnaspades: bool = True
    """DORMANT: rnaSPAdes is no longer in the active pipeline (Trinity covers de
    novo + genome-guided). Configures the kept-but-unwired rnaspades module; no
    effect unless it is re-wired into _steps_for."""
    min_sl_fragment: int = 25
    """Minimum length (nt) of a fragment kept after in-silico SL trans-splicing."""
    sl_sequence: str | None = None
    """Override the spliced-leader sequence used for SL detection (else taken from
    the read soft-clip verdict, or the dominant de novo insertion motif)."""
    combinr_path: Path | None = None
    """Explicit path to the combinr binary; resolved from PATH when unset."""
    combinr_stringent_overlap: float = 0.0
    """combinr ``--stringent-overlap`` for ``combinr assemble`` clustering: two
    transcripts cluster only when their genomic spans overlap >= this percent of the
    shorter span. 0 = off (any overlap clusters). Raise it (e.g. 30) to stop short or
    tip-overlapping transcripts pulling collinear neighbours into one cluster in dense
    / trans-spliced genomes."""
    sl_cluster_window: int = 5
    """Genomic window (bp) for consolidating SL acceptor sites per (chrom, strand)."""
    min_sl_clip_len: int = 8
    """Minimum soft-clip length (bp) considered as a spliced-leader acceptor signal."""
    min_sl_insertion_len: int = 10
    """Minimum internal-insertion length (bp) considered as a spliced-leader signal."""
    adapter_sequences: list[str] = Field(default_factory=list)
    """Extra sequencing-adapter sequences, added to the built-in Illumina/Nextera
    set (``eukan.assembly.sl_depletion._BUILTIN_ADAPTERS``), excluded from SL
    detection so residual adapter read-through isn't mistaken for a spliced leader."""
    sl_adapter_filter: bool = True
    """Screen recovered SL consensus / acceptor candidates against known adapters
    (the built-in set plus ``adapter_sequences``). Turn off to let adapter sequence
    be treated as a spliced leader — the escape hatch for the rare case where a
    genuine SL legitimately overlaps an adapter seed."""

    # --- Homology-based splice-strand correction (unstranded libraries) ---
    uniprot_db: Path | None = None
    """SwissProt FASTA (``uniprot_sprot.faa``) or prebuilt ``.dmnd`` enabling
    forward-frame ``diamond blastx`` strand correction. Unset (or with ``-S``)
    disables the strand_correct step."""
    min_strand_consensus: int = 50
    """Minimum confirmed introns required to trust the learned splice consensus;
    below this the correction falls back to canonical GT-AG."""
    strand_blastx_evalue: float = 1e-5
    """diamond blastx e-value cutoff for a transcript to count as homology-confirmed."""

    # --- Homology-grounded de-fusion (chimeric transcript splitting) ---
    defuse: bool = False
    """Split fused transcripts using protein homology: an ``--ultra-sensitive``
    ``diamond blastx`` vs SwissProt that finds >=2 distinct, non-overlapping protein
    hits on one transcript flags a chimera of two genes and cuts it at the inter-hit
    gap. Requires ``--uniprot`` (shares its DB); off by default."""
    defuse_overlap_tolerance: float = 0.10
    """Max fractional query overlap (of the shorter hit) for two protein hits to still
    count as *distinct, non-overlapping* evidence of separate genes (default 0.10)."""
    defuse_blastx_evalue: float = 1e-5
    """diamond blastx e-value cutoff for a protein hit to count toward de-fusion."""

    @computed_field  # type: ignore[prop-decorator]
    @property
    def name(self) -> str:
        return self.genome.stem

    @computed_field  # type: ignore[prop-decorator]
    @property
    def aligner_bam(self) -> str:
        """The read BAM downstream genome-guided steps (Trinity gg, SL read-side) use.

        minimap2 is the sole aligner, so this is a single constant; the escalation
        to non-canonical mapping overwrites this same BAM in place.
        """
        return "minimap2_Aligned.sortedByCoord.out.bam"

    @computed_field  # type: ignore[prop-decorator]
    @property
    def reads_args_minimap2(self) -> list[str]:
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
        default_factory=lambda: _default_assembly_memory_gb(),
        description="Assembly memory cap in GiB (Trinity --max_memory / rnaSPAdes -m).",
    )

    settings_customise_sources = _pyproject_settings_sources("assemble")


# ---------------------------------------------------------------------------
# Repeat-masking settings (eukan mask-repeats)
# ---------------------------------------------------------------------------


class RepeatsConfig(_StepRunSettings):
    """Configuration for the repeat-masking pipeline (eukan mask-repeats)."""

    model_config = SettingsConfigDict(
        env_prefix="EUKAN_REPEATS_",
        extra="ignore",
    )

    genome: Path
    lib: Path | None = None  # pre-built families library; skips RepeatModeler when set
    engine: str = "rmblast"

    @computed_field  # type: ignore[prop-decorator]
    @property
    def name(self) -> str:
        return self.genome.stem

    settings_customise_sources = _pyproject_settings_sources("mask-repeats")


# ---------------------------------------------------------------------------
# Functional annotation settings (eukan func-annot)
# ---------------------------------------------------------------------------


class FunctionalConfig(_StepRunSettings):
    """Configuration for functional annotation."""

    model_config = SettingsConfigDict(
        env_prefix="EUKAN_FUNC_",
        extra="ignore",
    )

    proteins: Path
    homology_db: str = "uniprot"  # "uniprot" → phmmer vs SwissProt; "kofam" → hmmscan vs KOfam
    uniprot_db: Path = Path("databases/uniprot_sprot.faa")
    kofam_db: Path = Path("databases/kofam_eukaryote.hmm")
    ko_list_path: Path = Path("databases/ko_list.tsv")
    pfam_db: Path = Path("databases/Pfam-A.hmm")
    gff3_path: Path | None = None
    evalue: str = "1e-1"

    @model_validator(mode="after")
    def _validate_homology_db(self) -> FunctionalConfig:
        if self.homology_db not in ("uniprot", "kofam"):
            raise ValueError(
                f"homology_db must be 'uniprot' or 'kofam', got {self.homology_db!r}"
            )
        return self

    settings_customise_sources = _pyproject_settings_sources("func-annot")


# ---------------------------------------------------------------------------
# Submission-prep settings (eukan prep-submission)
# ---------------------------------------------------------------------------


class SubmissionConfig(_StepRunSettings):
    """Configuration for NCBI submission preparation (eukan prep-submission).

    Wraps NCBI's table2asn validator. Auto-discovers ``genome`` from
    ``eukan-run.json`` and ``gff3`` from ``work_dir`` (preferring
    ``final.mod.gff3`` over ``final.gff3``) when not given explicitly.
    """

    model_config = SettingsConfigDict(
        env_prefix="EUKAN_SUBMIT_",
        extra="ignore",
    )

    # --- Inputs (auto-discoverable) ---
    genome: Path | None = None
    gff3: Path | None = None
    template: Path  # required, no auto-discovery

    # --- Source qualifiers (one of organism+isolate or source_info required) ---
    organism: str | None = None
    isolate: str | None = None
    source_info: str | None = None  # raw -j override; supersedes organism/isolate

    # --- table2asn flags (defaults match the standard NCBI submission recipe) ---
    cleanup: str = "befw"
    mode: str = "n"
    assembly_type: str = "r10k"
    locus_tag_prefix: str | None = None

    # --- Output ---
    output_file: Path | None = None  # default: <output_dir>/<genome-stem>.sqn
    output_dir: Path = Field(default_factory=lambda: Path.cwd() / "submission")
    extra_args: list[str] = Field(default_factory=list)

    # --- GFF3 preprocessing ---
    cleanup_gff3: bool = True  # strip UniProt cruft, drop CDS-less mRNAs, etc.

    @model_validator(mode="after")
    def _discover_inputs(self) -> SubmissionConfig:
        log = get_logger(__name__)
        from eukan.infra import artifacts
        from eukan.infra.artifacts import Artifact
        from eukan.infra.layout import sibling_step_dir
        from eukan.infra.manifest import load_manifest

        if self.genome is None:
            # Manifest lives at the run-dir root (manifest_dir); fall back
            # to legacy per-step locations for older runs.
            for candidate_dir in (
                self.manifest_dir,
                sibling_step_dir(self.work_dir, "annotate"),
                self.work_dir,
            ):
                manifest = load_manifest(candidate_dir)
                if manifest and manifest.genome:
                    discovered = Path(manifest.genome)
                    object.__setattr__(self, "genome", discovered)
                    log.info("Auto-discovered genome from eukan-run.json: %s", discovered)
                    break

        if self.gff3 is None:
            preferred = artifacts.find(self.work_dir, Artifact.FINAL_FUNC_GFF3)
            fallback = artifacts.find(self.work_dir, Artifact.FINAL_GFF3)
            if preferred is not None:
                object.__setattr__(self, "gff3", preferred)
                log.info("Auto-discovered annotated GFF3: %s", preferred)
            elif fallback is not None:
                object.__setattr__(self, "gff3", fallback)
                log.warning(
                    "Using %s (no functional annotations). "
                    "Run `eukan func-annot` first for product names.",
                    fallback,
                )

        if self.genome is None:
            raise ValueError(
                "genome not set and not discoverable. "
                "Pass --genome or run from a directory containing eukan-run.json."
            )
        if self.gff3 is None:
            raise ValueError(
                "gff3 not set and not discoverable. "
                "Pass --gff3 or run `eukan annotate` (and optionally `eukan func-annot`) first."
            )

        if self.output_file is None:
            object.__setattr__(self, "output_file", self.output_dir / f"{self.genome.stem}.sqn")

        return self

    settings_customise_sources = _pyproject_settings_sources("prep-submission")
