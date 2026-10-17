"""Annotation pipeline orchestration: step ordering, concurrency, and manifest tracking."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

from eukan.annotation.alignment import align_proteins
from eukan.annotation.augustus import run_augustus
from eukan.annotation.consensus import build_consensus_models
from eukan.annotation.genemark import run_genemark
from eukan.annotation.orf import create_transcriptome_orf_db
from eukan.annotation.snap import run_codingquarry, run_snap
from eukan.gff.io import count_gff3_features, featuredb2gff3_file
from eukan.infra.artifacts import Artifact, find, masked_genome
from eukan.infra.concurrency import parallel_map
from eukan.infra.logging import get_logger
from eukan.infra.manifest import (
    ANNOTATION,
    RunManifest,
    get_or_create_manifest,
    save_manifest,
    step_key,
)
from eukan.infra.pipeline import run_orchestrated_step
from eukan.infra.steps import step_dir
from eukan.settings import PipelineConfig
from eukan.validation import sanitize_genome_fasta, validate_fasta

log = get_logger(__name__)


# ---------------------------------------------------------------------------
# ORF finding step
# ---------------------------------------------------------------------------


def find_orfs(config: PipelineConfig, trans_gff3: Path) -> Path:
    """Find ORFs in transcript assemblies."""
    from eukan.validation import validate_gff

    output = "transcript_orfs.gff3"
    sdir = step_dir(config.work_dir, "orf_finder")
    log.info("Finding ORFs in transcript assemblies...")
    validate_gff(trans_gff3)

    orfs = create_transcriptome_orf_db(str(trans_gff3), str(config.genome), genetic_code=int(config.genetic_code))
    featuredb2gff3_file(orfs, sdir / output)
    return sdir / output


# ---------------------------------------------------------------------------
# Step helpers
# ---------------------------------------------------------------------------


def _run_step(
    config: PipelineConfig,
    manifest: RunManifest,
    name: str,
    fn: Callable,
    *args,
    forced: frozenset[str] = frozenset(),
    **kwargs,
) -> Path:
    """Run an annotation step via the shared pipeline helper.

    The step re-runs from scratch when its manifest key is in *forced*. All
    annotation steps produce a GFF3 output, so the result is always a Path
    (never None).
    """
    key = step_key(ANNOTATION, name)
    result = run_orchestrated_step(
        config.manifest_dir, manifest, key,
        fn, config, *args,
        step_dir=config.work_dir / name,
        force=key in forced,
        **kwargs,
    )
    assert result is not None, f"annotation step {name!r} returned no output"
    return result


def _log_prediction_count(label: str, gff3_path: Path) -> None:
    """Log the number of gene predictions in a GFF3 file."""
    log.info("%s: %d gene predictions", label, count_gff3_features(gff3_path))


def _run_concurrent_steps(
    config: PipelineConfig,
    manifest: RunManifest,
    tasks: list[tuple[str, Callable, tuple, dict]],
    forced: frozenset[str],
) -> dict[str, Path]:
    """Run multiple independent steps concurrently."""

    def run(task: tuple[str, Callable, tuple, dict]) -> Path:
        name, fn, args, kwargs = task
        return _run_step(config, manifest, name, fn, *args, forced=forced, **kwargs)

    paths = parallel_map(run, tasks)
    return {name: path for (name, *_), path in zip(tasks, paths, strict=True)}


# ---------------------------------------------------------------------------
# Pipeline orchestration
# ---------------------------------------------------------------------------


def run_annotation_pipeline(
    config: PipelineConfig,
    force_steps: list[str] | None = None,
) -> Path:
    """Run the full annotation pipeline.

    Features:
    - Writes eukan-run.json manifest for tracking/reproducibility
    - Skips completed steps on resume
    - Cleans up interrupted steps automatically
    - Runs independent steps concurrently where possible

    Args:
        config: Pipeline configuration.
        force_steps: Optional manifest keys of steps to re-run from scratch,
            even if previously completed (see ``run_orchestrated_step``'s
            ``force``).

    Returns:
        Path to the final GFF3 output.
    """
    validate_fasta(config.genome)

    # If the user ran `eukan mask-repeats` previously and is now passing the
    # unmasked genome, point them at the masked sibling. Don't auto-swap —
    # explicit input is safer than a silent path substitution.
    if not config.genome.name.endswith(".masked.fasta"):
        masked_sibling = masked_genome(config.work_dir, config.genome.stem)
        if masked_sibling.exists() and masked_sibling != config.genome:
            log.info(
                "Note: found %s alongside the input genome — pass it via -g for "
                "repeat-aware prediction.",
                masked_sibling.name,
            )

    # Sanitize genome headers (strip descriptions that break GFF tools)
    sanitized_genome = sanitize_genome_fasta(config.genome, config.work_dir)
    if sanitized_genome != config.genome:
        config = config.model_copy(update={"genome": sanitized_genome})

    if config.has_transcripts:
        # has_transcripts is True iff all three are non-None
        assert (
            config.transcripts_fasta is not None
            and config.transcripts_gff is not None
            and config.rnaseq_hints is not None
        )
        log.info("Transcript evidence: %s, %s, %s",
                 config.transcripts_fasta.name, config.transcripts_gff.name, config.rnaseq_hints.name)
    else:
        log.info("No transcript evidence found. Running protein-only annotation.")

    # Load or create shared manifest (used by all pipelines in this work_dir).
    # Each step decides on its own record whether it still has to run, so a
    # resume with nothing left to do skips every step.
    manifest = get_or_create_manifest(config.manifest_dir, config)
    save_manifest(config.manifest_dir, manifest)

    try:
        result = _execute_steps(config, manifest, frozenset(force_steps or ()))
        manifest.status = "completed"
        manifest.finished_at = manifest.steps[
            max(manifest.steps, key=lambda k: manifest.steps[k].finished_at or "")
        ].finished_at
        save_manifest(config.manifest_dir, manifest)
        return result
    except Exception:
        manifest.status = "failed"
        save_manifest(config.manifest_dir, manifest)
        raise


# Single source of truth: step name (no prefix) -> --run-* CLI flag that
# (a) forces re-run when the user passes the flag, and
# (b) is shown in error messages when a step's output is broken.
# Multiple step names mapping to the same flag means the flag forces all
# of them as a group (e.g. --run-genemark also re-runs orf_finder).
_ANNOTATION_STEP_FLAGS: dict[str, str] = {
    "genemark":             "--run-genemark",
    "orf_finder":           "--run-genemark",
    "prot_align":           "--run-prot-align",
    "prot_align_ssp":       "--run-prot-align",
    "augustus":             "--run-augustus",
    "snap":                 "--run-snap",
    "codingquarry":         "--run-snap",
    "consensus_models": "--run-consensus",
}


def force_steps_from_run_flags(
    *,
    spaln_ssp: bool = False,
    run_genemark: bool = False,
    run_prot_align: bool = False,
    run_augustus: bool = False,
    run_snap: bool = False,
    run_consensus: bool = False,
) -> list[str]:
    """Translate per-flag booleans into manifest step keys to force.

    Steps grouped under the same flag are all forced together, except
    that prot_align/prot_align_ssp is selected based on ``spaln_ssp``.
    """
    flag_states = {
        "--run-genemark":   run_genemark,
        "--run-prot-align": run_prot_align,
        "--run-augustus":   run_augustus,
        "--run-snap":       run_snap,
        "--run-consensus":  run_consensus,
    }
    forced: list[str] = []
    for name, flag in _ANNOTATION_STEP_FLAGS.items():
        if not flag_states.get(flag, False):
            continue
        if name == "prot_align" and spaln_ssp:
            continue
        if name == "prot_align_ssp" and not spaln_ssp:
            continue
        forced.append(step_key(ANNOTATION, name))
    return forced


def _execute_steps(
    config: PipelineConfig,
    manifest: RunManifest,
    forced: frozenset[str] = frozenset(),
) -> Path:
    """Execute annotation steps in dependency phases.

    Each phase reads inputs from ``ev`` and returns the keys it produced. A
    step whose manifest key is in *forced* re-runs from scratch. Every step
    declares the files it reads as its ``input_files``, so a step re-runs when
    one changed since it completed — an upstream step's rewritten output, an
    edited protein set, new RNA-seq hints — and that re-run cascades on through
    the steps that read its output.
    Behavior matrix preserved verbatim from the previous nested form:

      - has_transcripts + is_fungus
            ORF || GeneMark, then spaln (intron-hinted), AUGUSTUS,
            SNAP || CodingQuarry, consensus(spaln, augustus, snap, cq, trans)
      - has_transcripts + not is_fungus
            ORF || GeneMark, spaln, AUGUSTUS, *no SNAP/CodingQuarry*,
            consensus(spaln, augustus, trans)
      - no transcripts + (is_fungus | is_protist)
            GeneMark, spaln, AUGUSTUS, SNAP || CodingQuarry,
            consensus(spaln, augustus, snap, cq)
      - no transcripts + neither
            GeneMark, spaln, AUGUSTUS, SNAP, consensus(spaln, augustus, snap)
    """
    ev: dict[str, Path] = {}
    ev = _phase_orf_and_genemark(config, manifest, ev, forced)
    ev = _phase_protein_alignment(config, manifest, ev, forced)
    ev = _phase_augustus(config, manifest, ev, forced)
    ev = _phase_snap_codingquarry(config, manifest, ev, forced)
    return _phase_consensus(config, manifest, ev, forced)


def _phase_orf_and_genemark(
    config: PipelineConfig, manifest: RunManifest, ev: dict[str, Path],
    forced: frozenset[str],
) -> dict[str, Path]:
    """Phase 1: ORF finding (transcripts only) + GeneMark, run concurrently."""
    genome = config.genome
    if config.has_transcripts:
        assert config.transcripts_gff is not None and config.rnaseq_hints is not None
        trans, hints = config.transcripts_gff, config.rnaseq_hints
        results = _run_concurrent_steps(config, manifest, [
            ("orf_finder", find_orfs, (trans,), {"input_files": [genome, trans]}),
            ("genemark", run_genemark, (hints,), {"input_files": [genome, hints]}),
        ], forced)
        ev = {**ev, "transcriptORFs": results["orf_finder"], "genemark": results["genemark"]}
    else:
        ev = {
            **ev,
            "genemark": _run_step(
                config, manifest, "genemark", run_genemark,
                input_files=[genome], forced=forced,
            ),
        }
    _log_prediction_count("GeneMark", ev["genemark"])
    return ev


def _phase_protein_alignment(
    config: PipelineConfig, manifest: RunManifest, ev: dict[str, Path],
    forced: frozenset[str],
) -> dict[str, Path]:
    """Phase 2: spliced protein alignment via spaln (intron-hinted) or gth."""
    prot_step = "prot_align_ssp" if config.spaln_ssp else "prot_align"
    # The GeneMark models pick spaln or gth and, without intron hints, stand in for them.
    inputs = [config.genome, *config.proteins, ev["genemark"]]
    spaln_extra: tuple = ()
    if config.has_transcripts:
        intron_hints = config.work_dir / "genemark" / "introns.gff"
        spaln_extra = (intron_hints if intron_hints.exists() else None,)
        inputs.append(intron_hints)
    if config.spaln_ssp and config.transcripts_fasta:
        inputs.append(config.transcripts_fasta)  # --spsp trains spaln on the transcripts
    spaln_path = _run_step(
        config, manifest, prot_step, align_proteins,
        ev["genemark"], config.proteins, *spaln_extra,
        input_files=inputs, forced=forced,
    )
    _log_prediction_count("spaln", spaln_path)
    return {**ev, "spaln": spaln_path}


def _phase_augustus(
    config: PipelineConfig, manifest: RunManifest, ev: dict[str, Path],
    forced: frozenset[str],
) -> dict[str, Path]:
    """Phase 3: AUGUSTUS training and prediction."""
    evidence: tuple[Path, ...] = (ev["genemark"], ev["spaln"])
    if config.has_transcripts:
        evidence += (ev["transcriptORFs"],)
    # Besides its evidence, AUGUSTUS reads the hint files: the protein hints beside
    # the alignment, the RNA-seq hints, and the RepeatMasker hints and splice-site
    # summary when another pipeline left them.
    inputs = [config.genome, *evidence, ev["spaln"].parent / "hints_protein.gff"]
    if config.rnaseq_hints:
        inputs.append(config.rnaseq_hints)
    for artifact in (Artifact.REPEATMASK_HINTS, Artifact.SPLICE_SUMMARY):
        found = find(config.work_dir, artifact)
        if found is not None:
            inputs.append(found)
    aug_path = _run_step(
        config, manifest, "augustus", run_augustus, *evidence,
        input_files=inputs, forced=forced,
    )
    _log_prediction_count("AUGUSTUS", aug_path)
    return {**ev, "augustus": aug_path}


def _phase_snap_codingquarry(
    config: PipelineConfig, manifest: RunManifest, ev: dict[str, Path],
    forced: frozenset[str],
) -> dict[str, Path]:
    """Phase 4: SNAP, plus CodingQuarry for fungus/protist; skipped for has_t & non-fungus."""
    has_t = config.has_transcripts

    genome = config.genome
    if has_t and config.is_fungus:
        assert config.transcripts_gff is not None
        snap_ev: tuple[Path, ...] = (ev["augustus"], ev["spaln"], ev["transcriptORFs"])
        cq_ev = config.transcripts_gff
        results = _run_concurrent_steps(config, manifest, [
            ("snap", run_snap, snap_ev, {"input_files": [genome, *snap_ev]}),
            ("codingquarry", run_codingquarry, (cq_ev,), {"input_files": [genome, cq_ev]}),
        ], forced)
        ev = {**ev, "snap": results["snap"], "codingquarry": results["codingquarry"]}
        _log_prediction_count("SNAP", ev["snap"])
        _log_prediction_count("CodingQuarry", ev["codingquarry"])
    elif not has_t and (config.is_fungus or config.is_protist):
        snap_ev = (ev["augustus"], ev["spaln"])
        cq_ev = ev["augustus"]
        results = _run_concurrent_steps(config, manifest, [
            ("snap", run_snap, snap_ev, {"input_files": [genome, *snap_ev]}),
            ("codingquarry", run_codingquarry, (cq_ev,), {"input_files": [genome, cq_ev]}),
        ], forced)
        ev = {**ev, "snap": results["snap"], "codingquarry": results["codingquarry"]}
        _log_prediction_count("SNAP", ev["snap"])
        _log_prediction_count("CodingQuarry", ev["codingquarry"])
    elif not has_t:
        snap_ev = (ev["augustus"], ev["spaln"])
        snap_path = _run_step(
            config, manifest, "snap", run_snap, *snap_ev,
            input_files=[genome, *snap_ev], forced=forced,
        )
        ev = {**ev, "snap": snap_path}
        _log_prediction_count("SNAP", snap_path)
    # else: has_t but not fungus -- SNAP/CodingQuarry are skipped entirely.

    return ev


def _consensus_scalars(config: PipelineConfig) -> list[str]:
    """Config values that change combinr consensus output, folded into the consensus
    step's resume fingerprint so a change re-runs it on resume (in addition to the
    explicit ``--run-consensus``). ``combinr_stringent_overlap`` tunes the
    ``--alt-splice`` isoform grouping; weights and the genetic code feed the consensus
    engine directly."""
    return [
        f"weights={config.weights}",
        f"genetic_code={config.genetic_code}",
        f"combinr_stringent_overlap={config.combinr_stringent_overlap}",
    ]


def _phase_consensus(
    config: PipelineConfig, manifest: RunManifest, ev: dict[str, Path],
    forced: frozenset[str],
) -> Path:
    """Phase 5: combinr consensus. Argument order varies with which evidence ran."""
    consensus_args: list[Path] = [ev["spaln"], ev["augustus"]]
    if "snap" in ev:
        consensus_args.append(ev["snap"])
    if "codingquarry" in ev:
        consensus_args.append(ev["codingquarry"])

    transcripts: Path | None = None
    if config.has_transcripts:
        assert config.transcripts_gff is not None
        transcripts = config.transcripts_gff

    # The consensus also patches in the transcript ORFs, whenever orf_finder left them.
    inputs = [
        config.genome, *consensus_args, config.work_dir / "orf_finder" / "transcript_orfs.gff3",
    ]
    if transcripts is not None:
        inputs.append(transcripts)
    return _run_step(
        config, manifest, "consensus_models", build_consensus_models,
        *consensus_args, transcripts=transcripts,
        input_files=inputs, input_scalars=_consensus_scalars(config), forced=forced,
    )
