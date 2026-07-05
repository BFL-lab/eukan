"""Assembly pipeline: read mapping → de novo + genome-guided assembly → SL
trans-splice cut → combinr consolidation."""

from __future__ import annotations

from pathlib import Path

from eukan.assembly.combinr import run_combinr
from eukan.assembly.defuse import run_defuse
from eukan.assembly.jaccard import run_jaccard
from eukan.assembly.max_intron import run_max_intron_split
from eukan.assembly.minimap2 import map_reads_minimap2, map_transcripts_minimap2
from eukan.assembly.sl_acceptors import detect_sl_acceptors
from eukan.assembly.sl_cut import run_sl_cut
from eukan.assembly.strand_correction import run_strand_correction
from eukan.assembly.trinity import run_trinity
from eukan.infra.artifacts import Artifact
from eukan.infra.manifest import ASSEMBLY, step_key
from eukan.infra.pipeline import (
    StepSpec,
    run_simple_pipeline,
)
from eukan.settings import AssemblyConfig

# ---------------------------------------------------------------------------
# Per-step input resolvers (for fingerprint-based resume invalidation)
# ---------------------------------------------------------------------------
# Each names the files a step consumes from upstream steps; the driver records
# their fingerprint on completion and re-runs the step on resume if they change
# (eukan/infra/steps.py::fingerprint_inputs). This is what makes a non-forced
# resume self-heal when, say, the assembly is rebuilt: map_transcripts re-maps
# the new transcripts, and the changed BAM cascades through strand_correct →
# sl_detect → sl_cut → combinr, each re-running as its input flips.


def _resolved_transcript_queries(config: AssemblyConfig) -> list[Path]:
    """The transcript FASTAs map_transcripts maps (the jaccard sibling if present)."""
    from eukan.assembly.bam_utils import _TRANSCRIPT_SETS, _resolve_query

    return [_resolve_query(config.work_dir, name) for name, _ in _TRANSCRIPT_SETS]


def _denovo_genome_bams(config: AssemblyConfig) -> list[Path]:
    """The de novo transcript→genome BAMs map_transcripts produces."""
    from eukan.assembly.sl_cut import _DENOVO_BAMS

    return [config.work_dir / b for b in _DENOVO_BAMS]


def _strand_correct_inputs(config: AssemblyConfig) -> list[Path]:
    # strand_correct converts each Trinity transcript->genome BAM to models, so a
    # changed BAM (re-mapping) re-runs it.
    return _denovo_genome_bams(config)


def _sl_detect_inputs(config: AssemblyConfig) -> list[Path]:
    # The read aligner BAM also feeds SL detection, but it is large and only
    # changes when the aligner is re-run (which cascades via the run flags), so
    # only the Trinity track BAMs — the map_transcripts result — are fingerprinted.
    return _denovo_genome_bams(config)


def _defuse_inputs(config: AssemblyConfig) -> list[Path]:
    # The models de-fusion reads per track: the homology-stranded models when
    # strand_correct produced them, else the raw genome GFF3. Whichever exist are
    # fingerprinted, so re-running strand_correct re-runs defuse.
    from eukan.assembly.tracks import mapped_transcript_stems

    wd = config.work_dir
    paths: list[Path] = []
    for stem in mapped_transcript_stems():
        paths += [wd / f"{stem}.stranded.gff3", wd / f"{stem}.gff3"]
    return paths


def _max_intron_split_inputs(config: AssemblyConfig) -> list[Path]:
    # Per track, max_intron_split resolves defuse > stranded > raw genome GFF3.
    # Fingerprint all candidates so any upstream change (strand_correct/defuse) re-runs it.
    from eukan.assembly.tracks import mapped_transcript_stems

    wd = config.work_dir
    paths: list[Path] = []
    for stem in mapped_transcript_stems():
        paths += [
            wd / f"{stem}.defuse.gff3",
            wd / f"{stem}.stranded.gff3",
            wd / f"{stem}.gff3",
        ]
    return paths


def _sl_cut_inputs(config: AssemblyConfig) -> list[Path]:
    # sl_cut reads each track's max-intron-split models plus the SL acceptor sites.
    from eukan.assembly.tracks import mapped_transcript_stems

    wd = config.work_dir
    paths: list[Path] = [
        wd / f"{stem}.maxintron.gff3" for stem in mapped_transcript_stems()
    ]
    paths.append(wd / Artifact.SL_ACCEPTORS.value)
    return paths


def _combinr_inputs(config: AssemblyConfig) -> list[Path]:
    from eukan.assembly.combinr import _CUT_MODELS

    return [config.work_dir / m for m in _CUT_MODELS]


# Scalar (non-file) inputs folded into a step's resume fingerprint, so changing the
# value re-runs the step even when its input files are byte-identical. ``-M`` is
# enforced post-mapping (StringTie reads a bounded BAM; max_intron_split splits
# models; combinr passes --max-intron), so tightening it must re-run exactly those
# steps — but NOT the read/transcript mappers (segemehl ignores --max-intron
# natively, so a scalar there would force a multi-hour re-map for an identical BAM).


def _max_intron_split_scalars(config: AssemblyConfig) -> list[str]:
    # max_intron_split hard-imposes -M on the transcript models, and min_sl_fragment
    # gates the tiny-fragment drop; changing either re-runs it (and cascades to
    # sl_cut/combinr). Its declared output makes this fingerprint live on resume.
    return [
        f"max_intron_len={config.max_intron_len}",
        f"min_sl_fragment={config.min_sl_fragment}",
    ]


def _combinr_scalars(config: AssemblyConfig) -> list[str]:
    # combinr passes --max-intron and --stringent-overlap; changing either re-runs the
    # consolidation (the cut models are byte-identical, so only a scalar invalidates it).
    return [
        f"max_intron_len={config.max_intron_len}",
        f"combinr_stringent_overlap={config.combinr_stringent_overlap}",
    ]


def _sl_cut_scalars(config: AssemblyConfig) -> list[str]:
    return [f"min_sl_fragment={config.min_sl_fragment}"]


def _defuse_scalars(config: AssemblyConfig) -> list[str]:
    # A changed de-fusion knob re-runs defuse (and cascades to sl_cut/combinr).
    return [
        f"defuse={config.defuse}",
        f"defuse_overlap_tolerance={config.defuse_overlap_tolerance}",
        f"defuse_blastx_evalue={config.defuse_blastx_evalue}",
        f"min_sl_fragment={config.min_sl_fragment}",
    ]


def _read_mapping_step() -> StepSpec:
    """The read-mapping StepSpec: minimap2 (splice:sr), with the non-canonical
    escalation handled inside :func:`minimap2.map_reads_minimap2`."""
    return StepSpec(
        "minimap2", map_reads_minimap2,
        "minimap2_Aligned.sortedByCoord.out.bam", "-A / --run-minimap2",
    )


def _steps_for() -> list[StepSpec]:
    """Assembly steps: minimap2 → trinity (de novo + genome-guided) → jaccard
    → map_transcripts → strand_correct → defuse → max_intron_split → sl_detect
    → sl_cut → combinr."""
    return [
        _read_mapping_step(),
        # Trinity assembles both de novo and genome-guided transcript FASTAs;
        # genome-guided reads the aligner BAM (config.aligner_bam) only to cluster
        # reads per locus, so its output is still a transcript FASTA mapped back to
        # the genome by map_transcripts. Declared output is trinity-denovo.fasta
        # (produced last), so its presence marks both modes done on resume.
        StepSpec("trinity", run_trinity, "trinity-denovo.fasta", "-T / --run-trinity"),
        # No declared output: jaccard rewrites each Trinity FASTA into a
        # ``.jaccard.fasta`` sibling, but on single-end input (or zero clips) it
        # legitimately writes nothing, so stale-output validation must not fire.
        # Having no declared output it always re-runs on resume, so a changed
        # ``jaccard_greediness`` re-clips and cascades to map_transcripts (which
        # fingerprints the rewritten ``.jaccard.fasta``) — a scalar here would be
        # inert, since is_step_complete never reaches the fingerprint check.
        StepSpec("jaccard", run_jaccard, None, "--run-jaccard"),
        StepSpec(
            "map_transcripts", map_transcripts_minimap2,
            "trinity-denovo.genome.bam", "--run-map-transcripts",
            inputs=_resolved_transcript_queries,
        ),
        # No declared output: strand_correct always converts each Trinity track BAM
        # to <track>.genome.gff3 (for the SL cut), but only writes the *.stranded.gff3
        # homology-corrected models when --uniprot is given on an unstranded library,
        # so stale-output validation must not fire when it's a no-op.
        StepSpec(
            "strand_correct", run_strand_correction, None, "--run-strand-correct",
            inputs=_strand_correct_inputs,
        ),
        # Homology de-fusion: splits chimeric transcripts flagged by >=2 distinct
        # non-overlapping SwissProt hits. Always writes <track>.genome.defuse.gff3
        # when it runs (copy-through if nothing split), so it can declare an output +
        # scalars for precise resume; skipped entirely (run_assembly skip) unless
        # --defuse + --uniprot are set.
        StepSpec(
            "defuse", run_defuse, "trinity-denovo.genome.defuse.gff3", "--run-defuse",
            inputs=_defuse_inputs, scalars=_defuse_scalars,
        ),
        # Max-intron split: strand-agnostic, SL-independent sanitization that breaks
        # any model carrying an intron > -M into separate genes (the de novo segemehl
        # path has no native intron bound). Always writes <track>.maxintron.gff3
        # (copy-through when nothing is cut), so — like defuse — it declares an output
        # + scalars for precise resume: a changed -M re-runs it (and cascades on).
        StepSpec(
            "max_intron_split", run_max_intron_split,
            "trinity-denovo.genome.maxintron.gff3", "--run-max-intron-split",
            inputs=_max_intron_split_inputs, scalars=_max_intron_split_scalars,
        ),
        # No declared output for sl_detect/sl_cut: with no SL signal sl_detect writes a
        # header-only sl_acceptors.gff3 (zero features), and sl_cut may legitimately
        # write nothing for a track with no models (per-stem skip), so it always
        # re-runs on resume — its scalar fingerprint is thus inert but kept for intent.
        StepSpec(
            "sl_detect", detect_sl_acceptors, None, "--run-sl-detect",
            inputs=_sl_detect_inputs,
        ),
        StepSpec(
            "sl_cut", run_sl_cut, None, "--run-sl-cut",
            inputs=_sl_cut_inputs, scalars=_sl_cut_scalars,
        ),
        StepSpec(
            "combinr", run_combinr,
            Artifact.NR_TRANSCRIPTS_FASTA.value, "--run-combinr",
            inputs=_combinr_inputs, scalars=_combinr_scalars,
        ),
    ]


# Within-pipeline data dependencies: re-running a step (via its --run-* flag)
# invalidates every step that consumes its output, so those are forced too. The
# read aligner is cascaded because Trinity genome-guided (config.aligner_bam) and
# SL read-side detection read its BAM. Each edge below is "producer -> direct
# consumers".
_DOWNSTREAM: dict[str, tuple[str, ...]] = {
    "minimap2": ("trinity", "sl_detect"),
    "trinity": ("jaccard",),
    # jaccard rewrites the Trinity FASTAs into .jaccard.fasta siblings that
    # map_transcripts maps, so re-running it re-runs map_transcripts (which
    # cascades on to strand_correct etc.).
    "jaccard": ("map_transcripts",),
    "map_transcripts": ("strand_correct", "sl_detect"),
    # strand_correct feeds defuse, which feeds max_intron_split, which feeds sl_cut.
    # When defuse is disabled it is skipped but stays in the cascade, so the transitive
    # closure still reaches max_intron_split (which falls back to the stranded models).
    "strand_correct": ("defuse",),
    "defuse": ("max_intron_split",),
    "max_intron_split": ("sl_cut",),
    "sl_detect": ("sl_cut",),
    "sl_cut": ("combinr",),
}


def _expand_downstream(selected: set[str]) -> set[str]:
    """*selected* plus the transitive closure of their ``_DOWNSTREAM`` dependents."""
    out = set(selected)
    stack = list(selected)
    while stack:
        for dep in _DOWNSTREAM.get(stack.pop(), ()):
            if dep not in out:
                out.add(dep)
                stack.append(dep)
    return out


def force_steps_from_run_flags(
    *,
    run_minimap2: bool = False,
    run_trinity: bool = False,
    run_jaccard: bool = False,
    run_map_transcripts: bool = False,
    run_strand_correct: bool = False,
    run_defuse: bool = False,
    run_max_intron_split: bool = False,
    run_sl_detect: bool = False,
    run_sl_cut: bool = False,
    run_combinr: bool = False,
    force: bool = False,
) -> list[str]:
    """Translate ``--run-X`` / ``--force`` flags into manifest keys to force.

    ``--force`` re-runs every step. An individual ``--run-X`` re-runs that step
    *and every downstream step that consumes its output* (``_DOWNSTREAM``): e.g.
    ``--run-map-transcripts`` re-maps the Trinity transcripts and then re-runs
    strand_correct, defuse, max_intron_split, sl_detect, sl_cut, and combinr, which
    read the new BAM (directly or transitively). Returns keys in pipeline order.
    """
    steps = _steps_for()
    valid = {s.name for s in steps}
    flags = {
        "minimap2": run_minimap2, "trinity": run_trinity,
        "jaccard": run_jaccard,
        "map_transcripts": run_map_transcripts,
        "strand_correct": run_strand_correct, "defuse": run_defuse,
        "max_intron_split": run_max_intron_split,
        "sl_detect": run_sl_detect, "sl_cut": run_sl_cut, "combinr": run_combinr,
    }
    selected = {name for name, on in flags.items() if on and name in valid}
    # Individual --run-X flags take precedence over --force (scope to the
    # selected steps + their downstream); --force alone re-runs everything.
    if selected:
        expanded = _expand_downstream(selected)
        return [step_key(ASSEMBLY, s.name) for s in steps if s.name in expanded]
    if force:
        return [step_key(ASSEMBLY, s.name) for s in steps]
    return []


def run_assembly(
    config: AssemblyConfig,
    *,
    force_steps: list[str] | None = None,
) -> None:
    """Run the assembly pipeline with manifest tracking."""
    run_simple_pipeline(
        ASSEMBLY, _steps_for(), config,
        force_steps=force_steps,
        skip=lambda s: (s.name == "jaccard" and not config.jaccard_clip)
        or (s.name == "defuse" and not (config.defuse and config.uniprot_db is not None)),
    )
