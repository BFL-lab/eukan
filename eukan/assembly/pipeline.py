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
from eukan.assembly.tracks import (
    CUT_GFF3,
    DEFUSE_GFF3,
    MAXINTRON_GFF3,
    RAW_GFF3,
    STRANDED_GFF3,
    TRANSCRIPT_SETS,
    jaccard_output_name,
    mapped_transcript_stems,
    resolve_query,
)
from eukan.assembly.trinity import run_trinity
from eukan.infra.artifacts import Artifact
from eukan.infra.logging import get_logger
from eukan.infra.manifest import ASSEMBLY
from eukan.infra.pipeline import (
    StepSpec,
    run_simple_pipeline,
)
from eukan.infra.pipeline import (
    force_steps_from_run_flags as _force_steps_from_run_flags,
)
from eukan.settings import AssemblyConfig

log = get_logger(__name__)

# ---------------------------------------------------------------------------
# Per-step input resolvers (for fingerprint-based resume invalidation)
# ---------------------------------------------------------------------------
# Each names the files a step consumes from upstream steps; the driver records
# their fingerprint on completion and re-runs the step on resume if they change
# (eukan/infra/steps.py::fingerprint_inputs). This is what makes a non-forced
# resume self-heal when, say, the assembly is rebuilt: map_transcripts re-maps
# the new transcripts, and the changed BAM cascades through strand_correct →
# sl_detect → sl_cut → combinr, each re-running as its input flips. A step with no
# declared output resumes on this fingerprint alone, so its resolver also names the
# files the step writes: one removed or rewritten since then re-runs it.


def _resolved_transcript_queries(config: AssemblyConfig) -> list[Path]:
    """The transcript FASTAs map_transcripts maps (the jaccard sibling if present)."""
    return [resolve_query(config.work_dir, name) for name, _ in TRANSCRIPT_SETS]


def _denovo_genome_bams(config: AssemblyConfig) -> list[Path]:
    """The de novo transcript→genome BAMs map_transcripts produces."""
    return [config.work_dir / bam for _, bam in TRANSCRIPT_SETS]


def _jaccard_inputs(config: AssemblyConfig) -> list[Path]:
    # jaccard clips each Trinity FASTA into the .jaccard.fasta sibling it writes
    # (run_assembly removes those while jaccard does not run).
    wd = config.work_dir
    return [
        wd / name
        for query, _ in TRANSCRIPT_SETS
        for name in (query, jaccard_output_name(query))
    ]


def _strand_correct_inputs(config: AssemblyConfig) -> list[Path]:
    # strand_correct converts each Trinity transcript->genome BAM to models (and,
    # with --uniprot, corrects them against the protein DB), so a changed BAM
    # (re-mapping) or DB re-runs it; it writes <stem>.gff3 and <stem>.stranded.gff3.
    wd = config.work_dir
    paths = _denovo_genome_bams(config)
    for stem in mapped_transcript_stems():
        paths += [wd / f"{stem}{RAW_GFF3}", wd / f"{stem}{STRANDED_GFF3}"]
    if config.uniprot_db is not None:
        paths.append(config.uniprot_db)
    return paths


def _sl_detect_inputs(config: AssemblyConfig) -> list[Path]:
    # SL detection scans the read BAM and the Trinity track BAMs, takes the read-side
    # consensus from the soft-clip summary, and writes sl_acceptors.gff3. The read BAM
    # is large, but hashing it costs far less than the walk over it a re-run repeats.
    wd = config.work_dir
    return [
        wd / config.aligner_bam,
        *_denovo_genome_bams(config),
        wd / Artifact.SOFTCLIP_DIAGNOSTIC.value,
        wd / Artifact.SL_ACCEPTORS.value,
    ]


def _defuse_inputs(config: AssemblyConfig) -> list[Path]:
    # The models de-fusion reads per track: the homology-stranded models when
    # strand_correct produced them, else the raw genome GFF3. Whichever exist are
    # fingerprinted, so re-running strand_correct re-runs defuse.
    wd = config.work_dir
    paths: list[Path] = []
    for stem in mapped_transcript_stems():
        paths += [wd / f"{stem}{STRANDED_GFF3}", wd / f"{stem}{RAW_GFF3}"]
    return paths


def _max_intron_split_inputs(config: AssemblyConfig) -> list[Path]:
    # Per track, max_intron_split resolves defuse > stranded > raw genome GFF3.
    # Fingerprint all candidates so any upstream change (strand_correct/defuse) re-runs it.
    wd = config.work_dir
    paths: list[Path] = []
    for stem in mapped_transcript_stems():
        paths += [
            wd / f"{stem}{DEFUSE_GFF3}",
            wd / f"{stem}{STRANDED_GFF3}",
            wd / f"{stem}{RAW_GFF3}",
        ]
    return paths


def _sl_cut_inputs(config: AssemblyConfig) -> list[Path]:
    # sl_cut reads each track's max-intron-split models plus the SL acceptor sites
    # and writes each track's <stem>.cut.gff3.
    wd = config.work_dir
    paths: list[Path] = []
    for stem in mapped_transcript_stems():
        paths += [wd / f"{stem}{MAXINTRON_GFF3}", wd / f"{stem}{CUT_GFF3}"]
    paths.append(wd / Artifact.SL_ACCEPTORS.value)
    return paths


def _combinr_inputs(config: AssemblyConfig) -> list[Path]:
    # combinr consolidates each track's SL-cut models.
    return [config.work_dir / f"{stem}{CUT_GFF3}" for stem in mapped_transcript_stems()]


# Scalar (non-file) inputs folded into a step's resume fingerprint, so changing the
# value re-runs the step even when its input files are byte-identical. Every step
# that consumes ``-M`` carries it: minimap2 maps reads and transcripts with ``-G``,
# Trinity genome-guided takes it as ``--genome_guided_max_intron``, max_intron_split
# splits models at it, and combinr passes ``--max-intron``. A changed ``-M`` thus
# re-runs the read mapping, the assembly and the transcript mapping (hours on a
# large library) along with the model steps, so a looser limit actually recovers
# longer introns instead of reusing BAMs and assemblies capped at the old one.


def _max_intron_scalars(config: AssemblyConfig) -> list[str]:
    # minimap2 (reads and transcripts) and Trinity genome-guided bound introns at -M.
    # The trinity step re-runs both modes, though only genome-guided takes -M.
    return [f"max_intron_len={config.max_intron_len}"]


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


def _jaccard_scalars(config: AssemblyConfig) -> list[str]:
    # The clip-detection knobs and the floor below which a clipped piece is dropped.
    return [
        f"jaccard_greediness={config.jaccard_greediness}",
        f"jaccard_max_trough={config.jaccard_max_trough}",
        f"jaccard_min_delta={config.jaccard_min_delta}",
        f"jaccard_max_adaptive_trough={config.jaccard_max_adaptive_trough}",
        f"min_sl_fragment={config.min_sl_fragment}",
    ]


def _strand_correct_scalars(config: AssemblyConfig) -> list[str]:
    # A stranded library skips the correction; the rest tune the blastx and the vote.
    return [
        f"strand_specific={config.strand_specific}",
        f"genetic_code={config.genetic_code}",
        f"strand_blastx_evalue={config.strand_blastx_evalue}",
        f"min_strand_consensus={config.min_strand_consensus}",
    ]


def _sl_detect_scalars(config: AssemblyConfig) -> list[str]:
    return [
        f"sl_sequence={config.sl_sequence}",
        f"sl_adapter_filter={config.sl_adapter_filter}",
        f"adapter_sequences={config.adapter_sequences}",
        f"min_sl_clip_len={config.min_sl_clip_len}",
        f"min_sl_insertion_len={config.min_sl_insertion_len}",
        f"sl_cluster_window={config.sl_cluster_window}",
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


# The assembly steps, in execution order: minimap2 → trinity (de novo +
# genome-guided) → jaccard → map_transcripts → strand_correct → defuse →
# max_intron_split → sl_detect → sl_cut → combinr.
STEPS: tuple[StepSpec, ...] = (
    # Read mapping: minimap2 (splice:sr), with the non-canonical escalation handled
    # inside map_reads_minimap2. The soft-clip walk is skipped while its summary
    # exists, which lets the non-canonical re-map inside one run reuse the verdict. A
    # re-run of the completed step clears the summary so the walk rebuilds the
    # verdict, and the poly-A "reads" section it writes, from the new BAM.
    # polyA_diagnostic.json itself stays: map_transcripts writes its transcript
    # sections there.
    StepSpec(
        "minimap2", map_reads_minimap2,
        "minimap2_Aligned.sortedByCoord.out.bam", "-A / --run-minimap2",
        scalars=_max_intron_scalars,
        clears=(Artifact.SOFTCLIP_DIAGNOSTIC.value,),
    ),
    # Trinity assembles both de novo and genome-guided transcript FASTAs;
    # genome-guided reads the aligner BAM (config.aligner_bam) only to cluster
    # reads per locus, so its output is still a transcript FASTA mapped back to
    # the genome by map_transcripts. Declared output is trinity-denovo.fasta
    # (produced last), so its presence marks both modes done on resume.
    # run_trinity skips a mode whose FASTA exists (an interrupted run keeps the
    # finished mode), so a re-run of the completed step clears trinity-gg.fasta too.
    StepSpec(
        "trinity", run_trinity, "trinity-denovo.fasta", "-T / --run-trinity",
        clears=("trinity-gg.fasta",), scalars=_max_intron_scalars,
    ),
    # No declared output: jaccard rewrites each Trinity FASTA into a
    # ``.jaccard.fasta`` sibling, but on single-end input it legitimately writes
    # nothing, so stale-output validation must not fire. It resumes on its input
    # fingerprint instead: a changed knob re-clips and cascades to map_transcripts
    # (which fingerprints the rewritten ``.jaccard.fasta``).
    StepSpec(
        "jaccard", run_jaccard, None, "--run-jaccard",
        inputs=_jaccard_inputs, scalars=_jaccard_scalars,
    ),
    # The mapper reuses any per-set BAM that passes quickcheck, so a re-run of the
    # completed step clears the genome-guided BAM too, not only the declared one.
    StepSpec(
        "map_transcripts", map_transcripts_minimap2,
        "trinity-denovo.genome.bam", "--run-map-transcripts",
        inputs=_resolved_transcript_queries, scalars=_max_intron_scalars,
        clears=("trinity-gg.genome.bam",),
    ),
    # No declared output: strand_correct always converts each Trinity track BAM
    # to <track>.genome.gff3 (for the SL cut), but only writes the *.stranded.gff3
    # homology-corrected models when --uniprot is given on an unstranded library,
    # so stale-output validation must not fire when it's a no-op. It resumes on
    # its input fingerprint instead.
    StepSpec(
        "strand_correct", run_strand_correction, None, "--run-strand-correct",
        inputs=_strand_correct_inputs, scalars=_strand_correct_scalars,
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
    # any model carrying an intron > -M into separate genes (minimap2 also applies
    # -M as -G at mapping time). Always writes <track>.maxintron.gff3 (copy-through
    # when nothing is cut), so — like defuse — it declares an output + scalars for
    # precise resume: a changed -M re-runs it (and cascades on).
    StepSpec(
        "max_intron_split", run_max_intron_split,
        "trinity-denovo.genome.maxintron.gff3", "--run-max-intron-split",
        inputs=_max_intron_split_inputs, scalars=_max_intron_split_scalars,
    ),
    # No declared output for sl_detect/sl_cut: with no SL signal sl_detect writes a
    # header-only sl_acceptors.gff3 (zero features), and sl_cut may legitimately
    # write nothing for a track with no models (per-stem skip). They resume on
    # their input fingerprints instead.
    StepSpec(
        "sl_detect", detect_sl_acceptors, None, "--run-sl-detect",
        inputs=_sl_detect_inputs, scalars=_sl_detect_scalars,
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
)


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
    read the new BAM (directly or transitively). Individual ``--run-X`` flags
    take precedence over ``--force``. Returns keys in pipeline order.
    """
    return _force_steps_from_run_flags(
        ASSEMBLY, STEPS, force=force, downstream=_DOWNSTREAM,
        run_minimap2=run_minimap2, run_trinity=run_trinity, run_jaccard=run_jaccard,
        run_map_transcripts=run_map_transcripts, run_strand_correct=run_strand_correct,
        run_defuse=run_defuse, run_max_intron_split=run_max_intron_split,
        run_sl_detect=run_sl_detect, run_sl_cut=run_sl_cut, run_combinr=run_combinr,
    )


def run_assembly(
    config: AssemblyConfig,
    *,
    force_steps: list[str] | None = None,
) -> None:
    """Run the assembly pipeline with manifest tracking."""
    jaccard_on = config.jaccard_clip
    defuse_on = config.defuse and config.uniprot_db is not None

    # map_transcripts maps <track>.jaccard.fasta and max_intron_split reads
    # <track>.genome.defuse.gff3 whenever the file exists, so one left by an earlier
    # run would shadow the fresh models while its step does not run. Remove those
    # files first; the input fingerprints then re-run the steps that consumed them.
    stale: list[tuple[str, str]] = []
    if not (jaccard_on and config.left_reads and config.right_reads):
        # skipped, or run_jaccard's no-op on single-end reads
        stale += [(jaccard_output_name(query), "jaccard") for query, _ in TRANSCRIPT_SETS]
    if not defuse_on:
        stale += [(f"{stem}{DEFUSE_GFF3}", "defuse") for stem in mapped_transcript_stems()]
    for name, step in stale:
        path = config.work_dir / name
        if path.exists():
            log.info("Removing %s from an earlier run: %s does not run now.", name, step)
            path.unlink()

    run_simple_pipeline(
        ASSEMBLY, STEPS, config,
        force_steps=force_steps,
        skip=lambda s: (s.name == "jaccard" and not jaccard_on)
        or (s.name == "defuse" and not defuse_on),
    )
