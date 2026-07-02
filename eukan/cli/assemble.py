"""eukan assemble — transcriptome assembly pipeline."""

from __future__ import annotations

from pathlib import Path

import click
from click_option_group import optgroup

from eukan.cli._framework import (
    FULL_CODE_TABLE,
    PreformattedEpilogCommand,
    code_option,
    drop_none,
    force_option,
    genome_option,
    numcpu_option,
    resolve_optional_path,
)


@click.command(cls=PreformattedEpilogCommand, epilog=FULL_CODE_TABLE)
@optgroup.group("Required input")
@genome_option("Genome FASTA file.")
@optgroup.option("--left", "-l", type=click.Path(exists=True, path_type=Path), help="Left paired-end reads.")
@optgroup.option("--right", "-r", type=click.Path(exists=True, path_type=Path), help="Right paired-end reads.")
@optgroup.option("--single", "-s", type=click.Path(exists=True, path_type=Path), help="Single-end reads.")
@optgroup.group("Pipeline parameters")
@numcpu_option
@optgroup.option(
    "--strand-specific", "-S", type=click.Choice(["RF", "FR", "R", "F"]), default=None,
    help="Strand-specific library type.",
)
@optgroup.option(
    "--aligner", type=click.Choice(["auto", "star", "segemehl"]),
    default="auto", show_default=True,
    help="Read aligner. 'auto' (default) maps with STAR, then re-maps with "
    "splice-agnostic segemehl when the diagnostic finds extensive non-canonical "
    "splicing (so StringTie/hints aren't biased by STAR's canonical alignment). "
    "'star' skips that escalation; 'segemehl' maps with segemehl from the start.",
)
@optgroup.option(
    "--align-mode", "-t", type=click.Choice(["EndToEnd", "Local"]),
    default="Local", show_default=True,
    help="STAR read alignment mode (end-to-end vs soft-clipped local). "
    "STAR only; ignored when --aligner segemehl.",
)
@optgroup.option(
    "--splice-permissive", is_flag=True, default=False,
    help="Allow non-canonical splice sites (GC-AG, AT-AC). "
    "Sets splice boundary stringency to 0 and retains non-canonical junctions.",
)
@optgroup.option(
    "--diagnose-softclips/--no-diagnose-softclips", default=True,
    show_default=True,
    help="Run the soft-clip + intron diagnostic after STAR. "
    "Detects trans-splicing (via de novo splice-leader clusters) and "
    "non-canonical splice prevalence; surfaces both as INFO/WARNING.",
)
@code_option(default=1)
@optgroup.option("--min-intron", "-m", type=int, default=20, show_default=True, help="Minimum intron length.")
@optgroup.option(
    "--max-intron", "-M", type=int, default=5000, show_default=True,
    help="Maximum intron length, hard-imposed: transcript models are split at any "
    "longer intron (the max_intron_split step) and combinr enforces it, and Trinity "
    "genome-guided uses it as --genome_guided_max_intron. Changing it on a resumed "
    "run re-runs max_intron_split/combinr automatically; recovering longer introns "
    "from the assembly or mapping also needs --run-trinity/--run-map-transcripts. "
    "Set 0 to disable the model split.",
)
@optgroup.option("--phred", type=click.Choice(["33", "64"]), default="33", show_default=True, help="Phred quality score.")
@optgroup.option(
    "--jaccard-clip/--no-jaccard-clip", "-j", default=None,
    help="In-house jaccard clipping of fused Trinity transcripts (both de novo and "
    "genome-guided), splitting two adjacent loci joined into one contig. On by "
    "default; needs paired reads. --no-jaccard-clip disables; tune via the "
    "--jaccard-* knobs below.",
)
@optgroup.option(
    "--jaccard-greediness", type=float, default=None,
    help="Jaccard low-coverage sensitivity (default 1.5). Coverage-adaptive slack "
    "that lets faint fusion troughs in low-depth regions be split: 0 = Trinity's "
    "fixed floor (most stringent), higher = more aggressive in low coverage. "
    "Only used with -j.",
)
@optgroup.option(
    "--jaccard-max-trough", type=float, default=None,
    help="Jaccard trough-depth gate (default 0.05): a junction's jaccard must dip to "
    "<= this to be cut. LOWER for more stringent clipping (a deeper trough required). "
    "Only used with -j.",
)
@optgroup.option(
    "--jaccard-min-delta", type=float, default=None,
    help="Jaccard flanking-hill rise (default 0.35): the jaccard must climb this far "
    "above the trough on both sides for a cut. RAISE for more stringent clipping. "
    "Only used with -j.",
)
@optgroup.option(
    "--jaccard-max-adaptive-trough", type=float, default=None,
    help="Ceiling on the coverage-adaptive trough gate (default 0.30): caps how far "
    "--jaccard-greediness can relax --jaccard-max-trough in low coverage. LOWER to "
    "keep low-coverage clipping stringent. Only used with -j.",
)
@optgroup.option(
    "--combinr-stringent-overlap", type=float, default=None,
    help="combinr --stringent-overlap: two "
    "transcripts cluster only when their span overlap is >= this percent of the "
    "shorter transcript (default 0 = any overlap). Raise it (e.g. 30) to stop short "
    "or tip-overlapping transcripts welding collinear neighbours into one model.",
)
@optgroup.option(
    "--sl-sequence", default=None,
    help="Override the spliced-leader sequence for trans-splice acceptor detection "
    "(else recovered from the read soft-clip verdict or de novo insertions).",
)
@optgroup.option(
    "--sl-cluster-window", type=int, default=5, show_default=True,
    help="Genomic window (bp) for consolidating SL acceptor sites.",
)
@optgroup.option(
    "--adapter-sequence", "adapter_sequence", multiple=True,
    help="Extra sequencing-adapter sequence to exclude from SL detection "
    "(repeatable). Added to the built-in Illumina/Nextera set; use for "
    "platform-specific adapters (e.g. MGI/BGI).",
)
@optgroup.option(
    "--sl-adapter-filter/--no-sl-adapter-filter", default=None,
    help="Screen the recovered SL consensus / acceptor candidates against known "
    "sequencing adapters so residual Illumina read-through (e.g. AGATCGGAAGAGC) "
    "isn't mistaken for a spliced leader. On by default; --no-sl-adapter-filter "
    "disables it.",
)
@optgroup.option(
    "--combinr-path", type=click.Path(exists=True, path_type=Path), default=None,
    help="Path to the combinr binary (else resolved from PATH or "
    "$EUKAN_ASSEMBLE_COMBINR_PATH).",
)
@optgroup.option(
    "--uniprot", type=click.Path(exists=True, path_type=Path), default=None,
    help="SwissProt FASTA (uniprot_sprot.faa) or prebuilt diamond .dmnd. Enables "
    "homology-based splice-strand correction on unstranded libraries (skipped "
    "when -S/--strand-specific is given), and homology de-fusion with --defuse. "
    "Without it, both are no-ops.",
)
@optgroup.option(
    "--defuse", is_flag=True, default=False,
    help="Split chimeric (fused) transcripts using protein homology: an "
    "ultra-sensitive diamond blastx vs SwissProt that finds >=2 distinct, "
    "non-overlapping hits on one transcript cuts it at the inter-hit gap. "
    "Requires --uniprot.",
)
@optgroup.option(
    "--defuse-overlap-tolerance", type=float, default=None,
    help="Max fractional query overlap (of the shorter hit) for two protein hits to "
    "still count as distinct evidence of separate genes (default 0.10).",
)
@optgroup.option(
    "--memory-gb", type=int, default=None,
    help="Assembly memory cap in GiB (Trinity --max_memory / rnaSPAdes -m). "
         "Defaults to 60 percent of currently-available memory (floored at 4 GiB).",
)
@optgroup.group("Re-run steps")
@optgroup.option("--run-star", "-A", is_flag=True, help="Force re-run STAR read mapping.")
@optgroup.option("--run-segemehl", is_flag=True, help="Force re-run segemehl read mapping.")
@optgroup.option(
    "--run-trinity", "-T", is_flag=True,
    help="Force re-run Trinity de novo + genome-guided assembly.",
)
@optgroup.option(
    "--run-jaccard", is_flag=True,
    help="Force re-run jaccard clipping of fused transcripts.",
)
@optgroup.option(
    "--run-map-transcripts", is_flag=True,
    help="Force re-run STAR transcript→genome mapping.",
)
@optgroup.option(
    "--run-strand-correct", is_flag=True,
    help="Force re-run homology-based splice-strand correction.",
)
@optgroup.option(
    "--run-defuse", is_flag=True,
    help="Force re-run homology-based transcript de-fusion.",
)
@optgroup.option(
    "--run-max-intron-split", is_flag=True,
    help="Force re-run the max-intron split of transcript models.",
)
@optgroup.option(
    "--run-sl-detect", is_flag=True,
    help="Force re-run SL trans-splice acceptor detection.",
)
@optgroup.option(
    "--run-sl-cut", is_flag=True,
    help="Force re-run the genomic SL cut of transcript models.",
)
@optgroup.option(
    "--run-combinr", is_flag=True,
    help="Force re-run combinr transcript consolidation.",
)
@force_option
def assemble(
    genome: Path,
    left: Path | None,
    right: Path | None,
    single: Path | None,
    min_intron: int,
    max_intron: int,
    phred: str,
    numcpu: int,
    strand_specific: str | None,
    aligner: str,
    align_mode: str,
    run_star: bool,
    run_segemehl: bool,
    run_trinity: bool,
    run_jaccard: bool,
    run_map_transcripts: bool,
    run_strand_correct: bool,
    run_defuse: bool,
    run_max_intron_split: bool,
    run_sl_detect: bool,
    run_sl_cut: bool,
    run_combinr: bool,
    jaccard_clip: bool | None,
    jaccard_greediness: float | None,
    jaccard_max_trough: float | None,
    jaccard_min_delta: float | None,
    jaccard_max_adaptive_trough: float | None,
    combinr_stringent_overlap: float | None,
    defuse: bool,
    defuse_overlap_tolerance: float | None,
    sl_sequence: str | None,
    sl_cluster_window: int,
    adapter_sequence: tuple[str, ...],
    sl_adapter_filter: bool | None,
    combinr_path: Path | None,
    uniprot: Path | None,
    splice_permissive: bool,
    diagnose_softclips: bool,
    code: int,
    memory_gb: int | None,
    force: bool,
) -> None:
    """Assemble transcriptome from RNA-seq reads.

    \b
    Provide either paired-end reads (--left and --right together) or
    single-end reads (--single). If using paired-end reads, both --left
    and --right are required.
    """
    from eukan.assembly import run_assembly
    from eukan.assembly.pipeline import force_steps_from_run_flags
    from eukan.infra.layout import step_work_dir
    from eukan.settings import AssemblyConfig

    if not left and not right and not single:
        raise click.UsageError("Provide --left/--right (paired) or --single reads.")
    if (left or right) and not (left and right):
        raise click.UsageError("Paired-end mode requires both --left and --right.")

    if strand_specific:
        if single and strand_specific in ("RF", "FR"):
            raise click.UsageError(
                "Paired-end strand types (RF/FR) cannot be used with single-end reads."
            )
        if (left or right) and strand_specific in ("R", "F"):
            raise click.UsageError(
                "Single-end strand types (R/F) cannot be used with paired-end reads."
            )

    if memory_gb is not None and memory_gb < 1:
        raise click.UsageError("--memory-gb must be at least 1 GiB.")

    if defuse and uniprot is None:
        raise click.UsageError("--defuse needs a protein DB; supply --uniprot.")

    config = AssemblyConfig(**drop_none(
        genome=genome.resolve(),
        work_dir=step_work_dir("assemble"),
        manifest_dir=Path.cwd(),
        min_intron_len=min_intron,
        max_intron_len=max_intron,
        phred_quality=int(phred),
        num_cpu=numcpu,
        aligner=aligner,
        align_mode=align_mode,
        jaccard_clip=jaccard_clip,
        jaccard_greediness=jaccard_greediness,
        jaccard_max_trough=jaccard_max_trough,
        jaccard_min_delta=jaccard_min_delta,
        jaccard_max_adaptive_trough=jaccard_max_adaptive_trough,
        combinr_stringent_overlap=combinr_stringent_overlap,
        defuse=defuse,
        defuse_overlap_tolerance=defuse_overlap_tolerance,
        sl_sequence=sl_sequence,
        sl_cluster_window=sl_cluster_window,
        adapter_sequences=list(adapter_sequence),
        sl_adapter_filter=sl_adapter_filter,
        combinr_path=resolve_optional_path(combinr_path),
        uniprot_db=resolve_optional_path(uniprot),
        splice_permissive=splice_permissive,
        diagnose_softclips=diagnose_softclips,
        genetic_code=str(code),
        left_reads=resolve_optional_path(left),
        right_reads=resolve_optional_path(right),
        single_reads=resolve_optional_path(single),
        strand_specific=strand_specific,
        memory_gb=memory_gb,
    ))

    force_steps = force_steps_from_run_flags(
        aligner=aligner,
        run_star=run_star, run_segemehl=run_segemehl,
        run_trinity=run_trinity, run_jaccard=run_jaccard,
        run_map_transcripts=run_map_transcripts,
        run_strand_correct=run_strand_correct, run_defuse=run_defuse,
        run_max_intron_split=run_max_intron_split,
        run_sl_detect=run_sl_detect, run_sl_cut=run_sl_cut,
        run_combinr=run_combinr, force=force,
    )
    run_assembly(config, force_steps=force_steps or None)
    click.echo("Done.")
