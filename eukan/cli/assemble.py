"""eukan assemble — transcriptome assembly pipeline."""

from __future__ import annotations

from pathlib import Path
from typing import Any

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
)
from eukan.infra.layout import step_work_dir


@click.command(cls=PreformattedEpilogCommand, epilog=FULL_CODE_TABLE)
@optgroup.group("Required input")
@genome_option("Genome FASTA file.")
@optgroup.option("--left", "-l", "left_reads", type=click.Path(exists=True, resolve_path=True, path_type=Path), help="Left paired-end reads.")
@optgroup.option("--right", "-r", "right_reads", type=click.Path(exists=True, resolve_path=True, path_type=Path), help="Right paired-end reads.")
@optgroup.option("--single", "-s", "single_reads", type=click.Path(exists=True, resolve_path=True, path_type=Path), help="Single-end reads.")
@optgroup.group("Pipeline parameters")
@numcpu_option
@optgroup.option(
    "--strand-specific", "-S", type=click.Choice(["RF", "FR", "R", "F"]), default=None,
    help="Strand-specific library type.",
)
@optgroup.option(
    "--non-canonical", "non_canonical", type=click.Choice(["auto", "force", "off"]),
    default="auto", show_default=True,
    help="Non-canonical splice mapping. minimap2 maps with the splice preset; "
    "'auto' (default) re-maps with the non-canonical flags "
    "(-J 0 -C 3 --splice-flank=no) when the soft-clip diagnostic finds extensive "
    "non-canonical splicing (so hints/assembly aren't biased by the canonical "
    "alignment). 'force' always applies them; 'off' never does.",
)
@optgroup.option(
    "--diagnose-softclips/--no-diagnose-softclips", default=True,
    show_default=True,
    help="Run the soft-clip + intron diagnostic after read mapping. "
    "Detects trans-splicing (via de novo splice-leader clusters) and "
    "non-canonical splice prevalence; surfaces both as INFO/WARNING and gates "
    "the --non-canonical auto re-map.",
)
@code_option(default=1)
@optgroup.option("--min-intron", "-m", "min_intron_len", type=int, default=20, show_default=True, help="Minimum intron length.")
@optgroup.option(
    "--max-intron", "-M", "max_intron_len", type=int, default=5000, show_default=True,
    help="Maximum intron length, hard-imposed: minimap2 maps with it (-G), Trinity "
    "genome-guided uses it as --genome_guided_max_intron, transcript models are split "
    "at any longer intron (the max_intron_split step) and combinr enforces it. "
    "Changing it on a resumed run re-runs the read and transcript mapping, the Trinity "
    "assembly, the split and combinr automatically. Set 0 to disable the model split.",
)
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
    "--adapter-sequence", "adapter_sequences", multiple=True,
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
    "--combinr-path", type=click.Path(exists=True, resolve_path=True, path_type=Path), default=None,
    help="Path to the combinr binary (else resolved from PATH or "
    "$EUKAN_ASSEMBLE_COMBINR_PATH).",
)
@optgroup.option(
    "--uniprot", "uniprot_db", type=click.Path(exists=True, resolve_path=True, path_type=Path), default=None,
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
    help="Assembly memory cap in GiB (Trinity --max_memory). "
         "Defaults to 60 percent of currently-available memory (floored at 4 GiB).",
)
@optgroup.group("Re-run steps")
@optgroup.option(
    "--run-minimap2", "-A", is_flag=True, help="Force re-run minimap2 read mapping."
)
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
    help="Force re-run minimap2 transcript→genome mapping.",
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
def assemble(**opts: Any) -> None:
    """Assemble transcriptome from RNA-seq reads.

    \b
    Provide either paired-end reads (--left and --right together) or
    single-end reads (--single). If using paired-end reads, both --left
    and --right are required.
    """
    from eukan.assembly import run_assembly  # lazy: startup cost
    from eukan.assembly.pipeline import force_steps_from_run_flags  # lazy: startup cost
    from eukan.settings import AssemblyConfig  # lazy: startup cost

    if not opts["left_reads"] and not opts["right_reads"] and not opts["single_reads"]:
        raise click.UsageError("Provide --left/--right (paired) or --single reads.")

    if opts["memory_gb"] is not None and opts["memory_gb"] < 1:
        raise click.UsageError("--memory-gb must be at least 1 GiB.")

    if opts["defuse"] and opts["uniprot_db"] is None:
        raise click.UsageError("--defuse needs a protein DB; supply --uniprot.")

    run_flags = {k: opts.pop(k) for k in list(opts) if k.startswith("run_") or k == "force"}
    config = AssemblyConfig(**drop_none(
        **opts, work_dir=step_work_dir("assemble"), manifest_dir=Path.cwd(),
    ))

    force_steps = force_steps_from_run_flags(**run_flags)
    run_assembly(config, force_steps=force_steps or None)
    click.echo("Done.")
