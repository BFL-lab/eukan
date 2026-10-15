"""combinr transcript consolidation.

Runs the external ``combinr assemble`` over the **cut** transcript models
(``{stem}.cut.gff3`` from :mod:`eukan.assembly.sl_cut` — max-intron-split and,
when an SL signal is present, SL-cut de novo transcript→genome GFF3s) to build a
non-redundant set of transcript models, then emits the artifacts the annotation
pipeline consumes as transcript evidence:

* ``nr_transcripts.gff3`` — flat ``exon`` features grouped by ``Parent`` (the
  ``--transcript-alignments`` contract), source ``combinr-assembly``;
* ``nr_transcripts.fasta`` — spliced transcript sequences read off the genome;
* ``hints_rnaseq.gff`` — transcript exon hints concatenated with the intron and
  coverage hints already written by the aligner step.

All inputs are already in genome coordinates, jaccard-clipped, and cut, so
combinr simply consolidates them in one pass; ``combinr assemble`` auto-detects
each input's format (GFF3/GTF/BAM) by extension.
"""

from __future__ import annotations

from pathlib import Path

from eukan.gff.transcripts import Transcript, read_transcripts, stitch
from eukan.infra.artifacts import Artifact
from eukan.infra.genome import ContigIndex
from eukan.infra.logging import get_logger
from eukan.infra.runner import run_cmd
from eukan.infra.utils import concat_files
from eukan.settings import AssemblyConfig

log = get_logger(__name__)

# Cut transcript models from the sl_cut step (genome coordinates), one per mapped
# Trinity track (de novo + genome-guided). Named generically (`.cut.gff3`) because
# the file is a pass-through copy of the max-intron-split models when no SL signal
# is present, so it is not always SL-cut.
_CUT_MODELS = (
    "trinity-denovo.genome.cut.gff3",
    "trinity-gg.genome.cut.gff3",
)
# Source token written into nr_transcripts.gff3; the consensus weights.txt picks
# this up as the TRANSCRIPT evidence source (eukan/annotation/evidence.py::_first_source_token).
_SOURCE = "combinr-assembly"


def _combinr_bin(config: AssemblyConfig) -> str:
    """The combinr executable: explicit ``combinr_path`` or ``combinr`` on PATH."""
    return str(config.combinr_path) if config.combinr_path else "combinr"


def _run_combinr_assemble(config: AssemblyConfig, inputs: list[Path], out_gff: Path) -> None:
    """Run ``combinr assemble`` over *inputs* (GFF3/GTF/BAM), GFF3 stdout to *out_gff*."""
    cmd = [_combinr_bin(config), "assemble"]
    for path in inputs:
        cmd += ["-i", str(path)]
    cmd += [
        "--format", "gff3",
        "-t", str(config.num_cpu),
        "--max-intron", str(config.max_intron_len),
        # --stringent-overlap: 0 = off (any-overlap, the default).
        "--stringent-overlap", str(config.combinr_stringent_overlap),
    ]
    run_cmd(cmd, cwd=config.work_dir, out_file=out_gff.name)


def _write_transcript_alignments_and_hints(
    transcripts: list[Transcript], gff_out: Path, hints_out: Path
) -> None:
    """Emit the transcript-alignments GFF3 + AUGUSTUS exon hints.

    Both are flat ``exon`` features grouped by ``Parent`` (the
    ``--transcript-alignments`` contract), tagged with the ``combinr-assembly`` source.
    """
    n = 0
    with open(gff_out, "w") as gff, open(hints_out, "w") as hints:
        for tx in transcripts:
            for lend, rend in tx.exons:
                n += 1
                base = f"{tx.chrom}\t{_SOURCE}\texon\t{lend}\t{rend}\t.\t{tx.strand}\t."
                gff.write(f"{base}\tID={tx.tid}:exon:{n};Parent={tx.tid}\n")
                hints.write(f"{base}\tpri=3;src=E;group={tx.tid}\n")


def _write_transcript_fasta(transcripts: list[Transcript], genome: Path, out: Path) -> None:
    """Write spliced transcript sequences (exons joined, RC on '-') from the genome."""
    ordered = sorted(transcripts, key=lambda t: (t.chrom, t.start))  # keep ContigIndex warm
    with ContigIndex(genome) as contigs, open(out, "w") as fh:
        for tx in ordered:
            fh.write(f">{tx.tid}\n{stitch(tx, contigs)}\n")


def run_combinr(config: AssemblyConfig) -> None:
    """Consolidate the SL-cut transcript models with combinr."""
    wd = config.work_dir
    inputs = [
        wd / f for f in _CUT_MODELS if (wd / f).exists() and (wd / f).stat().st_size > 0
    ]
    if not inputs:
        raise FileNotFoundError(
            "No cut transcript models found for combinr; run the sl_cut step first."
        )

    all_gff = wd / "combinr_all.gff3"
    _run_combinr_assemble(config, inputs, all_gff)
    transcripts = read_transcripts(all_gff)
    log.info(
        "combinr consolidated %d transcripts from %d input(s): %s.",
        len(transcripts), len(inputs), ", ".join(p.name for p in inputs),
    )

    if not transcripts:
        log.warning("combinr produced no transcript models — transcript evidence is empty.")

    _write_transcript_alignments_and_hints(
        transcripts, wd / Artifact.NR_TRANSCRIPTS_GFF, wd / "hints_transcripts.gff"
    )
    _write_transcript_fasta(transcripts, config.genome, wd / Artifact.NR_TRANSCRIPTS_FASTA)
    concat_files(
        [
            wd / hf
            for hf in ("hints_transcripts.gff", "hints_introns.gff", "hints_coverage.gff")
            if (wd / hf).exists()
        ],
        wd / Artifact.RNASEQ_HINTS,
    )
