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

from dataclasses import dataclass, field
from pathlib import Path

from Bio.Seq import Seq

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


@dataclass
class _Transcript:
    """A combinr transcript model: an mRNA id with its genomic exon blocks."""

    tid: str
    chrom: str
    strand: str
    exons: list[tuple[int, int]] = field(default_factory=list)  # (lend, rend), 1-based

    @property
    def start(self) -> int:
        return self.exons[0][0]

    @property
    def end(self) -> int:
        return self.exons[-1][1]


def _combinr_bin(config: AssemblyConfig) -> str:
    """The combinr executable: explicit ``combinr_path`` or ``combinr`` on PATH."""
    return str(config.combinr_path) if config.combinr_path else "combinr"


def _parse_attrs(col9: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for part in col9.split(";"):
        part = part.strip()
        if "=" in part:
            key, val = part.split("=", 1)
            attrs[key.strip()] = val.strip()
    return attrs


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


def _parse_combinr_gff3(path: Path) -> list[_Transcript]:
    """Parse combinr's gene/mRNA/exon GFF3 into transcript records (exon-sorted)."""
    by_id: dict[str, _Transcript] = {}
    order: list[str] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            attrs = _parse_attrs(cols[8])
            if cols[2] == "mRNA":
                tid = attrs.get("ID", "")
                if tid:
                    by_id[tid] = _Transcript(tid, cols[0], cols[6])
                    order.append(tid)
            elif cols[2] == "exon":
                tx = by_id.get(attrs.get("Parent", ""))
                if tx is not None:
                    tx.exons.append((int(cols[3]), int(cols[4])))
    result = []
    for tid in order:
        tx = by_id[tid]
        if tx.exons:
            tx.exons.sort()
            result.append(tx)
    return result


def _write_transcript_alignments_and_hints(
    transcripts: list[_Transcript], gff_out: Path, hints_out: Path
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


def _write_transcript_fasta(transcripts: list[_Transcript], genome: Path, out: Path) -> None:
    """Write spliced transcript sequences (exons joined, RC on '-') from the genome."""
    ordered = sorted(transcripts, key=lambda t: (t.chrom, t.start))  # keep ContigIndex warm
    with ContigIndex(genome) as contigs, open(out, "w") as fh:
        for tx in ordered:
            seq = "".join(str(contigs[tx.chrom][lend - 1 : rend].seq) for lend, rend in tx.exons)
            if tx.strand == "-":
                seq = str(Seq(seq).reverse_complement())
            fh.write(f">{tx.tid}\n{seq}\n")


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
    transcripts = _parse_combinr_gff3(all_gff)
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
