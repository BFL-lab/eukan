"""Genome-guided transcript assembly with StringTie.

StringTie assembles transcript models directly from the read→genome BAM
(``config.aligner_bam``), emitting a genomic GTF. It replaces Trinity's
genome-guided mode: the output is already in genome coordinates, so it needs no
transcript→genome mapping. The library is treated as unstranded (no
``--rf``/``--fr``); the SL-cut step (:mod:`eukan.assembly.sl_cut`) later imposes
strand at trans-splice acceptor sites for any ``.``-strand models.
"""

from __future__ import annotations

from eukan.assembly.bam_introns import split_long_introns
from eukan.infra.logging import get_logger
from eukan.infra.runner import run_cmd
from eukan.settings import AssemblyConfig

log = get_logger(__name__)

_GTF = "stringtie.gtf"
# A max-intron-bounded copy of the segemehl read BAM, fed to StringTie so it can't
# fuse distant loci across an over-long intron. Disposable: regenerated each run.
_BOUNDED_BAM = "stringtie_input.bam"


def run_stringtie(config: AssemblyConfig) -> None:
    """Genome-guided assembly off the aligner BAM → ``stringtie.gtf``.

    Skips when ``stringtie.gtf`` already exists (mirrors the rnaSPAdes
    idempotence guard). The coordinate-sorted aligner BAM is read directly; no
    BAM index is required by StringTie.

    On the segemehl path (STAR is already ``--alignIntronMax``-bounded) the read
    BAM can carry over-long introns, so StringTie reads a max-intron-bounded copy
    (:func:`eukan.assembly.bam_introns.split_long_introns`) instead — splitting
    those reads stops StringTie fusing distant loci. The shared aligner BAM is
    left untouched: SL acceptor detection and the non-canonical diagnostic read it
    and must see the true alignments.
    """
    wd = config.work_dir
    out = wd / _GTF
    if out.exists():
        return

    bam = wd / config.aligner_bam
    if config.aligner_bam.startswith("segemehl") and config.max_intron_len:
        bounded = wd / _BOUNDED_BAM
        n_split = split_long_introns(
            bam, bounded, max_intron_len=config.max_intron_len, num_cpu=config.num_cpu
        )
        log.info(
            "Bounded read BAM for StringTie: split %d read(s) at introns > %d nt.",
            n_split, config.max_intron_len,
        )
        bam = bounded

    log.info("Running StringTie genome-guided assembly...")
    run_cmd(
        [
            "stringtie",
            str(bam),
            "-p", str(config.num_cpu),
            # Stringency above StringTie's defaults (-c 1, -f 0.01): drop
            # low-coverage spurious models and minor noise isoforms so the
            # genome-guided set fed to combinr is cleaner. -j (default 1) is
            # exposed so single-read spurious junctions can be filtered too.
            "-c", str(config.stringtie_min_coverage),
            "-f", str(config.stringtie_min_isoform_fraction),
            "-j", str(config.stringtie_min_junction_coverage),
            "-o", _GTF,
        ],
        cwd=wd,
    )
    (wd / _BOUNDED_BAM).unlink(missing_ok=True)
