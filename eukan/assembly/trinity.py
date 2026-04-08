"""Trinity genome-guided and de novo assembly."""

from __future__ import annotations

import shutil
from pathlib import Path

from eukan.infra.runner import run_cmd
from eukan.infra.utils import step_done
from eukan.infra.logging import get_logger
from eukan.settings import AssemblyConfig

log = get_logger(__name__)


def run_trinity(config: AssemblyConfig, force: bool = False) -> None:
    """Run genome-guided and de novo Trinity assembly."""
    wd = config.work_dir

    if not force and step_done(wd, ["trinity-gg.fasta", "trinity-denovo.fasta"]):
        log.info("[run_trinity] Already complete, skipping. Use force=True to re-run.")
        return

    lib_type_args = (
        ["--SS_lib_type", config.strand_specific] if config.strand_specific else []
    )
    jaccard_args = ["--jaccard_clip"] if config.jaccard_clip else []

    # Genome-guided assembly
    gg_fasta = wd / "trinity-gg.fasta"
    if not gg_fasta.exists():
        log.info("Running genome-guided Trinity assembly...")
        run_cmd(
            [
                "Trinity",
                "--genome_guided_bam", "STAR_Aligned.sortedByCoord.out.bam",
                "--genome_guided_max_intron", str(config.max_intron_len),
                "--max_memory", config.memory_gb,
                "--CPU", str(config.num_cpu),
                "--full_cleanup",
                "--output", "trinity-gg",
                *lib_type_args,
                *jaccard_args,
            ],
            cwd=wd,
        )
        # --full_cleanup puts output at {output_prefix}.Trinity-GG.fasta
        # alongside (not inside) the output dir
        trinity_gg_output = wd / "trinity-gg.Trinity-GG.fasta"
        if not trinity_gg_output.exists():
            # Fallback: check inside the output dir (no --full_cleanup)
            trinity_gg_output = wd / "trinity-gg" / "Trinity-GG.fasta"
        if trinity_gg_output.exists():
            shutil.move(str(trinity_gg_output), str(gg_fasta))
        shutil.rmtree(wd / "trinity-gg", ignore_errors=True)

    # De novo assembly
    dn_fasta = wd / "trinity-denovo.fasta"
    if not dn_fasta.exists():
        log.info("Running de novo Trinity assembly...")
        run_cmd(
            [
                "Trinity",
                "--seqType", "fq",
                "--max_memory", config.memory_gb,
                *config.reads_args_trinity,
                "--CPU", str(config.num_cpu),
                "--output", "trinity-denovo",
                "--full_cleanup",
                *lib_type_args,
                *jaccard_args,
            ],
            cwd=wd,
        )
        # --full_cleanup puts output at {output_prefix}.Trinity.fasta
        trinity_dn_output = wd / "trinity-denovo.Trinity.fasta"
        if not trinity_dn_output.exists():
            # Fallback: check inside the output dir (no --full_cleanup)
            trinity_dn_output = wd / "trinity-denovo" / "Trinity.fasta"
        if trinity_dn_output.exists():
            shutil.move(str(trinity_dn_output), str(dn_fasta))
        shutil.rmtree(wd / "trinity-denovo", ignore_errors=True)
