"""EVidenceModeler consensus gene model building."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from eukan.infra.runner import run_cmd, run_shell
from eukan.infra.steps import step_dir
from eukan.infra.utils import symlink
from eukan.infra.logging import get_logger
from eukan.settings import PipelineConfig

log = get_logger(__name__)


def run_evm(config: PipelineConfig, evidence: list[Path]) -> Path:
    """Run EVidenceModeler to build consensus gene models."""
    sdir = step_dir(config.work_dir, "evm_consensus_models")
    log.info("Running EVidenceModeler consensus building...")
    run_cmd(["cdbfasta", str(config.genome)], cwd=sdir)

    w = [str(d) for d in config.weights]
    weights_map = {
        "prot.gff3": ["PROTEIN", "prot_align", w[0]],
        "augustus.gff3": ["ABINITIO_PREDICTION", "augustus", w[1]],
        "snap.gff3": ["ABINITIO_PREDICTION", "snap", w[1]],
        "genemark.gff3": ["ABINITIO_PREDICTION", "genemark", w[1]],
        "codingquarry.gff3": ["ABINITIO_PREDICTION", "codingquarry", w[1]],
        "nr_transcripts.gff3": ["TRANSCRIPT", "PASA-assembly", w[2]],
    }

    with open(sdir / "weights.txt", "w") as wf, \
         open(sdir / "gene_predictions.gff3", "w") as pf:
        for ev in evidence:
            ev_name = ev.name
            symlink(ev, sdir / ev_name)
            if ev_name in weights_map:
                wf.write("\t".join(weights_map[ev_name]) + "\n")
            if ev_name not in ("nr_transcripts.gff3", "prot.gff3"):
                pf.write(ev.read_text())

    # Partition and run EVM
    run_cmd(
        [
            "partition_EVM_inputs.pl",
            "--genome", str(config.genome),
            "--gene_predictions", "gene_predictions.gff3",
            "--transcript_alignments", "nr_transcripts.gff3",
            "--protein_alignments", "prot.gff3",
            "--segmentSize", "100000",
            "--overlapSize", "10000",
            "--partition_listing", "partitions_list.out",
            "--partition_dir", str(sdir / "partitions"),
        ],
        cwd=sdir,
    )

    stop_codons = ",".join(config.genetic_code_obj.stop_codons)
    run_cmd(
        [
            "write_EVM_commands.pl",
            "--genome", str(config.genome),
            "--weights", str(sdir / "weights.txt"),
            "--gene_predictions", "gene_predictions.gff3",
            "--protein_alignments", "prot.gff3",
            "--transcript_alignments", "nr_transcripts.gff3",
            "--output_file_name", "consensus_models.out",
            "--stop_codons", stop_codons,
            "--partitions", "partitions_list.out",
        ],
        cwd=sdir,
        out_file="commands.list",
    )

    # Run EVM commands in parallel batches
    # Commands contain shell redirections (> and 2>) so run via shell
    evm_cmds = [
        line.rstrip()
        for line in (sdir / "commands.list").read_text().splitlines()
        if line.strip()
    ]

    def _run_evm_cmd(cmd_str: str) -> None:
        run_shell(cmd_str, cwd=sdir)

    with ThreadPoolExecutor(max_workers=config.num_cpu) as pool:
        list(pool.map(_run_evm_cmd, evm_cmds))

    # Recombine
    run_cmd(
        [
            "recombine_EVM_partial_outputs.pl",
            "--partitions", "partitions_list.out",
            "--output_file_name", "consensus_models.out",
        ],
        cwd=sdir,
    )
    run_cmd(
        [
            "convert_EVM_outputs_to_GFF3.pl",
            "--partitions", "partitions_list.out",
            "--output", "consensus_models.out",
            "--genome", str(config.genome),
        ],
        cwd=sdir,
    )

    # Gather all consensus GFF3 files
    cons_files = sorted(sdir.rglob("consensus_models.out.gff3"))
    with open(sdir / "consensus_models.gff3", "w") as outfile:
        for f in cons_files:
            outfile.write(f.read_text())

    return sdir / "consensus_models.gff3"
