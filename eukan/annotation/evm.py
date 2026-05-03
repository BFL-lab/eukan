"""EVidenceModeler consensus gene model building."""

from __future__ import annotations

import shlex
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from eukan.infra.logging import get_logger
from eukan.infra.runner import run_cmd, run_shell
from eukan.infra.steps import step_dir
from eukan.infra.utils import symlink
from eukan.settings import PipelineConfig

log = get_logger(__name__)


def _parse_evm_command(
    cmd_str: str,
) -> tuple[list[str], Path | None, str | None, str | None] | None:
    """Parse a single line of EVM's commands.list.

    Returns ``(argv, cwd, stdout_file, stderr_file)`` if the line can be
    represented as a plain run_cmd call (handles trailing ``> OUT`` and
    ``2> ERR`` redirects and an optional leading ``cd PATH &&``).  Returns
    ``None`` for anything else, signalling the caller to fall back to a
    shell.
    """
    try:
        tokens = shlex.split(cmd_str)
    except ValueError:
        return None

    if not tokens:
        return None

    cwd: Path | None = None
    # Leading "cd PATH && rest..." → strip and set cwd.
    if len(tokens) >= 3 and tokens[0] == "cd" and tokens[2] == "&&":
        cwd = Path(tokens[1])
        tokens = tokens[3:]

    # Walk from the right end, peeling off `> file` and `2> file`.
    stdout_file: str | None = None
    stderr_file: str | None = None
    while len(tokens) >= 2:
        op = tokens[-2]
        target = tokens[-1]
        if op == ">":
            if stdout_file is not None:
                return None
            stdout_file = target
            tokens = tokens[:-2]
        elif op == "2>":
            if stderr_file is not None:
                return None
            stderr_file = target
            tokens = tokens[:-2]
        else:
            break

    # Reject anything that still contains shell metacharacters we can't handle.
    if any(t in {"&&", "||", "|", ";", "<", ">>", "2>>"} for t in tokens):
        return None
    if not tokens:
        return None

    return tokens, cwd, stdout_file, stderr_file


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

    # Run EVM commands in parallel.
    # Each line in commands.list looks roughly like
    #   evidence_modeler.pl --genome ... > consensus.out 2> consensus.err
    # so we parse out > / 2> redirects and call run_cmd directly to avoid
    # forking a shell per partition (thousands for large genomes).
    evm_cmds = [
        line.rstrip()
        for line in (sdir / "commands.list").read_text().splitlines()
        if line.strip()
    ]

    def _run_evm_cmd(cmd_str: str) -> None:
        parsed = _parse_evm_command(cmd_str)
        if parsed is None:
            # Fallback for any line we can't safely strip-and-tokenize
            run_shell(cmd_str, cwd=sdir)
            return
        argv, run_cwd, stdout_file, stderr_file = parsed
        run_cmd(
            argv,
            cwd=run_cwd or sdir,
            out_file=stdout_file,
            err_file=stderr_file,
        )

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
