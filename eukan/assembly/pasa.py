"""PASA spliced alignment and transcript hint generation."""

from __future__ import annotations

import os
import shutil
from pathlib import Path

from eukan.infra.runner import run_cmd
from eukan.infra.utils import step_done
from eukan.infra.logging import get_logger
from eukan.settings import AssemblyConfig

log = get_logger(__name__)


def run_pasa(config: AssemblyConfig, force: bool = False) -> None:
    """Run PASA to assemble spliced alignments from transcriptome assemblies."""
    wd = config.work_dir

    if not force and step_done(wd, ["nr_transcripts.fasta", "nr_transcripts.gff3", "hints_rnaseq.gff"]):
        log.info("[run_pasa] Already complete, skipping. Use force=True to re-run.")
        return

    log.info("Running PASA spliced alignment...")

    pasahome = os.environ.get("PASAHOME", "/opt/PASApipeline.v2.4.1")
    db_path = wd / f"{config.name}.sqlite"

    # Write PASA configs
    with open(wd / "annotCompare.config", "w") as f:
        f.write(f"DATABASE={db_path}\n")
    with open(wd / "alignAssembly.config", "w") as f:
        f.write(f"DATABASE={db_path}\n")
        f.write("validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=95\n")
        f.write("validate_alignments_in_db.dbi:--MIN_AVG_PER_ID=95\n")
        f.write(f"validate_alignments_in_db.dbi:--NUM_BP_PERFECT_SPLICE_BOUNDARY={config.splice_boundary_stringency}\n")
        f.write("subcluster_builder.dbi:-m=50\n")

    # Extract de novo accessions
    with open(wd / "tdn.accs", "w") as f:
        for line in (wd / "trinity-denovo.fasta").read_text().splitlines():
            if line.startswith(">"):
                f.write(line[1:].split()[0] + "\n")

    # Concatenate assemblies
    comprehensive = wd / "trinity-comprehensive.fasta"
    with open(comprehensive, "w") as out:
        for fa in ["trinity-denovo.fasta", "trinity-gg.fasta"]:
            path = wd / fa
            if path.exists():
                out.write(path.read_text())

    # Clean sequences
    run_cmd(
        ["seqclean", "trinity-comprehensive.fasta", "-l", "90", "-c", "6"],
        cwd=wd,
    )

    # Build strand args
    strand_args = ["--transcribed_is_aligned_orient"] if config.strand_specific else []

    # Genetic code args
    gc_args = config.genetic_code_obj.pasa_flag

    # Run PASA
    run_cmd(
        [
            "Launch_PASA_pipeline.pl",
            "-c", "alignAssembly.config",
            "-C", "-r", "-R",
            "-g", str(config.genome),
            "-t", "trinity-comprehensive.fasta.clean",
            "-T",
            "-u", "trinity-comprehensive.fasta",
            "--ALIGNERS", "gmap,blat",
            "--CPU", str(config.num_cpu),
            "--TDN", "tdn.accs",
            "-I", str(config.max_intron_len),
            "--stringent_alignment_overlap", "30.0",
            *strand_args,
            *gc_args,
        ],
        cwd=wd,
    )

    # Build comprehensive transcriptome
    build_compreh = f"{pasahome}/scripts/build_comprehensive_transcriptome.dbi"
    run_cmd(
        [
            build_compreh,
            "-c", "alignAssembly.config",
            "-t", f"{config.name}.sqlite.assemblies.fasta",
            "--min_per_ID", "95",
            "--min_per_aligned", "95",
        ],
        cwd=wd,
    )

    # Process output: deduplicate FASTA and create GFF3 + hints
    compreh_dir = wd / "compreh_init_build"
    if compreh_dir.exists():
        shutil.copy(compreh_dir / "compreh_init_build.fasta", wd)
        shutil.copy(compreh_dir / "compreh_init_build.gff3", wd)

    # Deduplicate FASTA
    _deduplicate_fasta(
        wd / "compreh_init_build.fasta",
        wd / "nr_transcripts.fasta",
    )

    # Report non-redundant transcript count
    nr_path = wd / "nr_transcripts.fasta"
    if nr_path.exists():
        nr_count = sum(1 for line in open(nr_path) if line.startswith(">"))
        if nr_count < 1000:
            log.warning(
                "Only %d non-redundant transcripts assembled "
                "— transcript evidence may be insufficient",
                nr_count,
            )
        else:
            log.info("%d non-redundant transcripts assembled", nr_count)

    # Convert GFF3 to transcript hints
    _build_transcript_hints(wd)


def _deduplicate_fasta(input_fa: Path, output_fa: Path) -> None:
    """Remove duplicate sequences from a FASTA file."""
    seen: set[str] = set()
    with open(input_fa) as fin, open(output_fa, "w") as fout:
        write = False
        for line in fin:
            if line.startswith(">"):
                write = line not in seen
                if write:
                    seen.add(line)
            if write:
                fout.write(line)


def _build_transcript_hints(wd: Path) -> None:
    """Build transcript GFF3 and RNA-seq hints from PASA comprehensive build output."""
    gff3_in = wd / "compreh_init_build.gff3"
    if not gff3_in.exists():
        return

    count = 0
    with open(gff3_in) as fin, \
         open(wd / "nr_transcripts.gff3", "w") as gff_out, \
         open(wd / "hints_transcripts.gff", "w") as hints_out:
        for line in fin:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue

            # Parse the transcript ID from attributes
            attrs = cols[8]
            parts = attrs.split(";")
            transcript_id = ""
            for part in parts:
                if "=" in part:
                    key, val = part.split("=", 1)
                    if key.strip() == "ID":
                        transcript_id = val.strip()
                        break

            count += 1
            cols[1] = "PASA-assembly"
            cols[2] = "exon"
            cols[8] = f"ID={transcript_id}:exon:{count};Parent={transcript_id};"
            gff_out.write("\t".join(cols) + "\n")

            # Also write as hint
            hint_attrs = f"pri=3;src=E;group={transcript_id}"
            hint_cols = cols[:8] + [hint_attrs]
            hints_out.write("\t".join(hint_cols) + "\n")

    # Merge all hints
    with open(wd / "hints_rnaseq.gff", "w") as out:
        for hf in ["hints_transcripts.gff", "hints_introns.gff", "hints_coverage.gff"]:
            path = wd / hf
            if path.exists():
                out.write(path.read_text())
