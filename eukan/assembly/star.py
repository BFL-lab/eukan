"""STAR read mapping and hint generation."""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path

from eukan.exceptions import ExternalToolError
from eukan.infra.logging import get_logger
from eukan.infra.runner import run_cmd
from eukan.settings import AssemblyConfig

log = get_logger(__name__)


def _is_gzipped(path: Path) -> bool:
    """Check if a file is gzip-compressed by reading the magic bytes."""
    try:
        with open(path, "rb") as f:
            return f.read(2) == b"\x1f\x8b"
    except OSError:
        return False


def map_reads(config: AssemblyConfig) -> None:
    """Map RNA-seq reads to the genome using STAR."""
    wd = config.work_dir
    log.info("Running STAR read mapping...")

    index_dir = wd / "build-index"

    # Build genome index
    if not index_dir.exists():
        index_dir.mkdir()
        run_cmd(
            [
                "STAR",
                "--genomeSAindexNbases", "3",
                "--limitGenomeGenerateRAM", "40317074816",
                "--runThreadN", str(config.num_cpu),
                "--runMode", "genomeGenerate",
                "--genomeDir", str(index_dir),
                "--genomeFastaFiles", str(config.genome),
            ],
            cwd=wd,
        )

    # Detect compressed input
    reads = config.reads_args_star
    zcat_args = []
    first_read_file = Path(reads[0])
    if first_read_file.suffix in (".gz", ".gzip") or _is_gzipped(first_read_file):
        zcat_args = ["--readFilesCommand", "zcat"]

    quality_args = ["--outQSconversionAdd", "-31"] if config.phred_quality == 64 else []

    max_intron_args = (
        ["--alignIntronMax", str(config.max_intron_len)]
        if config.max_intron_len
        else []
    )

    star_cmd = [
        "STAR",
        "--runThreadN", str(config.num_cpu),
        "--genomeDir", str(index_dir),
        "--alignEndsType", config.align_mode,
        "--readFilesIn", *reads,
        "--outSAMtype", "BAM", "SortedByCoordinate",
        "--outSJfilterIntronMaxVsReadN", "100", "300", "500",
        "--alignIntronMin", str(config.min_intron_len),
        *max_intron_args,
        "--outFileNamePrefix", "STAR_",
        "--outSAMattributes", "All",
        "--outSAMattrIHstart", "0",
        "--outSAMstrandField", "intronMotif",
        "--limitBAMsortRAM", "27643756136",
        *zcat_args,
        *quality_args,
    ]

    try:
        run_cmd(star_cmd, cwd=wd)
    except ExternalToolError:
        log.warning("STAR failed, falling back to STARlong")
        star_long_cmd = ["STARlong", *star_cmd[1:]]
        run_cmd(star_long_cmd, cwd=wd)

    # Report mapping rate from STAR log
    _log_mapping_rate(wd)

    # Generate hints from splice junctions and analyze splice sites
    _generate_hints_from_star(wd, config.genome)

    # Cleanup build index (large)
    shutil.rmtree(index_dir, ignore_errors=True)


def _log_mapping_rate(wd: Path) -> None:
    """Parse STAR Log.final.out and log the overall mapping rate."""
    log_file = wd / "STAR_Log.final.out"
    if not log_file.exists():
        return

    unique_pct = 0.0
    multi_pct = 0.0
    for line in log_file.read_text().splitlines():
        if "Uniquely mapped reads %" in line:
            unique_pct = float(line.split("|")[-1].strip().rstrip("%"))
        elif "% of reads mapped to multiple loci" in line:
            multi_pct = float(line.split("|")[-1].strip().rstrip("%"))

    total_pct = unique_pct + multi_pct
    if total_pct < 75:
        log.warning(
            "Low read mapping rate: %.1f%% (%.1f%% unique, %.1f%% multi) "
            "— check genome/reads compatibility",
            total_pct, unique_pct, multi_pct,
        )
    else:
        log.info(
            "Read mapping rate: %.1f%% (%.1f%% unique, %.1f%% multi)",
            total_pct, unique_pct, multi_pct,
        )


# STAR motif codes → canonical splice site dinucleotide pairs
_STAR_MOTIF_NAMES: dict[int, str] = {
    1: "GT-AG",
    2: "CT-AC",
    3: "GC-AG",
    4: "CT-GC",
    5: "AT-AC",
    6: "GT-AT",
}


def _analyze_splice_sites(sj_file: Path, genome: Path, wd: Path) -> None:
    """Extract splice site dinucleotides from STAR junctions and write a summary.

    For each junction in SJ.out.tab, extracts the donor and acceptor
    dinucleotides from the genome FASTA.  Writes ``splice_site_summary.json``
    with per-type counts and read support.

    STAR SJ.out.tab columns (from STAR source, OutSJ.cpp):
      col2 = first base of intron (1-based)
      col3 = last base of intron (1-based)
      col5 = motif (0=non-canonical, 1=GT/AG, 2=CT/AC, 3=GC/AG, ...)
      col7 = unique reads, col8 = multi-mapping reads
    """
    from Bio import SeqIO as _SeqIO

    contigs = _SeqIO.to_dict(_SeqIO.parse(str(genome), "fasta"))

    # Tally: splice_type → {"count": int, "unique_reads": int}
    tallies: dict[str, dict[str, int]] = defaultdict(lambda: {"count": 0, "unique_reads": 0})

    with open(sj_file) as fin:
        reader = csv.reader(fin, delimiter="\t")
        for row in reader:
            chrom = row[0]
            intron_start = int(row[1])  # 1-based, first base of intron
            intron_end = int(row[2])    # 1-based, last base of intron
            motif = int(row[4])
            unique = int(row[6])

            if motif != 0:
                # Use STAR's motif classification for canonical/semi-canonical
                splice_type = _STAR_MOTIF_NAMES[motif]
            else:
                # Extract actual dinucleotides from the genome
                seq = contigs.get(chrom)
                if seq is None or intron_end > len(seq):
                    splice_type = "unknown"
                else:
                    genome_seq = seq.seq
                    donor = str(genome_seq[intron_start - 1 : intron_start + 1]).upper()
                    acceptor = str(genome_seq[intron_end - 2 : intron_end]).upper()
                    splice_type = f"{donor}-{acceptor}"

            tallies[splice_type]["count"] += 1
            tallies[splice_type]["unique_reads"] += unique

    summary = dict(sorted(tallies.items(), key=lambda kv: -kv[1]["count"]))
    with open(wd / "splice_site_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    # Log summary
    for stype, counts in summary.items():
        if stype in ("GT-AG", "CT-AC"):
            continue  # skip canonical in log — they dominate
        log.info(
            "Splice sites (%s): %d junctions, %d unique reads",
            stype, counts["count"], counts["unique_reads"],
        )


def _generate_hints_from_star(wd: Path, genome: Path) -> None:
    """Generate AUGUSTUS hints from STAR splice junction and coverage output."""
    # Parse splice junctions into intron hints
    sj_file = wd / "STAR_SJ.out.tab"
    if sj_file.exists():
        _analyze_splice_sites(sj_file, genome, wd)
    if sj_file.exists():
        strand_map = {"0": ".", "1": "+", "2": "-"}
        with open(sj_file) as fin, open(wd / "hints_introns.gff", "w") as fout:
            reader = csv.reader(fin, delimiter="\t")
            for row in reader:
                chrom, start, end = row[0], row[1], row[2]
                strand = strand_map.get(row[3], ".")
                unique = int(row[6]) + int(row[7])
                fout.write(
                    f"{chrom}\tSTAR\tintron\t{start}\t{end}\t{unique}\t"
                    f"{strand}\t.\tmult={unique};pri=4;src=E\n"
                )

    # Generate coverage hints from BAM
    bam = wd / "STAR_Aligned.sortedByCoord.out.bam"
    if bam.exists():
        run_cmd(
            ["samtools", "view", "-b", "-f", "0x10", str(bam)],
            cwd=wd, out_file="STAR_reverse.bam", binary=True,
        )
        run_cmd(
            ["samtools", "view", "-b", "-F", "0x10", str(bam)],
            cwd=wd, out_file="STAR_forward.bam", binary=True,
        )

        for direction, _strand, wig in [
            ("STAR_reverse.bam", "-", "minus.wig"),
            ("STAR_forward.bam", "+", "plus.wig"),
        ]:
            run_cmd(["bam2wig", direction], cwd=wd, out_file=wig)

        # Generate exonpart hints from coverage
        # wig2hints.pl reads from stdin, so we pipe the wig file in
        for wig_file, strand, hints_file in [
            ("minus.wig", "-", "hints.ep.minus.gff"),
            ("plus.wig", "+", "hints.ep.plus.gff"),
        ]:
            _run_wig2hints(wd, wig_file, strand, hints_file)

        # Merge coverage hints
        with open(wd / "hints_coverage.gff", "w") as out:
            for hf in ["hints.ep.minus.gff", "hints.ep.plus.gff"]:
                path = wd / hf
                if path.exists():
                    out.write(path.read_text())

        # Cleanup intermediate files
        for f in ["STAR_reverse.bam", "STAR_forward.bam", "minus.wig", "plus.wig",
                   "hints.ep.minus.gff", "hints.ep.plus.gff"]:
            (wd / f).unlink(missing_ok=True)


def _run_wig2hints(wd: Path, wig_file: str, strand: str, out_file: str) -> None:
    """Run wig2hints.pl with stdin/stdout redirection.

    wig2hints.pl reads from stdin and writes to stdout:
        wig2hints.pl [options] < input.wig > output.gff
    """
    cmd = [
        "wig2hints.pl",
        "--width=10", "--margin=10", "--minthresh=2",
        "--minscore=4", "--prune=0.1", "--src=W",
        "--type=exonpart", "--radius=4.5", "--pri=4",
        f"--strand={strand}",
    ]

    log.debug("Running: %s < %s > %s", " ".join(cmd), wig_file, out_file)

    with open(wd / wig_file) as stdin_f, open(wd / out_file, "w") as stdout_f:
        result = subprocess.run(
            cmd,
            cwd=wd,
            stdin=stdin_f,
            stdout=stdout_f,
            stderr=subprocess.PIPE,
            text=True,
        )

    if result.returncode != 0:
        raise ExternalToolError(
            "wig2hints.pl failed",
            tool="wig2hints.pl",
            returncode=result.returncode,
            cmd=cmd,
            stderr_snippet=result.stderr or "",
        )
