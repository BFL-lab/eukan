"""Spliced protein alignment via spaln (intron-rich) or GenomeThreader (intron-poor)."""

from __future__ import annotations

import re
import statistics
from collections import Counter
from pathlib import Path

import gffutils

from eukan.annotation.spaln_params import build_ssp
from eukan.gff.normalize import normalize_to_gff3
from eukan.gff.rows import iter_gff3_rows
from eukan.infra.logging import get_logger
from eukan.infra.runner import run_cmd
from eukan.infra.steps import step_dir
from eukan.infra.utils import concat_files, symlink
from eukan.settings import PipelineConfig
from eukan.validation import validate_fasta

log = get_logger(__name__)


def align_proteins(
    config: PipelineConfig,
    gff3: Path,
    proteins: list[Path],
    intron_hints: Path | None = None,
) -> Path:
    """Spliced alignment of protein sequences to the genome."""
    output = "prot.gff3"
    step_name = "prot_align_ssp" if config.spaln_ssp else "prot_align"
    sdir = step_dir(config.work_dir, step_name)
    log.info("Running protein alignment...")

    concat_files(proteins, sdir / "prots.faa")

    validate_fasta(sdir / "prots.faa")
    symlink(config.genome, sdir / "genome")
    symlink(config.genome, sdir / "genome.gf")

    # Decide between spaln (intron-rich) and gth (intron-poor). We count
    # `CDS` rows rather than `exon` because genemark.gtf — the most common
    # input here — emits CDS/intron/gene/mRNA but no `exon` rows; using
    # `exon` would always evaluate to 0 and force gth.
    lengths, intron_count = _gene_lengths_and_introns(gff3)
    if lengths and intron_count / len(lengths) > 0.25:
        _run_spaln(config, sdir, lengths, intron_hints or gff3)
    else:
        _run_gth(sdir)

    _write_protein_hints(sdir / output, sdir / "hints_protein.gff")
    return sdir / output


def _write_protein_hints(prot_gff3: Path, out: Path) -> None:
    """AUGUSTUS protein hints: each CDS as a CDSpart grouped by its Parent."""
    with open(prot_gff3) as fin, open(out, "w") as fout:
        for line in fin:
            cols = line.rstrip("\n").split("\t")
            if len(cols) == 9 and cols[2] == "CDS" and not line.startswith("#"):
                end = ";" if cols[8].endswith(";") else ""
                cols[2] = "CDSpart"
                cols[8] = f"pri=1;src=P;group={_attr(cols[8], 'Parent')}{end}"
                fout.write("\t".join(cols) + "\n")


def _attr(column9: str, key: str) -> str:
    """Raw (still URL-escaped) value of *key* in a GFF3 attribute column."""
    m = re.search(rf"(?:^|;){key}=([^;]*)", column9)
    return m.group(1) if m else ""


def _gene_lengths_and_introns(gff3: Path) -> tuple[list[int], int]:
    """Gene lengths and intron count of *gff3*, by a line scan.

    A gene's introns are its CDS (under the gene or its mRNAs) minus one.
    """
    lengths: list[int] = []
    gene_of: dict[str, str] = {}
    cds_per_parent: Counter[str] = Counter()
    with open(gff3) as fh:
        for line in fh:
            cols = line.split("\t")
            if len(cols) != 9 or line.startswith("#"):
                continue
            if cols[2] == "gene":
                lengths.append(abs(int(cols[4]) - int(cols[3])))
            elif cols[2] == "mRNA":
                gene_of[_attr(cols[8], "ID")] = _attr(cols[8], "Parent")
            elif cols[2] == "CDS":
                cds_per_parent[_attr(cols[8], "Parent")] += 1
    cds_per_gene: Counter[str] = Counter()
    for parent, n in cds_per_parent.items():
        cds_per_gene[gene_of.get(parent, parent)] += n
    return lengths, sum(n - 1 for n in cds_per_gene.values())


def _run_spaln(
    config: PipelineConfig,
    sdir: Path,
    lengths: list[int],
    intron_hints_path: Path,
) -> None:
    """Run spaln protein alignment (for intron-rich genomes)."""
    mean_len = int(sum(lengths) / len(lengths))
    sd_len = int(statistics.stdev(lengths)) if len(lengths) >= 2 else 0
    maxlen = mean_len + 2 * sd_len

    run_cmd(["makdbs", "-KD", "genome"], cwd=sdir)
    run_cmd(["makblk.pl", "-Wgenome.bkp", "genome.gf"], cwd=sdir)

    use_ssp = config.spaln_ssp and config.transcripts_fasta
    if use_ssp:
        ssp_name = build_ssp(config, sdir)
    else:
        int_dist = _run_fitild(intron_hints_path, sdir)

    run_cmd(
        ["spaln", "-Wgenome.bkp", f"-XG{maxlen}", "-KP", "genome.gf"],
        cwd=sdir,
    )

    # Remove stale spaln output files to avoid interactive overwrite prompts
    for ext in (".grd", ".erd", ".qrd"):
        (sdir / f"prots{ext}").unlink(missing_ok=True)

    if use_ssp:
        run_cmd(
            [
                "spaln", "-O12", "-Q7", "-LS",
                f"-C{config.genetic_code}", f"-t{config.num_cpu}",
                f"-T{ssp_name}", "-dgenome", "prots.faa",
            ],
            cwd=sdir,
        )
    else:
        run_cmd(
            [
                "spaln", "-O12", "-Q7", "-yX", "-LS",
                f"-C{config.genetic_code}", f"-t{config.num_cpu}",
                f"-yI{int_dist}", "-dgenome", "prots.faa",
            ],
            cwd=sdir,
        )

    run_cmd(
        ["sortgrcd", "-F0", "-C60", "-O0", "prots.grd"],
        cwd=sdir, out_file="temp.gff3",
    )

    normalize_to_gff3(
        sdir / "temp.gff3", sdir / "prot.gff3",
        parse_transform=_spaln_fix_cds_featuretype,
        post_transform=_spaln_fix_ids,
    )


def _run_fitild(intron_hints_path: Path, sdir: Path) -> str:
    """Compute intron length distribution via fitild. Returns the -yI parameter string."""
    counts: Counter[int] = Counter()
    for cols in iter_gff3_rows(intron_hints_path):
        try:
            length = int(cols[4]) - int(cols[3])
        except ValueError:
            continue
        counts[length] += 1
    with open(sdir / "introns.ild", "w") as fout:
        for length in sorted(counts):
            fout.write(f"{length} {counts[length]}\n")

    run_cmd(["fitild", "-d", "IldModel.txt", "introns.ild"], cwd=sdir, out_file="fitild.out")

    int_dist_parts: list[str] = []
    with open(sdir / "fitild.out") as fp:
        for line in fp:
            if re.search(r"^introns\s", line):
                parts = line.replace("\t", " ").replace("\n", " ").split()[3:-3]
                int_dist_parts = [f"{float(x):f}" for x in parts]
    if int_dist_parts:
        del int_dist_parts[1]
    return " ".join(int_dist_parts)


def _run_gth(sdir: Path) -> None:
    """Run GenomeThreader protein alignment (for intron-poor genomes)."""
    run_cmd(
        [
            "gth", "-genomic", "genome.gf", "-protein", "prots.faa",
            "-intermediate", "-xmlout", "-gzip", "-o", "gth.xml.gz", "-force",
        ],
        cwd=sdir,
    )
    run_cmd(
        [
            "gthconsensus", "-mincoverage", "0.70", "-minalignmentscore", "0.30",
            "-intermediate", "-gff3out", "-force", "-o", "gth.gff3", "gth.xml.gz",
        ],
        cwd=sdir,
    )
    normalize_to_gff3(
        sdir / "gth.gff3", sdir / "prot.gff3",
        parse_transform=_gth_filter,
        post_transform=_gth_fix_ids,
    )


# ---------------------------------------------------------------------------
# spaln output transforms
# ---------------------------------------------------------------------------


def _spaln_fix_cds_featuretype(f: gffutils.Feature) -> gffutils.Feature:
    """Fix spaln's lowercase ``cds`` featuretype."""
    if f.featuretype == "cds":
        f.attributes["ID"][0] = f"{f.attributes['Parent'][0]}:{f.attributes['ID'][0]}"
        f.featuretype = "CDS"
    return f


def _spaln_fix_ids(f: gffutils.Feature) -> gffutils.Feature:
    """Normalize spaln feature IDs and source."""
    if f.featuretype in ("CDS", "exon"):
        f.attributes["ID"][0] = f.id
    f.source = "prot_align"
    return f


# ---------------------------------------------------------------------------
# GenomeThreader output transforms
# ---------------------------------------------------------------------------


def _gth_filter(f: gffutils.Feature) -> gffutils.Feature | None:
    """Keep only gene and exon features; reparent exons to ``_mRNA``."""
    if f.featuretype in ("gene", "exon"):
        if f.featuretype == "exon":
            f.attributes["Parent"][0] += "_mRNA"
        return f
    return None


def _gth_fix_ids(f: gffutils.Feature) -> gffutils.Feature:
    """Normalize GenomeThreader feature IDs."""
    f.source = "prot_align"
    if f.featuretype == "exon":
        f.attributes["ID"] = [f.id]
    if f.featuretype == "mRNA":
        f.id = f.attributes["ID"][0]
    return f
