"""SNAP and CodingQuarry gene prediction."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import gffutils
from Bio import SeqIO

from eukan.annotation.training import build_training_set
from eukan.gff import GFF3_DIALECT, create_gff_db, transform_db
from eukan.gff import transforms as gff_transforms
from eukan.gff.normalize import normalize_to_gff3
from eukan.infra.logging import get_logger
from eukan.infra.runner import run_cmd
from eukan.infra.steps import step_dir
from eukan.infra.utils import symlink
from eukan.settings import PipelineConfig

log = get_logger(__name__)


def _read_snap_output(snap_gff: Path) -> gffutils.FeatureDB:
    """Parse SNAP's faux-GTF output into a normalised GFF3 FeatureDB.

    SNAP emits a flat exon list keyed by gene group; we first relabel
    every line as ``exon`` with synthetic IDs (the per-call counter in
    ``_make_snap_featuretype_transform``), then re-parse as GTF so
    ``Parent=`` is treated as a transcript_id and the gene/mRNA hierarchy
    is materialised.
    """
    snap = create_gff_db(
        snap_gff, transform=_make_snap_featuretype_transform(),
    )
    dialect = gffutils.DataIterator(str(snap_gff)).dialect
    dialect["fmt"] = "gtf"
    snap = gffutils.create_db(
        snap, ":memory:", dialect=dialect,
        gtf_transcript_key="Parent", gtf_gene_key="Parent",
    )
    snap.dialect = GFF3_DIALECT
    return transform_db(snap, gff_transforms.gff3_it)


def run_snap(config: PipelineConfig, *evidence: Path) -> Path:
    """Train and run SNAP gene prediction."""
    output = "snap.gff3"
    sdir = step_dir(config.work_dir, "snap")
    log.info("Running SNAP training and prediction...")
    symlink(config.genome, sdir / "genome.dna")

    build_training_set(config, evidence, sdir)
    _gff_to_zff(sdir / "training_set.gff3", str(config.genome), output_dir=sdir)

    run_cmd(["fathom", "-categorize", "1000", "genome.ann", "genome.dna"], cwd=sdir)
    run_cmd(["fathom", "-export", "1000", "-plus", "uni.ann", "uni.dna"], cwd=sdir)
    run_cmd(["forge", "export.ann", "export.dna"], cwd=sdir)
    run_cmd(["hmm-assembler.pl", "snap", "."], cwd=sdir, out_file="snap.hmm")
    run_cmd(
        ["snap", "-gff", "snap.hmm", str(config.genome)],
        cwd=sdir, out_file="snap.gff",
    )

    return normalize_to_gff3(
        _read_snap_output(sdir / "snap.gff"),
        sdir / output,
        post_transform=_snap_homogenize_source,
    )


def run_codingquarry(config: PipelineConfig, evidence: Path) -> Path:
    """Run CodingQuarry gene prediction (fungi/protists)."""
    output = "codingquarry.gff3"
    sdir = step_dir(config.work_dir, "codingquarry")
    log.info("Running CodingQuarry gene prediction...")

    flag = "-t" if config.transcripts_gff else "-a"
    run_cmd(
        [
            "CodingQuarry",
            "-p", str(config.num_cpu),
            "-f", str(config.genome),
            flag, str(evidence),
        ],
        cwd=sdir,
    )

    return normalize_to_gff3(
        sdir / "out" / "PredictedPass.gff3", sdir / output,
        parse_transform=_codingquarry_add_mRNA,
        post_transform=_codingquarry_fix_dup_ids,
    )


# ---------------------------------------------------------------------------
# SNAP output transforms
# ---------------------------------------------------------------------------


def _make_snap_featuretype_transform() -> Callable[[gffutils.Feature], gffutils.Feature]:
    """Return a fresh transform that normalizes SNAP output features.

    Each call creates an independent counter so that multiple transform
    passes don't share state (which would cause ID collisions).
    """
    counter: dict[str, int] = {}

    def fix_snap_featuretype(f: gffutils.Feature) -> gffutils.Feature:
        f.source = "snap"
        f.featuretype = "exon"
        # SNAP puts the gene group name as the first attribute key
        parent = next(iter(f.attributes))
        counter[parent] = counter.get(parent, 0) + 1
        exon_id = f"{parent}_exon_{counter[parent]}"
        f.attributes = {"ID": [exon_id], "Parent": [parent]}
        return f

    return fix_snap_featuretype


def _snap_homogenize_source(f: gffutils.Feature) -> gffutils.Feature:
    """Set ``source`` to ``snap`` for all features."""
    f.source = "snap"
    return f


# ---------------------------------------------------------------------------
# CodingQuarry output transforms
# ---------------------------------------------------------------------------


def _codingquarry_add_mRNA(f: gffutils.Feature) -> gffutils.Feature:
    """Fix CodingQuarry's incomplete output by attaching mRNA parents."""
    f.source = "codingquarry"
    if f.featuretype == "CDS":
        f.attributes["Parent"][0] += "_mRNA"
    return f


def _codingquarry_fix_dup_ids(f: gffutils.Feature) -> gffutils.Feature:
    """Fix duplicate IDs in CodingQuarry output."""
    if f.featuretype == "mRNA":
        f.id = f.attributes["ID"][0]
    else:
        f.attributes["ID"] = [f.id]
    return f


# ---------------------------------------------------------------------------
# GFF3 → ZFF conversion (SNAP's training input format)
# ---------------------------------------------------------------------------


def _gff_to_zff(
    gff3: str | Path, fasta: str | Path, output_dir: Path | None = None,
) -> None:
    """Convert GFF3 to SNAP's ZFF annotation format (``genome.ann``)."""
    featuredb = create_gff_db(gff3)
    contig_list = [f.id for f in SeqIO.parse(str(fasta), "fasta")]
    contigs: dict[str, list[list]] = {k: [] for k in contig_list}

    for mRNA in featuredb.features_of_type("mRNA"):
        exons = list(featuredb.children(mRNA, featuretype="exon"))
        tot_exons = len(exons)
        parent = mRNA.id

        for i, exon in enumerate(exons, start=1):
            if tot_exons == 1:
                feat = "Esngl"
            elif i == 1:
                feat = "Einit" if mRNA.strand == "+" else "Eterm"
            elif i == tot_exons:
                feat = "Eterm" if mRNA.strand == "+" else "Einit"
            else:
                feat = "Exon"

            start, end = (
                (exon.start, exon.end)
                if mRNA.strand == "+"
                else (exon.end, exon.start)
            )
            contigs[mRNA.chrom].append([feat, start, end, parent])

    out_path = (output_dir / "genome.ann") if output_dir else Path("genome.ann")
    with open(out_path, "w") as fout:
        for contig, entries in contigs.items():
            fout.write(f">{contig}\n")
            for feat, start, end, parent in entries:
                fout.write(f"{feat}  {start} {end} {parent}\n")
