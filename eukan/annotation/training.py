"""Evidence-driven training-set construction for gene predictors.

Used by both AUGUSTUS and SNAP — extracts gene models that two or more
prediction sources agree on, then materialises them in a format the
downstream trainers can read (GenBank flat file, randomly partitioned
into training and test sets).
"""

from __future__ import annotations

from pathlib import Path

from eukan.gff import concordance
from eukan.gff.rows import iter_gff3_rows
from eukan.infra.runner import run_cmd
from eukan.settings import PipelineConfig


def build_training_set(
    config: PipelineConfig, evidence: tuple[Path, ...], sdir: Path
) -> None:
    """Build evidence-driven training set for gene predictors.

    Writes ``training_set.gff3``, ``genbank.gb``, ``genbank.gb.train``,
    and ``genbank.gb.test`` into *sdir*. The training-set size targets
    a quarter of the concordant gene count, with a flank derived from
    the mean gene length of the first evidence source.
    """
    training_db = concordance.extract_supported_models(*evidence, output_dir=sdir)

    num_training = max(
        1, len([f.id for f in training_db.features_of_type("gene")]) // 4
    )

    gene_lengths = [
        int(cols[4]) - int(cols[3]) for cols in iter_gff3_rows(evidence[0]) if cols[2] == "gene"
    ]
    flank = int(sum(gene_lengths) / len(gene_lengths) / 2)

    run_cmd(
        ["gff2gbSmallDNA.pl", "training_set.gff3", str(config.genome), str(flank), "genbank.gb"],
        cwd=sdir,
    )
    run_cmd(["randomSplit.pl", "genbank.gb", str(num_training)], cwd=sdir)
