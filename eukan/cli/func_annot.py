"""eukan func-annot — UniProt + Pfam functional annotation."""

from __future__ import annotations

from pathlib import Path

import click
from click_option_group import optgroup

from eukan.cli._framework import (
    PreformattedEpilogCommand,
    drop_none,
    force_option,
    numcpu_option,
)


@click.command("func-annot", cls=PreformattedEpilogCommand)
@optgroup.group("Pipeline parameters")
@numcpu_option
@optgroup.option(
    "--homology-db", type=click.Choice(["uniprot", "kofam"], case_sensitive=False),
    default="uniprot", show_default=True,
    help="Homology source: 'uniprot' runs phmmer vs SwissProt (broad coverage); "
         "'kofam' runs hmmscan vs the KOfam HMM database with per-KO bit-score "
         "thresholds (KEGG-pathway focused). Pfam hmmscan runs in both modes.",
)
@optgroup.option("--evalue", "-e", type=str, default="1e-1", show_default=True, help="E-value cutoff.")
@optgroup.group("Override options")
@optgroup.option(
    "--proteins", "-p", type=click.Path(exists=True, path_type=Path),
    help="Amino acid sequences in FASTA format.",
)
@optgroup.option(
    "--uniprot", type=click.Path(exists=True, path_type=Path),
    default=None, help="UniProt-SwissProt database FASTA.",
)
@optgroup.option(
    "--kofam", type=click.Path(exists=True, path_type=Path),
    default=None, help="KOfam pressed HMM database.",
)
@optgroup.option(
    "--ko-list", type=click.Path(exists=True, path_type=Path),
    default=None, help="KOfam ko_list TSV (per-KO thresholds + definitions).",
)
@optgroup.option(
    "--pfam", type=click.Path(exists=True, path_type=Path),
    default=None, help="Pfam HMM database.",
)
@optgroup.option(
    "--gff3", type=click.Path(exists=True, path_type=Path),
    default=None, help="GFF3 file to annotate with functional info.",
)
@force_option
def func_annot(
    proteins: Path,
    homology_db: str,
    uniprot: Path | None,
    kofam: Path | None,
    ko_list: Path | None,
    pfam: Path | None,
    gff3: Path | None,
    numcpu: int,
    evalue: str,
    force: bool,
) -> None:
    """Add functional annotations (UniProt/KOfam + Pfam) to proteins.

    \b
    When run after `eukan annotate` and `eukan db-fetch`, the predicted
    protein sequences and homology/Pfam databases are discovered
    automatically. Use the override options to point to different files
    or to run functional annotation independently of the main pipeline.

    \b
    --homology-db uniprot (default) runs phmmer against UniProt-SwissProt
    and emits inference=similar to AA sequence:UniProtKB:... per hit.
    --homology-db kofam runs hmmscan against KEGG's KOfam HMM database
    and emits product=<KO definition>, ec_number=<EC>, Dbxref=KEGG:K...,
    inference=protein motif:KOFAM:K... when score >= the per-KO threshold.
    """
    from eukan.functional import run_functional_annotation
    from eukan.infra.layout import step_work_dir
    from eukan.settings import FunctionalConfig

    if proteins is None:
        raise click.UsageError(
            "No protein file found. Provide --proteins or run `eukan annotate` first."
        )

    config = FunctionalConfig(**drop_none(
        proteins=proteins.resolve(),
        work_dir=step_work_dir("func-annot"),
        manifest_dir=Path.cwd(),
        num_cpu=numcpu,
        evalue=evalue,
        homology_db=homology_db.lower(),
        uniprot_db=uniprot.resolve() if uniprot else None,
        kofam_db=kofam.resolve() if kofam else None,
        ko_list_path=ko_list.resolve() if ko_list else None,
        pfam_db=pfam.resolve() if pfam else None,
        gff3_path=gff3.resolve() if gff3 else None,
    ))
    run_functional_annotation(config, force=force)
    click.echo("Done.")
