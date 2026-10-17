"""eukan func-annot — UniProt + Pfam functional annotation."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import click
from click_option_group import optgroup

from eukan.cli._framework import (
    PreformattedEpilogCommand,
    drop_none,
    force_option,
    numcpu_option,
)
from eukan.infra.databases import HOMOLOGY_DBS
from eukan.infra.layout import step_work_dir


@click.command("func-annot", cls=PreformattedEpilogCommand)
@optgroup.group("Pipeline parameters")
@numcpu_option
@optgroup.option(
    "--homology-db", type=click.Choice(HOMOLOGY_DBS, case_sensitive=False),
    default="uniprot", show_default=True,
    help="Homology source: 'uniprot' runs phmmer vs SwissProt (broad coverage); "
         "'kofam' runs hmmscan vs the KOfam HMM database with per-KO bit-score "
         "thresholds (KEGG-pathway focused). Pfam hmmscan runs in both modes.",
)
@optgroup.option("--evalue", "-e", type=str, default="1e-1", show_default=True, help="E-value cutoff.")
@optgroup.group("Override options")
@optgroup.option(
    "--proteins", "-p", type=click.Path(exists=True, resolve_path=True, path_type=Path),
    help="Amino acid sequences in FASTA format.",
)
@optgroup.option(
    "--uniprot", "uniprot_db", type=click.Path(exists=True, resolve_path=True, path_type=Path),
    default=None, help="UniProt-SwissProt database FASTA.",
)
@optgroup.option(
    "--kofam", "kofam_db", type=click.Path(exists=True, resolve_path=True, path_type=Path),
    default=None, help="KOfam pressed HMM database.",
)
@optgroup.option(
    "--ko-list", "ko_list_path", type=click.Path(exists=True, resolve_path=True, path_type=Path),
    default=None, help="KOfam ko_list TSV (per-KO thresholds + definitions).",
)
@optgroup.option(
    "--pfam", "pfam_db", type=click.Path(exists=True, resolve_path=True, path_type=Path),
    default=None, help="Pfam HMM database.",
)
@optgroup.option(
    "--gff3", "gff3_path", type=click.Path(exists=True, resolve_path=True, path_type=Path),
    default=None, help="GFF3 file to annotate with functional info.",
)
@force_option
def func_annot(**opts: Any) -> None:
    """Add functional annotations (UniProt/KOfam + Pfam) to proteins.

    \b
    The proteins are not discovered automatically: extract the predicted
    proteins of an `eukan annotate` run with `eukan gff3toseq` and pass
    them with --proteins; pass annotate/final.gff3 with --gff3 to get the
    annotated final.mod.gff3. The homology and Pfam databases default to
    where `eukan db-fetch` puts them (databases/ in the current directory).

    \b
    --homology-db uniprot (default) runs phmmer against UniProt-SwissProt
    and emits inference=similar to AA sequence:UniProtKB:... per hit.
    --homology-db kofam runs hmmscan against KEGG's KOfam HMM database
    and emits product=<KO definition>, ec_number=<EC>, Dbxref=KEGG:K...,
    inference=protein motif:KOFAM:K... when score >= the per-KO threshold.
    """
    from eukan.functional import run_functional_annotation  # lazy: startup cost
    from eukan.settings import FunctionalConfig  # lazy: startup cost

    if opts["proteins"] is None:
        raise click.UsageError(
            "No protein file given. Pass --proteins; extract the proteins an annotate run "
            "predicted with `eukan gff3toseq -g <genome> -i annotate/final.gff3 "
            "--output-format protein -o proteins.faa`."
        )

    force = opts.pop("force")
    config = FunctionalConfig(**drop_none(
        **opts, work_dir=step_work_dir("func-annot"), manifest_dir=Path.cwd(),
    ))
    run_functional_annotation(config, force=force)
    click.echo("Done.")
