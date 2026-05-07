"""eukan compare — compare predicted gene models against a reference."""

from __future__ import annotations

from pathlib import Path

import click
from click_option_group import optgroup

from eukan.cli._framework import PreformattedEpilogCommand


def _derive_compare_labels(paths: tuple[Path, ...]) -> list[str]:
    """Default labels for --predicted: file stems, dedup'd by numeric suffix."""
    from collections import Counter

    stems = [p.stem for p in paths]
    counts = Counter(stems)
    if all(c == 1 for c in counts.values()):
        return stems

    seen: Counter[str] = Counter()
    out: list[str] = []
    for stem in stems:
        if counts[stem] > 1:
            seen[stem] += 1
            out.append(f"{stem}_{seen[stem]}")
        else:
            out.append(stem)
    return out


@click.command(cls=PreformattedEpilogCommand)
@optgroup.group("Required input")
@optgroup.option(
    "--reference", "-r", required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Reference GFF3 file.",
)
@optgroup.option(
    "--predicted", "-p", required=True, multiple=True,
    type=click.Path(exists=True, path_type=Path),
    help="Predicted GFF3 file. Repeat to compare multiple predictions "
         "against the same reference.",
)
@optgroup.group("Pipeline parameters")
@optgroup.option(
    "--label", "-L", multiple=True,
    help="Short label per --predicted (must appear once per --predicted, "
         "or be omitted to use file stems).",
)
@optgroup.group("Output options")
@optgroup.option(
    "--output-file", "-o", type=click.Path(path_type=Path), default=None,
    help="Write per-feature details to a TSV file. In multi-prediction mode "
         "a leading 'prediction' column is prepended.",
)
def compare(
    reference: Path,
    predicted: tuple[Path, ...],
    label: tuple[str, ...],
    output_file: Path | None,
) -> None:
    """Compare predicted gene models against a reference GFF3. The reference
    can be manually curated models, or predictions from another tool or
    pipeline.

    \b
    Computes quality metrics at gene, mRNA, CDS, and intron levels. Repeat
    --predicted to evaluate multiple predictions against the same reference;
    this additionally tabulates, for each gene-level classification
    (match/missing/merged/fragmented), the powerset of predictions that
    agreed on that classification per reference gene.

    \b
    Use --output-file to write a per-feature TSV for further analysis.
    """
    from eukan.compare import (
        compare_annotations,
        compare_multiple,
        format_multi_results,
        format_results,
        write_details_tsv,
    )

    n = len(predicted)
    if label and len(label) != n:
        raise click.UsageError(
            f"--label specified {len(label)} times but --predicted "
            f"specified {n} times; counts must match"
        )

    # Single-prediction path: byte-identical to the legacy behaviour.
    if n == 1:
        if label:
            click.echo(
                "Note: --label is only used in multi-prediction mode; ignored.",
                err=True,
            )
        pred = predicted[0]
        single = compare_annotations(reference.resolve(), pred.resolve())
        click.echo(format_results(single))
        if output_file is not None:
            tsv_path = write_details_tsv(single, output_file.resolve())
            click.echo(f"\nDetails written to {tsv_path}")
        return

    labels = list(label) if label else _derive_compare_labels(predicted)
    if len(set(labels)) != len(labels):
        # Only reachable when user supplied colliding --label values.
        raise click.UsageError("--label values must be unique")
    if not label:
        # Auto-derived labels could collide on stems; warn if dedup happened.
        original_stems = [p.stem for p in predicted]
        if len(set(original_stems)) != len(original_stems):
            click.echo(
                "Note: --predicted files share stems; labels auto-numbered. "
                "Use --label to override.",
                err=True,
            )

    multi = compare_multiple(
        reference.resolve(),
        [(la, p.resolve()) for la, p in zip(labels, predicted, strict=True)],
    )
    click.echo(format_multi_results(multi))
    if output_file is not None:
        tsv_path = write_details_tsv(multi, output_file.resolve())
        click.echo(f"\nDetails written to {tsv_path}")
