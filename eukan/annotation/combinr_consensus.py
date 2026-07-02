"""Consensus gene models via the external ``combinr consensus`` engine.

The consensus engine. ``combinr consensus`` integrates weighted evidence — ab
initio predictions, protein alignments, and transcript alignments — into one
best-scoring coding model per locus, and folds UTRs and alternative isoforms in
from the transcript evidence (via ``--alt-splice``) in a single pass.

The evidence is staged into the ``consensus_models`` step dir and written to
``consensus_models.gff3``, so the shared tail in
:func:`eukan.annotation.consensus.build_consensus_models` (ORF patch +
prettification) consumes it directly.

Input contract (verified against combinr's ``src/consensus/evidence.rs``):

* ``--gene-predictions`` reads ``CDS`` rows grouped by ``Parent`` — eukan's
  concatenated ab initio ``gene_predictions.gff3`` is consumed as-is.
* ``--protein-alignments`` / ``--transcript-alignments`` read spliced **match**
  chains carrying ``Target=``; rows without ``Target=`` are silently skipped. An
  evidence chain's *class* (PROTEIN / TRANSCRIPT / ABINITIO_PREDICTION) comes
  from the weights file keyed on the GFF source token, not from the flag.
* eukan's ``prot.gff3`` (spaln/gth) is CDS-format and its ``nr_transcripts.gff3``
  is flat ``exon``/``Parent`` — neither carries ``Target=`` — so both are
  converted here into match chains via :func:`_chains_to_match`.
"""

from __future__ import annotations

from pathlib import Path

from eukan.annotation.evidence import EVIDENCE_ROLES, _first_source_token
from eukan.infra.logging import get_logger
from eukan.infra.runner import run_cmd
from eukan.settings import PipelineConfig

log = get_logger(__name__)


def _combinr_bin(config: PipelineConfig) -> str:
    """The combinr executable: explicit ``combinr_path`` or ``combinr`` on PATH."""
    return str(config.combinr_path) if config.combinr_path else "combinr"


def _parse_attrs(col9: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for part in col9.split(";"):
        part = part.strip()
        if "=" in part:
            key, val = part.split("=", 1)
            attrs[key.strip()] = val.strip()
    return attrs


def _chains_to_match(in_gff: Path, out_gff: Path, *, feature_type: str, match_type: str) -> int:
    """Rewrite ``feature_type`` rows grouped by ``Parent`` into ``Target=`` match chains.

    combinr's consensus alignment parser groups rows by their ``ID`` and requires a
    ``Target=`` attribute, so this turns eukan's CDS-based protein alignments
    (``feature_type="CDS"``) and flat-exon transcript alignments
    (``feature_type="exon"``) into the spliced-match form combinr expects. The
    source column (col 2) is preserved verbatim so it still matches the weights
    entry. Target coordinates are cumulative over the sorted blocks (combinr only
    requires ``Target`` to be present, not meaningful). Returns the chain count.
    """
    chains: dict[str, list] = {}  # key -> [chrom, source, strand, [(lend, rend)]]
    order: list[str] = []
    with open(in_gff) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != feature_type:
                continue
            attrs = _parse_attrs(cols[8])
            key = attrs.get("Parent") or attrs.get("ID")
            if not key:
                continue
            rec = chains.get(key)
            if rec is None:
                rec = [cols[0], cols[1], cols[6], []]
                chains[key] = rec
                order.append(key)
            rec[3].append((int(cols[3]), int(cols[4])))

    with open(out_gff, "w") as out:
        for key in order:
            chrom, source, strand, blocks = chains[key]
            blocks.sort()
            tpos = 1
            for lend, rend in blocks:
                tend = tpos + (rend - lend)
                out.write(
                    f"{chrom}\t{source}\t{match_type}\t{lend}\t{rend}\t.\t{strand}\t."
                    f"\tID={key};Target={key} {tpos} {tend}\n"
                )
                tpos = tend + 1
    return len(order)


def _stage_combinr_inputs(
    config: PipelineConfig,
    sdir: Path,
    evidence: list[Path],
    transcripts: Path | None,
) -> bool:
    """Stage gene_predictions.gff3, prot.match.gff3, transcripts.match.gff3, weights.txt.

    Uses the same ab initio concatenation and the same weight tokens (via the
    shared :data:`~eukan.annotation.evidence.EVIDENCE_ROLES`) the consensus engine
    scores. The protein and transcript files are converted to
    ``Target=`` match chains rather than symlinked. Returns ``True`` when
    transcript evidence was staged.
    """
    weights = [str(w) for w in config.weights]
    weight_lines: list[str] = []

    with open(sdir / "gene_predictions.gff3", "wb") as pf:
        for ev in evidence:
            if ev.name == "prot.gff3":
                _chains_to_match(
                    ev, sdir / "prot.match.gff3",
                    feature_type="CDS", match_type="nucleotide_to_protein_match",
                )
                weight_lines.append("\t".join(["PROTEIN", "prot_align", weights[0]]))
                continue
            role = EVIDENCE_ROLES.get(ev.name)
            if role is not None:
                _cls, token = role
                weight_lines.append("\t".join(["ABINITIO_PREDICTION", token, weights[1]]))
            with open(ev, "rb") as ef:
                pf.write(ef.read())

    have_transcripts = bool(config.has_transcripts and transcripts is not None)
    if have_transcripts:
        assert transcripts is not None
        source = _first_source_token(transcripts) or "transcript"
        _chains_to_match(
            transcripts, sdir / "transcripts.match.gff3",
            feature_type="exon", match_type="cDNA_match",
        )
        weight_lines.append("\t".join(["TRANSCRIPT", source, weights[2]]))

    (sdir / "weights.txt").write_text("\n".join(weight_lines) + "\n")
    return have_transcripts


def run_combinr_consensus(
    config: PipelineConfig,
    sdir: Path,
    evidence: list[Path],
    *,
    transcripts: Path | None = None,
) -> Path:
    """Build consensus gene models with ``combinr consensus`` into ``consensus_models.gff3``.

    Transcript evidence (when present) enables ``--alt-splice`` so combinr emits
    alternative isoforms with UTRs derived from the consensus CDS.
    """
    log.info("Running combinr consensus model building...")
    have_transcripts = _stage_combinr_inputs(config, sdir, evidence, transcripts)

    cmd = [
        _combinr_bin(config), "consensus",
        "--weights", "weights.txt",
        "--genome", str(config.genome),
        "--gene-predictions", "gene_predictions.gff3",
        "--genetic-code", str(config.genetic_code),
        "-t", str(config.num_cpu),
        "--format", "gff3",
    ]
    if (sdir / "prot.match.gff3").exists():
        cmd += ["--protein-alignments", "prot.match.gff3"]
    if have_transcripts:
        cmd += [
            "--transcript-alignments", "transcripts.match.gff3",
            "--alt-splice", "--events", "combinr.alt_splice_events.tsv",
            # --stringent-overlap tunes the isoform grouping; 0 = off (any overlap).
            "--stringent-overlap", str(config.combinr_stringent_overlap),
        ]

    run_cmd(cmd, cwd=sdir, out_file="consensus_models.gff3")
    return sdir / "consensus_models.gff3"
