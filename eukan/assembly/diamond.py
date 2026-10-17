"""``diamond blastx`` of the Trinity transcript models against the ``--uniprot`` DB.

strand_correct and defuse each run one blastx over the same query set against the
same DB. strand_correct asks for the best hit per transcript at default
sensitivity. defuse asks for every HSP with ``--ultra-sensitive``. The DB is
compiled from a FASTA once and cached. Each query is named ``{track}:{tid}`` and
written 5'->3' in the model's labelled orientation, which is what
:func:`coding_strand` reads a hit's frame against.
"""

from __future__ import annotations

import os
from pathlib import Path

from eukan.gff.transcripts import Transcript, stitch
from eukan.infra.genome import ContigIndex
from eukan.infra.logging import get_logger
from eukan.infra.runner import run_cmd
from eukan.settings import AssemblyConfig

log = get_logger(__name__)

# diamond memory caps for the 15 GB / 0-swap box (block-size in billions of query
# letters; index-chunks splits the reference load).
_BLOCK_SIZE = "1.0"
_INDEX_CHUNKS = "4"


def coding_strand(label: str, qframe: int) -> str:
    """The coding strand implied by a blastx hit's query frame on a labelled query.

    The query was stitched in its labelled orientation (``.`` forward, like ``+``),
    so a positive ``qframe`` means that orientation is coding and a negative one
    means its reverse complement is. Returns ``"+"`` or ``"-"``.
    """
    base = label if label in ("+", "-") else "+"
    if qframe > 0:
        return base
    return "-" if base == "+" else "+"


def _resolve_db(config: AssemblyConfig) -> str:
    """Return the diamond ``--db`` base (no ``.dmnd``); build it from a FASTA if needed.

    A ``.dmnd`` path is used directly. A FASTA is compiled with ``diamond makedb``,
    cached as ``<stem>.dmnd`` beside the FASTA when that dir is writable, else in the
    work dir, and reused while it is at least as new as the source.
    """
    src = config.uniprot_db
    assert src is not None  # gated by the caller
    if src.suffix == ".dmnd":
        return str(src.with_suffix(""))  # diamond --db wants the base name

    base = src.with_suffix("")
    out = Path(f"{base}.dmnd")
    if not (src.parent.exists() and os.access(src.parent, os.W_OK)):
        base = config.work_dir / "uniprot_sprot"
        out = Path(f"{base}.dmnd")
    if not (out.exists() and out.stat().st_mtime >= src.stat().st_mtime):
        log.info("Building diamond DB %s from %s ...", out.name, src.name)
        run_cmd(
            ["diamond", "makedb", "--in", str(src), "--db", str(base),
             "--threads", str(config.num_cpu), "--quiet"],
            cwd=config.work_dir,
        )
    return str(base)


def run_blastx(
    config: AssemblyConfig,
    models_by_tag: dict[str, list[Transcript]],
    contigs: ContigIndex,
    query: Path,
    out: Path,
    *,
    evalue: float,
    max_target_seqs: int,
    fields: tuple[str, ...],
    ultra_sensitive: bool = False,
) -> None:
    """Write every model to *query* and blastx it, both strands, into *out*.

    Each set's models go out sorted by chrom and start, as ``>{tag}:{tid}`` and the
    spliced sequence; *out* holds diamond's ``--outfmt 6`` *fields*.
    """
    db = _resolve_db(config)
    with open(query, "w") as fh:
        for tag, models in models_by_tag.items():
            for tx in sorted(models, key=lambda t: (t.chrom, t.exons[0][0])):
                fh.write(f">{tag}:{tx.tid}\n{stitch(tx, contigs)}\n")
    run_cmd(
        ["diamond", "blastx",
         "--db", db,
         "--query", str(query),
         "--out", str(out),
         *(["--ultra-sensitive"] if ultra_sensitive else []),
         "--strand", "both",
         "--query-gencode", str(config.genetic_code_obj.ncbi_id),
         "--evalue", f"{evalue:g}",
         "--max-target-seqs", str(max_target_seqs),
         "--outfmt", "6", *fields,
         "--block-size", _BLOCK_SIZE,
         "--index-chunks", _INDEX_CHUNKS,
         "--threads", str(config.num_cpu),
         "--quiet"],
        cwd=config.work_dir,
    )
