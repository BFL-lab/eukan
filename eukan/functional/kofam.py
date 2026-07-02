"""KOfam (KEGG Orthology HMM) search and KO-list parsing.

Adapts the kofam-scan algorithm to pyhmmer:

* ``parse_ko_list`` loads the per-KO threshold + score_type + definition.
* ``extract_ec_numbers`` peels ``[EC:...]`` tags out of KO definitions so
  they can be emitted as their own ``ec_number=`` GFF3 attribute.
* ``run_kofam_search`` runs ``hmmscan`` (proteome vs the pressed KOfam
  HMM database, like Pfam) and applies the per-KO bit-score threshold
  to flag confident assignments. Each hit's score/E-value is read from
  ``hit.score`` for full-score KOs or ``hit.best_domain.score`` for
  domain-score KOs (matching kofam-scan's behaviour).
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import pyhmmer

from eukan.functional.search import HitResults, _decode, _load_digital_sequences
from eukan.infra.logging import get_logger

log = get_logger(__name__)


# ---------------------------------------------------------------------------
# ko_list parsing
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class KOEntry:
    """A single KO's metadata loaded from ``ko_list``.

    ``threshold`` is None when the KO has no curated cutoff (the file
    stores ``-`` in that column). Hits against such KOs can still appear
    in results, but ``above_threshold`` is always ``False`` for them.
    """

    k_number: str
    threshold: float | None
    score_type: str  # "full" or "domain"
    definition: str


def parse_ko_list(ko_list_path: Path) -> dict[str, KOEntry]:
    """Parse the KOfam ``ko_list`` TSV into a dict keyed by K-number.

    The file has a header row and 12 columns; we only keep K-number,
    threshold, score_type, and definition. Lines with fewer columns
    (truncated, malformed) are skipped silently — the file is large and
    occasionally has stray rows.
    """
    entries: dict[str, KOEntry] = {}
    with open(ko_list_path) as f:
        f.readline()  # header
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 12:
                continue
            k_number, thr, score_type, _profile_type = cols[0], cols[1], cols[2], cols[3]
            definition = cols[11]
            threshold = float(thr) if thr not in ("-", "") else None
            entries[k_number] = KOEntry(
                k_number=k_number,
                threshold=threshold,
                score_type=score_type,
                definition=definition,
            )
    return entries


# ---------------------------------------------------------------------------
# EC number extraction
# ---------------------------------------------------------------------------


# Matches one bracketed EC tag, e.g. " [EC:3.2.1.92]" or "[EC:1.1.1.- 1.1.1.184]".
# Contents can be comma- or whitespace-separated codes (KEGG mixes both).
_EC_BRACKET_RE = re.compile(r"\s*\[EC:([0-9.\-,nN\s]+)\]")

# Splits the captured content on commas and whitespace.
_EC_TOKEN_RE = re.compile(r"[\s,]+")

# Valid EC qualifier values for NCBI: top class 1-7, four levels, each
# either digits, a dash (-), or a preliminary code (n<digits>). The KEGG
# placeholder "n.n.n.n" (each level just "n" with no digits) is filtered
# out because it isn't a real EC number.
_VALID_EC_RE = re.compile(r"^[1-7](?:\.(?:\d+|-|n\d+)){3}$")


def extract_ec_numbers(definition: str) -> tuple[str, list[str]]:
    """Strip ``[EC:...]`` tags from a KO definition; return cleaned text + EC list.

    Drops placeholder codes like ``1.n.n.n`` (KEGG's "no EC assigned"
    marker) — they are valid in the KEGG TSV but not as NCBI
    ``ec_number=`` qualifier values.
    """
    ec_numbers: list[str] = []

    def collect(m: re.Match[str]) -> str:
        for tok in _EC_TOKEN_RE.split(m.group(1).strip()):
            if tok and _VALID_EC_RE.match(tok):
                ec_numbers.append(tok)
        return ""

    cleaned = _EC_BRACKET_RE.sub(collect, definition).strip()
    return cleaned, list(dict.fromkeys(ec_numbers))


# ---------------------------------------------------------------------------
# Search via pyhmmer hmmscan (mirrors the Pfam search path)
# ---------------------------------------------------------------------------


def _hit_score_and_evalue(
    hit: pyhmmer.plan7.Hit, score_type: str,
) -> tuple[float, float]:
    """Pick full-sequence vs best-domain score per the KO's ``score_type``.

    Falls back to full-sequence numbers when ``score_type == "domain"``
    but the hit has no domain data (shouldn't happen in practice, but
    keeps the loader robust against odd HMM outputs).
    """
    if score_type == "domain" and hit.domains:
        best = hit.best_domain
        return float(best.score), float(best.i_evalue)
    return float(hit.score), float(hit.evalue)


def run_kofam_search(
    proteins: Path,
    kofam_db: Path,
    ko_list_path: Path,
    num_cpu: int,
    evalue_threshold: float,
) -> HitResults:
    """Run hmmscan of the proteome against the pressed KOfam HMM database.

    Returns ``query_id → {K-number → HitInfo}``. Hits below the per-KO
    bit-score threshold are kept but marked ``above_threshold=False`` so
    the GFF3 emitter can choose how to render them; the global
    ``evalue_threshold`` is applied as a coarse pre-filter (mirrors
    kofam-scan's ``-E``).
    """
    log.info("Loading KO list from %s", ko_list_path)
    ko_entries = parse_ko_list(ko_list_path)
    log.info("Loaded %d KO entries (with %d curated thresholds)",
             len(ko_entries),
             sum(1 for e in ko_entries.values() if e.threshold is not None))

    log.info("Loading proteome from %s", proteins)
    queries = _load_digital_sequences(proteins)

    log.info(
        "Running KOfam hmmscan (%d queries vs %d profiles, %d CPUs, E≤%s)...",
        len(queries),
        len(ko_entries),
        num_cpu,
        evalue_threshold,
    )

    results: HitResults = {}
    with pyhmmer.plan7.HMMFile(str(kofam_db)) as hf:
        for top_hits in pyhmmer.hmmer.hmmscan(
            queries, hf, cpus=num_cpu, E=evalue_threshold,
        ):
            query_name = _decode(top_hits.query.name)
            if not top_hits:
                continue

            for hit in top_hits:
                k_number = _decode(hit.name)
                ko = ko_entries.get(k_number)
                if ko is None:
                    # HMM not in ko_list (database mismatch or trimmed
                    # ko_list). Skip — without metadata we can't score it.
                    continue

                score, evalue = _hit_score_and_evalue(hit, ko.score_type)
                above = ko.threshold is not None and score >= ko.threshold

                cleaned_def, ec_numbers = extract_ec_numbers(ko.definition)

                per_query = results.setdefault(query_name, {})
                prev = per_query.get(k_number)
                if prev is not None and prev["score"] >= score:
                    continue
                per_query[k_number] = {
                    "description": cleaned_def,
                    "evalue": evalue,
                    "score": score,
                    "threshold": ko.threshold if ko.threshold is not None else 0.0,
                    "score_type": ko.score_type,
                    "above_threshold": above,
                    "ec_numbers": ec_numbers,
                }

    n_above = sum(
        1 for hits in results.values()
        if any(h.get("above_threshold") for h in hits.values())
    )
    log.info(
        "KOfam: %d queries with hits (%d above per-KO threshold)",
        len(results), n_above,
    )
    return results
