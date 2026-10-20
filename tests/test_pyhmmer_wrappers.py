"""Tests for the pyhmmer search wrappers in eukan.functional (real pyhmmer runs)."""

import random

import pyhmmer
import pytest

from eukan.functional.kofam import run_kofam_search
from eukan.functional.search import run_hmmscan_search, run_phmmer_search

AMINO = pyhmmer.easel.Alphabet.amino()

N_TARGETS = 20


def _protein(name: str, sequence: str) -> pyhmmer.easel.DigitalSequence:
    return pyhmmer.easel.TextSequence(name=name, sequence=sequence).digitize(AMINO)


@pytest.fixture(scope="module")
def related_first() -> tuple[str, list[pyhmmer.easel.DigitalSequence]]:
    """A query protein and N_TARGETS targets named like K-numbers: K00001 is
    a mutant of the query, the rest are unrelated random sequences."""
    rnd = random.Random(1)
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"

    def random_seq(n: int) -> str:
        return "".join(rnd.choice(amino_acids) for _ in range(n))

    query = random_seq(150)
    related = "".join(rnd.choice(amino_acids) if rnd.random() < 0.5 else aa for aa in query)
    seqs = [related] + [random_seq(150) for _ in range(N_TARGETS - 1)]
    return query, [_protein(f"K{i:05d}", s) for i, s in enumerate(seqs, start=1)]


def _straddling_cutoff(top_hits: pyhmmer.plan7.TopHits) -> float:
    """An E-value cutoff that K00001's hit passes in-scan, where pyhmmer
    gates it against the targets seen so far (it is first: Z=1), while its
    E-value against all N_TARGETS targets is twice the cutoff."""
    related = next(hit for hit in top_hits if hit.name == "K00001")
    return related.pvalue * N_TARGETS / 2


class TestEvalueCutoff:
    """No hit above the requested E-value survives. At 4x the straddling
    cutoff the related hit (E = 2x) is within it and must be kept."""

    def test_phmmer(self, related_first):
        query, targets = related_first
        q = _protein("q", query)
        cutoff = _straddling_cutoff(next(pyhmmer.hmmer.phmmer([q], targets, cpus=1)))

        assert run_phmmer_search([q], targets, 1, cutoff) == {}
        assert "K00001" in run_phmmer_search([q], targets, 1, 4 * cutoff)["q"]

    def test_hmmscan(self, related_first):
        query, targets = related_first
        q = _protein("q", query)
        builder, background = pyhmmer.plan7.Builder(AMINO), pyhmmer.plan7.Background(AMINO)
        hmms = [builder.build(t, background)[0] for t in targets]
        cutoff = _straddling_cutoff(next(pyhmmer.hmmer.hmmscan([q], hmms, cpus=1)))

        assert run_hmmscan_search([q], hmms, 1, cutoff) == {}
        assert "K00001" in run_hmmscan_search([q], hmms, 1, 4 * cutoff)["q"]

    @pytest.mark.parametrize("pressed", [False, True], ids=["unpressed", "pressed"])
    def test_kofam(self, related_first, tmp_path, pressed):
        # A pressed database is read straight into a block of profiles, an unpressed one
        # is configured into one by hmmscan.
        query, targets = related_first
        builder, background = pyhmmer.plan7.Builder(AMINO), pyhmmer.plan7.Background(AMINO)
        db = tmp_path / "kofam.hmm"
        with open(db, "wb") as fh:
            for t in targets:
                builder.build(t, background)[0].write(fh)
        if pressed:
            with pyhmmer.plan7.HMMFile(str(db)) as hf:
                pyhmmer.hmmer.hmmpress(hf, str(db))
        ko_list = tmp_path / "ko_list.tsv"
        ko_list.write_text("knum\tthreshold\tscore_type" + "\tx" * 9 + "\n" + "".join(
            f"{t.name}\t-\tfull" + "\t-" * 8 + "\tsynthetic KO\n" for t in targets
        ))
        proteins = tmp_path / "q.faa"
        proteins.write_text(f">q\n{query}\n")
        with pyhmmer.plan7.HMMFile(str(db)) as hf:
            cutoff = _straddling_cutoff(next(pyhmmer.hmmer.hmmscan([_protein("q", query)], hf, cpus=1)))

        assert run_kofam_search(proteins, db, ko_list, 1, cutoff) == {}
        assert "K00001" in run_kofam_search(proteins, db, ko_list, 1, 4 * cutoff)["q"]
