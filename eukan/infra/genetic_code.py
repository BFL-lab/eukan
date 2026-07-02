"""Centralised genetic code handling.

Provides a single :class:`GeneticCode` type that wraps an NCBI genetic code ID
and exposes codon-table lookups, tool-specific CLI flag generation, and
stop-codon lists.
"""

from __future__ import annotations

from Bio.Data.CodonTable import unambiguous_dna_by_id

from eukan.exceptions import InvalidOptionError

# Genetic codes supported by GeneMark-ES/ET (gmes_petap.pl --gcode).
# Code 26 is remapped to 6 internally by GeneMark.
_GENEMARK_CODES: set[int] = {1, 6, 26}


class GeneticCode:
    """Wrapper around an NCBI genetic code ID.

    Accepts ``int`` or ``str`` (the string is cast to ``int``).  Raises
    :class:`~eukan.exceptions.InvalidOptionError` if the code is not a
    recognised NCBI translation table.

    Examples::

        gc = GeneticCode(6)
        gc.genemark_flag    # ["--gcode=6"]
        gc.stop_codons      # ["TGA"]
        gc.codon_table      # Bio.Data.CodonTable object
    """

    __slots__ = ("_ncbi_id", "_table")

    def __init__(self, code: int | str) -> None:
        try:
            ncbi_id = int(code)
        except (ValueError, TypeError) as exc:
            raise InvalidOptionError(
                f"Invalid genetic code: {code!r}",
                hint="Provide a numeric NCBI translation table ID (e.g. 1, 6, 10, 12).",
            ) from exc

        if ncbi_id not in unambiguous_dna_by_id:
            raise InvalidOptionError(
                f"Unknown NCBI genetic code: {ncbi_id}",
                hint=f"Valid codes: {', '.join(str(k) for k in sorted(unambiguous_dna_by_id))}",
            )

        self._ncbi_id = ncbi_id
        self._table = unambiguous_dna_by_id[ncbi_id]

    # -- properties ----------------------------------------------------------

    @property
    def ncbi_id(self) -> int:
        return self._ncbi_id

    @property
    def codon_table(self):
        """The BioPython ``CodonTable`` for this code."""
        return self._table

    @property
    def stop_codons(self) -> list[str]:
        return self._table.stop_codons

    @property
    def genemark_flag(self) -> list[str]:
        """CLI args for GeneMark's ``--gcode`` option.

        Returns an empty list for code 1 (default) or unsupported codes.
        """
        if self._ncbi_id in _GENEMARK_CODES and self._ncbi_id != 1:
            return [f"--gcode={self._ncbi_id}"]
        return []

    @property
    def is_genemark_supported(self) -> bool:
        return self._ncbi_id in _GENEMARK_CODES

    # -- dunder --------------------------------------------------------------

    def __repr__(self) -> str:
        return f"GeneticCode({self._ncbi_id})"

    def __eq__(self, other: object) -> bool:
        if isinstance(other, GeneticCode):
            return self._ncbi_id == other._ncbi_id
        return NotImplemented

    def __hash__(self) -> int:
        return hash(self._ncbi_id)

    def __int__(self) -> int:
        return self._ncbi_id

    def __str__(self) -> str:
        return str(self._ncbi_id)
