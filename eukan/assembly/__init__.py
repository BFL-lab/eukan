"""Transcriptome assembly pipeline: read mapping, Trinity assembly (de novo + genome-guided), SL trans-splice cut, combinr consolidation."""

from eukan.assembly.pipeline import run_assembly

__all__ = ["run_assembly"]
