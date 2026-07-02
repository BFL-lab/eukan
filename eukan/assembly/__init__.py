"""Transcriptome assembly pipeline: read mapping, StringTie + rnaSPAdes assembly, SL trans-splice cut, combinr consolidation."""

from eukan.assembly.pipeline import run_assembly

__all__ = ["run_assembly"]
