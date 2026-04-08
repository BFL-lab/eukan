"""Transcriptome assembly pipeline: read mapping, Trinity assembly, and PASA alignment."""

from eukan.assembly.orchestrator import run_assembly

__all__ = ["run_assembly"]
