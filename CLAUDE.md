# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Eukan is a eukaryotic genome annotation pipeline that integrates ab initio gene prediction (GeneMark, AUGUSTUS, SNAP, CodingQuarry) with homology-based evidence (protein alignments via spaln/gth) and transcript assemblies to produce consensus gene models via EVidenceModeler (EVM). It optionally adds UTRs via PASA and functional annotation via phmmer/hmmscan.

## Build and Run

Uses Poetry for package management. The pipeline runs inside Docker with a GeneMark-ES/ET/EP+ license required.

```bash
# Install locally for development
poetry install

# Run tests
poetry run pytest tests/ -v

# CLI (all subcommands)
poetry run eukan --help
poetry run eukan annotate -g genome.fasta -p proteins.fasta --kingdom protist
poetry run eukan assemble -g genome.fasta -l left.fq -r right.fq -A -T -P
poetry run eukan func-annot -p proteins.faa --gff3 genes.gff3
poetry run eukan gff3toseq -g genome.fa -i genes.gff3 -o protein
poetry run eukan db-fetch -o databases/

# Dev tooling (not exposed via main CLI)
python tests/run_pipeline.py setup-test-data
python tests/run_pipeline.py test-pipeline --kingdom fungus -n 8
python tests/run_pipeline.py clean-test-data --all
python scripts/generate-env.py -o environment.yml

# Docker build and run
docker build -t eukan -f docker/Dockerfile .
./eukan-docker annotate -g genome.fasta -p proteins.fasta --kingdom protist
```

## Architecture

### CLI (`eukan/cli.py`)

Click-based CLI with subcommands: `annotate`, `assemble`, `func-annot`, `gff3toseq`, `db-fetch`, `check`, `status`. Entry point defined in `pyproject.toml` as `eukan = "eukan.cli:cli"`.

### Package Structure

```
eukan/
├── cli.py              # Click CLI entry points
├── settings.py         # PipelineConfig, AssemblyConfig, FunctionalConfig (pydantic-settings)
├── check.py            # Pre-flight checks for external tools and databases
│
├── infra/              # Runtime infrastructure
│   ├── runner.py       # run_cmd(), run_piped(), run_parallel() — subprocess execution
│   ├── manifest.py     # RunManifest, pipeline_step() — run tracking and reproducibility
│   ├── steps.py        # step_dir(), step_complete() — step directory management
│   └── logging.py      # get_logger(), setup_logging(), md5_file(), validate_gff()
│
├── gff/                # GFF3 format operations
│   ├── parser.py       # Transform callbacks for gffutils.create_db(transform=fn)
│   ├── intersecter.py  # Genomic interval operations (concordance, overlap, merging)
│   └── io.py           # featuredb2gff3_file(), extract_sequences()
│
├── annotation/         # Genome annotation pipeline
│   ├── orchestrator.py # run_annotation_pipeline(), step ordering and concurrency
│   ├── orf.py          # ORF identification in transcript assemblies
│   ├── genemark.py     # GeneMark-ES/ET gene prediction
│   ├── alignment.py    # Protein alignment via spaln (intron-rich) or gth (intron-poor)
│   ├── augustus.py      # AUGUSTUS training and prediction
│   ├── snap.py         # SNAP and CodingQuarry gene prediction
│   ├── evm.py          # EVidenceModeler consensus building
│   ├── consensus.py    # Final model building: EVM + PASA UTRs + prettification
│   └── validation.py   # FASTA/GFF3 validation and genome header sanitization
│
├── assembly/           # Transcriptome assembly pipeline
│   ├── orchestrator.py # run_assembly() dispatch
│   ├── star.py         # STAR read mapping and hint generation
│   ├── trinity.py      # Trinity genome-guided and de novo assembly
│   └── pasa.py         # PASA spliced alignment and transcript hints
│
└── functional/         # Functional annotation pipeline
    ├── orchestrator.py # run_functional_annotation()
    ├── search.py       # pyhmmer phmmer/hmmscan search and result annotation
    └── dbfetch.py      # UniProt/Pfam database download and integrity tracking
```

### Pipeline Flow

1. Find ORFs in transcripts (if provided)
2. GeneMark gene prediction (ES or ET mode depending on RNA-seq hints)
3. Protein alignment via spaln (intron-rich) or gth (intron-poor)
4. AUGUSTUS training and prediction using protein + RNA-seq hints
5. SNAP training and prediction (fungus/protist also runs CodingQuarry)
6. Consensus model building via EVM, weighted by evidence type
7. Optional UTR addition via PASA
8. Final GFF3 formatting with locus tags

### Conventions

- GFF3 manipulations chain through `gffutils.create_db(':memory:', transform=fn)` passes
- External commands use `run_cmd(["cmd", "arg"], cwd=step_dir)` — never shell strings
- `--kingdom` flag (fungus/protist/animal/plant) controls which predictors run
- All pipeline state is in `PipelineConfig` (pydantic-settings) — no mutable class state
- Each annotation step lives in its own file under `annotation/`, one file per tool
- The three pipelines (annotation, assembly, functional) share the same structure: `orchestrator.py` + tool-specific files
