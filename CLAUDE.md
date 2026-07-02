# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Eukan is a eukaryotic genome annotation pipeline that integrates ab initio gene prediction (GeneMark, AUGUSTUS, SNAP, CodingQuarry) with homology-based evidence (protein alignments via spaln/gth) and transcript assemblies to produce consensus gene models via the external `combinr consensus` engine, which folds in UTRs and alternative isoforms from the transcript evidence. It optionally adds functional annotation via phmmer-against-UniProt or hmmscan-against-KOfam (adapted from KofamKOALA), plus hmmscan-against-Pfam.

## Build and Run

Uses Poetry for package management. The pipeline runs inside Docker with a GeneMark-ES/ET/EP+ license required.

```bash
# Install locally for development (pulls in pytest, ruff, mypy)
poetry install --with dev

# Run tests
poetry run pytest tests/ -v
# The combinr-consensus end-to-end test runs the external combinr binary; it
# auto-skips unless combinr is on PATH or $COMBINR_BIN points at the binary.

# CLI (all subcommands)
poetry run eukan --help
poetry run eukan mask-repeats -g genome.fasta
poetry run eukan annotate -g genome.fasta -p proteins.fasta --kingdom protist
poetry run eukan assemble -g genome.fasta -l left.fq -r right.fq -A -T --run-combinr
poetry run eukan func-annot -p proteins.faa --gff3 genes.gff3
poetry run eukan prep-submission -t submission.sbt --organism "Genus species"
poetry run eukan gff3toseq -g genome.fa -i genes.gff3 --output-format protein -o proteins.faa
poetry run eukan db-fetch -o databases/
poetry run eukan compare -r ref.gff3 -p pred.gff3                # single
poetry run eukan compare -r ref.gff3 -p p1.gff3 -p p2.gff3 -p p3.gff3 \
    -o details.tsv                                                 # multi-pred + per-feature TSV

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

### CLI (`eukan/cli/`)

Click-based CLI package with subcommands: `annotate`, `assemble`, `mask-repeats`, `func-annot`, `prep-submission`, `gff3toseq`, `db-fetch`, `check`, `status`, `compare`. One file per subcommand under `eukan/cli/`; shared infrastructure (option-group rendering, `numcpu_option`, `genome_option`, error formatting) lives in `eukan/cli/_framework.py`. Entry point defined in `pyproject.toml` as `eukan = "eukan.cli:cli"`. All subcommands use harmonized option groups ("Required input", "Pipeline parameters", "Re-run steps"). Step re-run flags follow the `--run-*` pattern (e.g., `--run-genemark`, `--run-star`).

### Package Structure

```
eukan/
├── cli/                # Click CLI package (one file per subcommand)
│   ├── _framework.py   # Shared Click helpers: option groups, common options, error formatting
│   ├── annotate.py     # eukan annotate
│   ├── assemble.py     # eukan assemble
│   ├── mask_repeats.py # eukan mask-repeats
│   ├── func_annot.py   # eukan func-annot
│   ├── prep_submission.py # eukan prep-submission
│   ├── compare.py      # eukan compare
│   ├── db_fetch.py     # eukan db-fetch
│   ├── gff3toseq.py    # eukan gff3toseq
│   ├── check.py        # eukan check (pre-flight tool/database checks)
│   └── status.py       # eukan status (manifest reader)
├── settings.py         # PipelineConfig, AssemblyConfig, RepeatsConfig, FunctionalConfig, SubmissionConfig (pydantic-settings)
├── validation.py       # FASTA/GFF3 validation and genome header sanitization
├── exceptions.py       # ConfigurationError and friends
│
├── infra/              # Runtime infrastructure (cross-pipeline)
│   ├── runner.py       # run_cmd(), run_piped(), run_parallel() — subprocess execution
│   ├── concurrency.py  # parallel_map() and friends
│   ├── manifest.py     # RunManifest, pipeline_step() — run tracking and reproducibility
│   ├── steps.py        # step_dir(), step_complete(), validate_or_raise() — step dir mgmt
│   ├── pipeline.py     # StepSpec + run_simple_pipeline() / run_orchestrated_step() driver
│   ├── layout.py       # PIPELINE_SUBDIRS, step_work_dir(), sibling_step_dir() — per-step run dir layout
│   ├── artifacts.py    # Artifact enum + find() — cross-pipeline artifact registry
│   ├── logging.py      # get_logger(), setup_logging(), md5_file()
│   ├── genome.py       # ContigIndex, FASTA helpers
│   ├── genetic_code.py # Genetic code table abstractions
│   ├── health.py       # Tool/database probes for `eukan check`
│   ├── tools_registry.py # Tool metadata loaded from data/tools.toml
│   ├── conda_env.py    # Conda env var setup at CLI startup
│   └── environ.py      # Env var helpers
│
├── gff/                # GFF3 format operations
│   ├── transforms.py   # Transform callbacks for gffutils.create_db(transform=fn)
│   ├── concordance.py  # Genomic interval operations (concordance, overlap, merging)
│   ├── intervals.py    # Lower-level interval primitives
│   ├── hierarchy.py    # gene>mRNA>CDS hierarchy fixers, prettify_gff3()
│   ├── normalize.py    # GFF3 cleanup before downstream tools (e.g., table2asn)
│   └── io.py           # featuredb2gff3_file(), count_gff3_features(), iter_assembled_sequences()
│
├── annotation/         # Genome annotation pipeline
│   ├── pipeline.py     # run_annotation_pipeline(), phase ordering, prediction-count logging
│   ├── orf.py          # ORF identification in transcript assemblies
│   ├── genemark.py     # GeneMark-ES/ET gene prediction
│   ├── alignment.py    # Protein alignment via spaln (intron-rich) or gth (intron-poor)
│   ├── spaln_params.py # Experimental --spsp species-specific spaln parameter builder
│   ├── augustus.py     # AUGUSTUS training and prediction
│   ├── training.py     # Training-set construction shared across predictors
│   ├── snap.py         # SNAP and CodingQuarry gene prediction
│   ├── evidence.py     # Evidence-source role mapping + weights helpers for the consensus engine
│   ├── combinr_consensus.py # Consensus engine: external `combinr consensus` (weighted DP; folds in UTRs + alt isoforms)
│   └── consensus.py    # Final model building: combinr consensus + ORF patch + prettification
│
├── assembly/           # Transcriptome assembly pipeline
│   ├── pipeline.py     # run_assembly() dispatch (StepSpec-driven); _DOWNSTREAM cascade
│   ├── star.py         # STAR read mapping + STARlong spliced transcript→genome mapping (map_transcripts)
│   ├── segemehl.py     # segemehl read/transcript mapping (fallback + non-canonical primary); _TRANSCRIPT_SETS
│   ├── align_hints.py  # RNA-seq intron/coverage hint emission from the read BAM
│   ├── bam_diagnostic.py # Post-STAR soft-clip + intron BAM walk → trans-splicing / non-canonical-splice verdict
│   ├── bam_introns.py  # Max-intron-bounded BAM filtering primitive
│   ├── trinity.py      # Trinity genome-guided + de novo assembly (the active assembler, both modes)
│   ├── tracks.py       # Single source of truth for the mapped Trinity track stems + per-track filenames
│   ├── jaccard.py      # In-house jaccard fusion-clip (replaces Trinity --jaccard_clip; STAR-based, tunable)
│   ├── strand_correction.py # Homology-calibrated per-transcript splice-strand flip (diamond blastx, opt-in --uniprot)
│   ├── defuse.py       # Homology-grounded chimera splitting (opt-in --defuse + --uniprot)
│   ├── max_intron.py   # Strand-agnostic max-intron split of transcript models (SL-independent; runs before sl_cut)
│   ├── sl_acceptors.py # Spliced-leader trans-splice acceptor detection (read-side + de novo)
│   ├── sl_cut.py       # Genomic SL cut of transcript models at trans-splice acceptors → {stem}.cut.gff3
│   ├── polya.py        # Poly-A characterization + unmapped de-novo transcript output
│   ├── combinr.py      # combinr-assemble consolidation → nr_transcripts.{fasta,gff3} + hints_rnaseq.gff
│   ├── rnaspades.py    # DORMANT: rnaSPAdes de novo (kept, unwired — Trinity replaced it)
│   ├── stringtie.py    # DORMANT: StringTie genome-guided (kept, unwired — Trinity replaced it)
│   └── sl_depletion.py # SL-motif primitives (FASTA depletion path retired for the genomic SL cut)
│
├── repeats/            # Repeat masking pipeline
│   ├── pipeline.py     # run_repeats() (StepSpec-driven)
│   ├── modeler.py      # RepeatModeler family-library construction
│   └── masker.py       # RepeatMasker softmasking + AUGUSTUS hint emission
│
├── functional/         # Functional annotation pipeline
│   ├── pipeline.py     # run_functional_annotation() — branches on homology_db
│   ├── search.py       # pyhmmer phmmer/hmmscan + UniProt/KOfam/Pfam GFF3/FASTA emitters
│   ├── kofam.py        # KOfam (KofamKOALA-style): ko_list parsing, EC extraction, per-KO threshold scoring
│   └── dbfetch.py      # UniProt/KOfam/Pfam database download and integrity tracking
│
├── submission/         # NCBI submission prep
│   ├── pipeline.py     # table2asn wrapper (eukan prep-submission)
│   └── cleanup.py      # GFF3 pre-clean to preserve attributes through table2asn
│
├── compare/            # Annotation comparison (eukan compare)
│   ├── engine.py       # compare_annotations() single-pred + compare_multiple() driver
│   ├── format.py       # Terminal report + per-feature TSV writer (single + multi)
│   └── models.py       # ComparisonResult, MultiComparisonResult, FeatureRecord
│
└── data/               # Static data shipped with the package
    ├── tools.toml      # External-tool registry (versions, probe commands, env hints)
    └── configs/        # AUGUSTUS config template
```

### Pipeline Flow

1. Find ORFs in transcripts (if provided) — respects configured genetic code
2. GeneMark gene prediction (ES or ET mode depending on RNA-seq hints) — passes `--gcode` for codes 6/26
3. Protein alignment via spaln (intron-rich) or gth (intron-poor)
   - Default: fitild intron length distribution → spaln `-yI`
   - `--spsp`: species-specific parameters via `make_eij.pl`/`make_ssp.pl` → spaln `-T` (experimental, uses separate `prot_align_ssp/` step dir)
4. AUGUSTUS training and prediction using protein + RNA-seq hints — auto-allows non-canonical splice sites from STAR evidence (`splice_site_summary.json`); `--splice-permissive` lowers thresholds
5. SNAP training and prediction (fungus/protist also runs CodingQuarry)
6. Consensus model building, weighted by evidence type, via the external `combinr consensus` engine (weighted DP), genetic-code aware, which folds UTRs and alternative isoforms in from transcript evidence in a single pass. Protein (`prot.gff3`, CDS-format) and transcript (`nr_transcripts.gff3`, flat exon) evidence are converted to `Target=` match chains; ab initio predictions and the weights file are reused as-is. Isoform-grouping stringency is tunable via `--combinr-stringent-overlap`; the binary is resolved via `--combinr-path` or PATH.
7. Final GFF3 formatting with locus tags
8. Optional functional annotation via `func-annot` (UniProt-or-KOfam plus Pfam, selected by `--homology-db`) → `final.mod.gff3`. KOfam mode is an adaptation of KofamKOALA: per-KO bit-score thresholds from `ko_list`, full vs domain score selection per KO, EC numbers parsed out of `[EC:…]` tags into a dedicated `ec_number=` GFF3 attribute, KEGG accessions emitted as `Dbxref=KEGG:K…`
9. Optional `prep-submission` runs NCBI's table2asn validator over `final.mod.gff3` to produce a `.sqn` plus `.val/.dr/.stats` reports for iterative GFF3 refinement

### Per-step Run Directory Layout

Each subcommand resolves its own `work_dir` to a sibling subdir under the user's cwd, declared in `eukan/infra/layout.py::PIPELINE_SUBDIRS`:

```
<cwd>/
├── repeats/        # eukan mask-repeats
├── assemble/       # eukan assemble
├── annotate/       # eukan annotate
├── func-annot/     # eukan func-annot
└── submission/     # eukan prep-submission
```

Cross-pipeline artifact lookups go through `eukan/infra/artifacts.py`:

- `Artifact` (StrEnum) lists every file that crosses pipeline boundaries (e.g. `RNASEQ_HINTS`, `FINAL_GFF3`, `FINAL_FUNC_GFF3`).
- `_PRODUCER` declares which step writes each artifact (e.g. `FINAL_FUNC_GFF3 → "func-annot"`).
- `find(work_dir, Artifact.X)` returns the first existing match across the caller's own work_dir, then the sibling step dir. This is what lets `prep-submission` auto-discover `func-annot/final.mod.gff3` from `submission/` without an explicit path.

When adding a new cross-pipeline file, register it in both `Artifact` and `_PRODUCER` rather than hardcoding paths.

### Run Manifest

All pipelines share a single `eukan-run.json` per run dir. Step names are prefixed by pipeline (`annotation/genemark`, `assembly/star`, `functional/search`, `repeats/modeler`, `submission/table2asn`). The manifest tracks per-step status, timing, and output checksums for resume and integrity checking. Key functions: `get_or_create_manifest()`, `pipeline_step()`, `is_step_complete()` in `infra/manifest.py`. Each config has a `manifest_dir` field (defaults to `work_dir`) controlling where the manifest is written.

### Pipeline Driver

`infra/pipeline.py` provides two entry points:

- `run_simple_pipeline(steps: list[StepSpec], …)` — linear case used by `assemble` and `mask-repeats`. Each `StepSpec` declares name, function, output filename, and re-run flag display string.
- `run_orchestrated_step(work_dir, manifest, step_key, fn, *args, …)` — lower-level primitive used directly by the annotation and functional pipelines, whose execution graph isn't a straight line (annotation has fan-out phases, functional caches JSON between steps).

Both wrap step-dir setup, manifest updates, and force/skip logic so that step-level resume and integrity checking are uniform across pipelines.

### Conventions

- GFF3 manipulations chain through `gffutils.create_db(':memory:', transform=fn)` passes
- External commands use `run_cmd(["cmd", "arg"], cwd=step_dir)` — never shell strings
- `--kingdom` flag (fungus/protist/animal/plant) controls which predictors run
- All pipeline state is in a pydantic-settings config (`PipelineConfig`, `AssemblyConfig`, etc.) — no mutable class state
- Each pipeline package follows the same shape: `pipeline.py` (driver) + one file per tool
- Per-step output isolation: every subcommand writes under its own `<cwd>/<step-subdir>/`. Don't write into `cwd` directly — use `step_work_dir(step)` from `infra/layout.py`
- Cross-pipeline files go through `Artifact` + `find()` rather than hardcoded paths
- CLI option groups are harmonized across subcommands: "Required input", "Pipeline parameters", "Re-run steps". Step re-run flags follow the `--run-*` pattern (e.g., `--run-genemark`, `--run-star`). See "CLI conventions" below for the full canonical layout.
- Prediction-count logging: every per-tool step that emits a GFF3 calls `_log_prediction_count("Tool", path)` (annotation pipeline) or `count_gff3_features()` so the user sees gene counts at each phase
- 3-way concordance (`gff/concordance.extract_supported_models`) emits per-source counts at INFO and a WARNING when concordance falls below `WEAK_CONCORDANCE_THRESHOLD` (250 gene models) — used by both AUGUSTUS and SNAP training-set construction

### CLI conventions

Authoritative reference for shared flag spellings and option-group layout. New subcommands should follow this template.

**Canonical short flags** (the same letter means the same thing in every subcommand that uses it):

| Flag | Long form | Notes |
|------|-----------|-------|
| `-g` | `--genome` | Use `genome_option()` from `_framework.py`. |
| `-i` | `--gff3` | Input GFF3 file. |
| `-p` | `--proteins` | And `--predicted` in `compare` (same letter, related concept). |
| `-n` | `--numcpu` | Use `numcpu_option` from `_framework.py`. |
| `-f` | `--force` | **Reserved.** Use `force_option` (or `force_option(help_text=...)` for custom wording). |
| `-c` | `--code` | NCBI genetic code. Use `code_option(default=N)` from `_framework.py`. |
| `-o` | `--output-file` / `--output-dir` | Output destination — pick `-file` for a path, `-dir` for a directory. Never bare `--output`. |
| `-d` | directory option | `--db-dir`, `--work-dir`, `--output-dir`. |
| `-k` | `--kingdom` | annotate only. |
| `-l` / `-r` / `-s` | `--left` / `--right` / `--single` | assemble paired/single reads. |
| `-r` | `--reference` (compare), `--rnaseq-hints` (annotate) | Different commands, no real collision. |
| `-e` | `--evalue` | func-annot. |
| `-w` | `--weights` | annotate. |
| `-t` | `--template` (prep-submission), `--align-mode` (assemble) | Different commands. |

**Option groups**, in display order:

1. `Required input` — required paths (genome, gff3, proteins, reads).
2. `Source qualifiers` — prep-submission only (`--organism`, `--isolate`, `--source-info`, `--locus-tag-prefix`).
3. `Pipeline parameters` — tunables that affect behaviour but have defaults.
4. `Override options` — flags that override auto-discovered values from prior pipeline runs.
5. `Experimental` — opt-in flags for unstable / under-evaluation features.
6. `Output options` — output paths and dry-run/print-only modifiers.
7. `Re-run steps` — `--run-*` per-step flags plus the `-f, --force` master flag at the end.

**Shared helpers** (`eukan/cli/_framework.py`): always use these instead of redefining inline.

- `genome_option(help_text=...)` — required `--genome/-g`.
- `numcpu_option` — `--numcpu/-n` with `cpu_count()` default.
- `force_option` or `force_option(help_text="…")` — `--force/-f` flag.
- `code_option(default=N)` — `--code/-c` int with command-specific default.

Always use `@click.command(cls=PreformattedEpilogCommand, ...)` when the command uses option groups so the help text renders without the wrapping "Options:" header.
