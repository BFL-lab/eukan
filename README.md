# eukan: Eukaryotic Genome Annotation Pipeline

A comprehensive annotation pipeline tailored for eukaryotic genomes, particularly those from less well-studied organisms.

## Installation

Currently, Eukan installation is only supported via Docker and Conda.

### Docker

The Docker image installs all external tools via conda (from the same `environment.yml` used for local installs), then builds fitild from source and optionally includes GeneMark.

```bash
git clone https://github.com/BFL-lab/eukan.git
cd eukan
docker build -t eukan -f docker/Dockerfile .
```

To include **GeneMark** ([license required](https://topaz.gatech.edu/GeneMark/license_download.cgi)), place `gmes_linux_64*.tar.gz` and `gm_key_64.gz` in the project root before building. If omitted, the build will succeed but `eukan check` will report GeneMark as missing.

A separate **development image** adds test tooling (NCBI datasets CLI, procps):

```bash
docker build -t eukan-dev -f docker/Dockerfile.dev .
```

### Conda

Installs all external tools via bioconda and eukan itself via pip, in one step:

```bash
git clone https://github.com/BFL-lab/eukan.git
cd eukan
conda env create -f environment.yml
conda activate eukan
eukan check
```

The `eukan` CLI configures all required environment variables (e.g. `$ZOE`, `$ALN_TAB`, EVM paths) automatically at startup. If you need to run the underlying tools directly outside of `eukan`, install the optional activation script:

```bash
mkdir -p $CONDA_PREFIX/etc/conda/activate.d
cp conda-activate.sh $CONDA_PREFIX/etc/conda/activate.d/eukan.sh
conda deactivate && conda activate eukan
```

The `environment.yml` is auto-generated from `tools.toml`. To regenerate after modifying tool versions: `python scripts/generate-env.py`.

Two tools require manual installation after creating the environment. A helper script handles both:

```bash
# Download GeneMark (license required) from:
#   https://topaz.gatech.edu/GeneMark/license_download.cgi
# Place gmes_linux_64*.tar.gz and gm_key_64.gz in the current directory, then:
./scripts/install-extras.sh

# Or point to a GeneMark archive elsewhere:
./scripts/install-extras.sh --genemark-tar /path/to/gmes_linux_64_4.tar.gz
```

`eukan check` will tell you exactly what's missing and how to install it.

### Local development

Requires Python >= 3.10 and [Poetry](https://python-poetry.org/). External tools must be installed separately (via conda or Docker):

```bash
poetry install
poetry run eukan --help
```

### Dependencies

**Python** (managed by Poetry): click, gffutils, biopython, pandas, requests, pydantic-settings.

**External tools** (via Docker image or conda): AUGUSTUS, SNAP, CodingQuarry, spaln, GenomeThreader, EVidenceModeler, PASA, Trinity, STAR, samtools, BLAT, jellyfish, GMAP/GSNAP, fasta36, TRF.

**Manual install**:
- GeneMark-ES/ET/EP+:  [license required](https://topaz.gatech.edu/GeneMark/license_download.cgi)
- fitild:  [github.com/ogotoh/fitild](https://github.com/ogotoh/fitild) (only needed for spaln-based protein alignment in default fitild mode)
- spaln analysis utilities (`utn`, `npssm`, `exinpot`, etc.): built from [spaln source](https://github.com/ogotoh/spaln) via `make all` (only needed for `--spsp` species-specific parameter mode)

## Usage

### Docker

Use `eukan-docker` as a wrapper to run any subcommand inside the Docker container. It bind-mounts the current directory and runs as your user:

The typical workflow runs each subcommand from the same working directory, so later steps auto-discover outputs from earlier ones:

```bash
# 1. Download reference databases
./eukan-docker db-fetch

# 2. Assemble transcriptome from RNA-seq reads (optional but recommended)
./eukan-docker assemble \
    -g genome.fasta \
    -l left_reads.fastq -r right_reads.fastq \
    -S RF --kingdom protist

# 3. Annotate:  auto-discovers assembly outputs + databases from steps above
./eukan-docker annotate -g genome.fasta -p proteins.fasta --kingdom protist

# 4. Functional annotation:  auto-discovers predicted proteins + databases
./eukan-docker func-annot

# Extract sequences from GFF3
./eukan-docker gff3toseq -g genome.fasta -i genes.gff3 -o protein
```

Each auto-discovered input can be overridden explicitly (see subcommand docs below).

Set `EUKAN_IMAGE` to use a custom image name (default: `eukan`):

```bash
EUKAN_IMAGE=myregistry/eukan:latest ./eukan-docker annotate ...
```

### Local (development)

```bash
poetry run eukan annotate -g genome.fasta -p proteins.fasta --kingdom fungus
poetry run eukan assemble -g genome.fasta -l left.fq -r right.fq -S RF
```

## Subcommands

### `eukan annotate`

Run the genome annotation pipeline. When run in the same directory as `eukan assemble`, transcript evidence (FASTA, GFF3, RNA-seq hints), strand-specificity, and a PASA database for UTR addition are discovered automatically.

```
Usage: eukan annotate [OPTIONS]

Required input:
  -g, --genome PATH               Genome FASTA (no lower-case; pipeline soft-masks repeats). [required]
  -p, --proteins PATH             One or more protein FASTA files. [required]

Pipeline parameters:
  -k, --kingdom [fungus|protist|animal|plant]
                                   Target organism kingdom.
  -n, --numcpu INTEGER             Number of CPU threads. [default: all]
  --existing-augustus TEXT          Use pre-trained AUGUSTUS species parameters.
  -w, --weights INTEGER            Weights: protein, gene predictions, transcripts. [default: 2 1 3]
  -C, --code INTEGER               NCBI genetic code table number. [default: 11]

Override options:
  -tf, --transcripts-fasta PATH   Override auto-discovered transcript FASTA.
  -tg, --transcripts-gff PATH     Override auto-discovered transcript GFF3.
  -r, --rnaseq-hints PATH         Override auto-discovered RNA-seq hints GFF.
  --strand-specific                Transcripts are strand-oriented.
  --utrs PATH                      PASA SQLite database for UTR addition.
  --splice-permissive              Allow non-canonical splice sites (GC-AG, AT-AC).

Experimental:
  --spsp                           Build species-specific spaln parameters from transcripts
                                   (alternative to fitild). See "Protein alignment modes" below.

Re-run steps:
  --run-genemark                   Force re-run GeneMark gene prediction.
  --run-prot-align                 Force re-run protein alignment (spaln/gth).
  --run-augustus                    Force re-run AUGUSTUS training and prediction.
  --run-snap                       Force re-run SNAP (and CodingQuarry) prediction.
  --run-consensus                  Force re-run EVM consensus model building.
```

### `eukan assemble`

Assemble transcriptome from RNA-seq reads for use with `eukan annotate`. Provide either paired-end reads (`--left` and `--right`) or single-end reads (`--single`).

```
Usage: eukan assemble [OPTIONS]

Required input:
  -g, --genome PATH               Genome FASTA. [required]
  -l, --left PATH                 Left paired-end reads.
  -r, --right PATH                Right paired-end reads.
  -s, --single PATH               Single-end reads.

Pipeline parameters:
  -n, --numcpu INTEGER             Number of CPUs. [default: all]
  -S, --strand-specific [RF|FR|R|F]
                                   Strand-specific library type (RF/FR for paired, R/F for single).
  -t, --align-mode [EndToEnd|Local] STAR alignment mode. [default: Local]
  --splice-permissive              Allow non-canonical splice sites (GC-AG, AT-AC).
  -c, --genetic-code [1|6|10|12]   Genetic code for PASA. [default: 1]
  -m, --min-intron INTEGER         Min intron length. [default: 20]
  -M, --max-intron INTEGER         Max intron length. [default: 5000]
  --phred [33|64]                  Phred quality score. [default: 33]
  -j, --jaccard-clip               Enable jaccard clipping.

Re-run steps:
  -A, --run-star                   Force re-run STAR read mapping.
  -T, --run-trinity                Force re-run Trinity assembly.
  -P, --run-pasa                   Force re-run PASA alignment.
  -f, --force                      Force re-run all steps.
```

The pipeline runs STAR mapping, genome-guided + de novo Trinity assembly, and PASA alignment. STAR also profiles splice site types from junction evidence (`splice_site_summary.json`), which the annotation pipeline uses to allow non-canonical splice sites in AUGUSTUS. If no step flags (`-A`, `-T`, `-P`) are given, all steps run.

### `eukan func-annot`

Add functional annotations (UniProt + Pfam) to proteins. When run after `eukan annotate` and `eukan db-fetch`, the predicted protein sequences, UniProt, and Pfam databases are discovered automatically.

```
Usage: eukan func-annot [OPTIONS]

Pipeline parameters:
  -n, --numcpu INTEGER   Number of CPUs. [default: all]
  -e, --evalue TEXT      E-value cutoff. [default: 1e-1]

Override options:
  -p, --proteins PATH    Amino acid sequences FASTA.
  --uniprot PATH         UniProt-SwissProt database FASTA.
  --pfam PATH            Pfam HMM database.
  --gff3 PATH            GFF3 file to annotate with functional info.
  -f, --force            Re-run steps even if outputs exist.
```

Runs phmmer against UniProt and hmmscan against Pfam (via pyhmmer). Produces:
- `input.mod.faa`: annotated FASTA with functional descriptions in headers.
- `input.mod.gff3`: (if `--gff3` provided) GFF3 with `product` and `inference` attributes.

Hits with e-values between 1e-3 and the cutoff are reported as marginal.

### `eukan gff3toseq`

Extract protein or cDNA sequences from a GFF3 + genome.

```
Usage: eukan gff3toseq [OPTIONS]

Options:
  -g, --genome PATH           Genome FASTA. [required]
  -i, --gff3 PATH             GFF3 with gene models. [required]
  -o, --output [protein|cdna] Output type. [default: protein]
  -c, --code INTEGER          Genetic code table. [default: 1]
```

### `eukan db-fetch`

Download reference databases (UniProt, Pfam).

```
Usage: eukan db-fetch [OPTIONS]

Options:
  -o, --output-dir PATH   Directory to download into. [default: databases]
  -f, --force             Re-download even if databases are up to date.
  -d, --database [uniprot|pfam]
                          Specific database(s) to fetch. If omitted, fetch all.
```

Downloads and prepares:
- `uniprot_sprot.faa`:  UniProt-SwissProt protein sequences (converted from XML).
- `Pfam-A.hmm`:  Pfam HMM profiles (decompressed and pressed for hmmscan).

### `eukan compare`

Compare predicted gene models against a reference or previous annotation to assess annotation quality. Reports gene-level classification (exact, inexact, missing, merged, fragmented, novel), subfeature-level metrics (mRNA, CDS, intron), and overlap-based sensitivity/specificity/F1 scores.

```
Usage: eukan compare [OPTIONS]

Required input:
  -r, --reference PATH    Reference GFF3 file.
  -p, --predicted PATH    Predicted GFF3 file to evaluate.

Output options:
  -o, --output-file PATH  Write per-feature details to a TSV file.
```

The classification system and metrics are further described in the paper referenced in [Citation](#citation). Gene-level classifications:
- **exact**: prediction coordinates match reference exactly.
- **inexact**: prediction overlaps a single reference with boundary differences.
- **missing**: reference gene with no overlapping prediction (false negative).
- **merged**: prediction spans 2+ reference genes.
- **fragmented**: 2+ predictions each cover a single reference gene.
- **novel**: prediction with no reference overlap (possibly false positive).

For matched features, reports overlap-based sensitivity (overlap / reference length), specificity (overlap / prediction length), and F1 score. Boundary differences (5' and 3') are reported for inexact matches.

```bash
# Basic comparison
eukan compare -r reference.gff3 -p predicted.gff3

# With per-feature TSV output
eukan compare -r reference.gff3 -p predicted.gff3 -o details.tsv
```

### `eukan check`

Verify Python dependencies, external tools, and databases.

```
Usage: eukan check [OPTIONS]

Options:
  --for [annotate|assemble|func-annot|db-fetch]
      Only check tools needed by these subcommands. If omitted, check all.
  --db-dir PATH   Database directory to check. [default: databases]
```

Checks Python dependencies, probes each external tool with a version/help command, and verifies database integrity. Exits 0 if all checks pass, 1 if any fail.

```bash
# Check everything
eukan check

# Check only what's needed for annotation
eukan check --for annotate

# Check multiple subcommands
eukan check --for annotate --for assemble
```

Example output:
```
Checked 14 external tools:

  12 tools OK:
    ✓ samtools                       samtools 1.20
    ✓ AUGUSTUS                        AUGUSTUS (3.5.0)
    ...

  2 tools MISSING or BROKEN:
    ✗ codingquarry                   `CodingQuarry` not found on PATH; env not set: $QUARRY_PATH
      used by: annotate
    ✗ fitild                         `fitild` not found on PATH
      used by: annotate
      hint: Build from source: git clone https://github.com/ogotoh/fitild ...
```

## Pipeline Overview

The annotation pipeline (`eukan annotate`) runs the following steps:

1. **ORF finding**:  Identify ORFs in transcript assemblies (if provided). Uses the configured genetic code (`-C`/`--code`) so that alternative stop codons (e.g., code 6 where TAA/TAG encode glutamine) are handled correctly.
2. **GeneMark**:  Ab initio gene prediction (ES mode, or ET mode if RNA-seq intron hints are available with >= 150 introns). Passes `--gcode` for non-standard genetic codes (codes 6 and 26).
3. **Protein alignment**:  Spliced alignment via spaln (intron-rich genomes, > 25% introns/gene) or GenomeThreader (intron-poor). See [Protein alignment modes](#protein-alignment-modes).
4. **AUGUSTUS**:  Train species-specific parameters from concordant GeneMark/protein models, then predict genes using protein + RNA-seq hints. Non-canonical splice sites (e.g., AT-AC) are allowed automatically when supported by sufficient junction evidence from STAR; `--splice-permissive` lowers the evidence thresholds.
5. **SNAP**:  Train and predict (all kingdoms). **CodingQuarry** also runs for fungus/protist genomes.
6. **EVidenceModeler**:  Build weighted consensus gene models from all evidence sources.
7. **PASA UTRs**:  Add UTR annotations and model alternative splicing from the transcriptome database.
8. **Final output**:  Assign sequential locus tags and correct CDS phases. Non-overlapping transcript ORFs not captured by EVM are patched into the final model set.

Output: `final.gff3` in the working directory.

### Protein alignment modes

spaln protein alignment supports two modes for modeling intron structure:

**Default (fitild)**: Builds an intron length distribution from GeneMark predictions and feeds it to spaln via the `-yI` parameter. This is the established approach and requires the [fitild](https://github.com/ogotoh/fitild) tool.

**Species-specific parameters (`--spsp`)**: Uses transcript data to build full species-specific spaln parameters (splice site models, intron potential, codon usage) via spaln's `make_eij.pl` and `make_ssp.pl` scripts. This produces richer alignment parameters than fitild alone, at the cost of additional computation. Requires transcript data (from `eukan assemble` or provided via `--transcripts-fasta`).

```bash
# Default mode (fitild)
eukan annotate -g genome.fasta -p proteins.fasta

# Species-specific parameter mode (experimental)
eukan annotate -g genome.fasta -p proteins.fasta --spsp
```

When `--spsp` is used, protein alignment results are written to a separate directory (`prot_align_ssp/`) so both modes can coexist for comparison.

### Run tracking and resume

All pipelines write to a shared `eukan-run.json` manifest in the working directory, tracking per-step status, timing, and output checksums. This enables:

- **Resume**: Re-running a subcommand skips steps that already completed.
- **Selective re-run**: Use `--run-*` flags to force specific steps to re-execute.
- **Integrity checking**: On startup, completed steps are validated (output exists and is non-empty). If a discrepancy is found, the pipeline reports the issue and suggests the `--run-*` flag to fix it.
- **Progress monitoring**: `eukan status` reads the manifest and shows progress across all pipelines.

```bash
# View run status
eukan status

# Re-run only protein alignment
eukan annotate -g genome.fasta -p proteins.fasta --run-prot-align

# Re-run only PASA in assembly
eukan assemble -g genome.fasta -l left.fq -r right.fq -P
```

## Testing

### Unit tests

```bash
poetry run pytest tests/ -v
```

Unit tests cover GFF3 parsing/serialization, genomic interval operations, ORF finding, configuration validation, run manifest tracking, database integrity checks, and CLI wiring. They run without external tools or network access.

### Pipeline integration test

A development CLI at `tests/run_pipeline.py` drives a full end-to-end pipeline run using _S. pombe_ chromosome III as the test organism.

#### Prerequisites

- All external tools installed (Docker or conda environment; verify with `eukan check`)
- [NCBI datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/) and [SRA Toolkit](https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit) on PATH (for downloading test data)

When using Docker, build and use the dev image (`eukan-dev`) which includes the NCBI datasets CLI:

```bash
docker build -t eukan-dev -f docker/Dockerfile.dev .
```

#### 1. Download test data

```bash
python tests/run_pipeline.py setup-test-data [-o tests/data]
```

Downloads from NCBI:
- **Genome**: _S. pombe_ chromosome III (`NC_003424.3`, ~2.5 Mbp)
- **Proteins**: 10 close neighbor proteomes via `datasets`
- **RNA-seq**: 5 SRA paired-end runs via `prefetch` + `fasterq-dump`

Accession lists live in `tests/data/*.txt` and are never deleted by cleanup. Downloads are idempotent.

#### 2. Run the pipeline

```bash
# Full run: assembly + annotation (default kingdom: fungus)
python tests/run_pipeline.py test-pipeline --kingdom fungus -n 8

# Protein-only: skip transcriptome assembly
python tests/run_pipeline.py test-pipeline --protein-only -n 8

# Custom directories
python tests/run_pipeline.py test-pipeline -d tests/data -w tests/pipeline-run
```

The test pipeline runs:
1. **Transcriptome assembly**: STAR read mapping, Trinity (genome-guided + de novo), PASA alignment
2. **Genome annotation**:  GeneMark, protein alignment (spaln/gth), AUGUSTUS, SNAP, EVM consensus

If assembly fails, it falls back to protein-only annotation automatically.

Output lands in `tests/pipeline-run/` with subdirectories for assembly and annotation. View run details with `eukan status -d tests/pipeline-run/annotation`.

#### 3. Clean up

```bash
# Remove pipeline outputs only
python tests/run_pipeline.py clean-test-data

# Remove outputs + downloaded data (genome, proteins, reads)
python tests/run_pipeline.py clean-test-data --all
```

All three subcommands accept `-h` for detailed help.

## Citation

If you use eukan, please cite:

> Sarrasin M, Burger G, Lang BF. Eukan: a fully automated nuclear genome annotation pipeline for less studied and divergent eukaryotes. *NAR Genomics and Bioinformatics*. 2026 Mar;8(1):lqag003. doi:[10.1093/nargab/lqag003](https://doi.org/10.1093/nargab/lqag003)

A [CITATION.cff](CITATION.cff) file is included for automated citation tools.

## License

See [LICENSE.md](LICENSE.md).
