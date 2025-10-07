# eukan: Eukaryotic Genome Annotation Pipeline

A comprehensive annotation pipeline tailored for eukaryotic genomes, particularly those from less well-studied organisms. The pipeline integrates evidence-based gene prediction, homology mapping, and functional annotation to produce genome annotations that can be easily passed to downstream submission tools.

## Features

- **Genome Annotation**: Combines ab initio gene prediction with homology-based evidence to annotate protein-coding genes in eukaryotic genomes.
- **Evidence Integration**: Incorporates protein alignments, transcript assemblies, and RNA-seq hints for accurate gene models.
- **Transcript Assembly**: Preprocesses RNA-seq data into transcript assemblies for improved annotation accuracy.
- **Functional Annotation**: Annotates predicted proteins with functional information from UniProt-SwissProt and Pfam databases.
- **Flexible Configuration**: Supports various kingdoms (fungus, protist, animal, plant) and custom genetic codes.
- **Docker Deployment**: Simplifies installation and usage through containerization.

## Installation

The pipeline requires Docker for isolated and reproducible execution. Before building the image, obtain a license for GeneMark-ES/ET/EP+ from [topaz.gatech.edu/GeneMark/license_download.cgi](topaz.gatech.edu/GeneMark/license_download.cgi).

### Building the Docker Image

```bash
git clone https://github.com/BFL-lab/eukan.git
cd eukan
docker build -t eukan -f Dockerfile .
```

### Dependencies

- Docker
- GeneMark-ES/ET/EP+ license
- For functional annotation: HMMER suite, Python 3 with biopython and other dependencies

## Usage

### Genome Annotation

Use the provided `annot-docker` script as a wrapper to run the pipeline inside the Docker container.

```bash
# Display help
./annot-docker eukan -h
```

#### Basic Command

```bash
./annot-docker eukan --genome genome.fasta --proteins protein1.faa protein2.faa
```

#### Full Usage

```
usage: eukan [-h] --genome genome.fasta --proteins PROTEINS [PROTEINS ...] [--transcriptsFasta transcriptassembly.fasta] [--transcriptsGFF transcriptassembly.gff3] [--rnaseq_hints hints.gff] [--existing_augustus species]
                    [--strand_specific_transcripts] [--numcpu N] [--weights x y [z] [x y [z] ...]] [--code CODE] [--utrs UTRS] [--fungus] [--protist] [--animal] [--plant]

Annotates a eukaryotic genome.

required arguments:
  --genome genome.fasta       REQUIRED. Genome sequence in Fasta format. Ensure no lower-case nucleotides; the pipeline soft-masks repeats by converting to lower-case.
  --proteins PROTEINS [PROTEINS ...]
                               REQUIRED. One or more protein sequence Fasta files, separated by spaces.

optional arguments:
  --transcriptsFasta transcriptassembly.fasta
                               Assembled transcripts in Fasta format.
  --transcriptsGFF transcriptassembly.gff3
                               Assembled transcripts in GFF3 format.
  --rnaseq_hints hints.gff      GFF hints file generated from RNA-seq alignment.
  --existing_augustus species   Use pre-trained AUGUSTUS species parameters.
  --strand_specific_transcripts
                               Specify that assembled transcripts are strand-oriented.
  --numcpu N                    Number of CPU threads to use (default: all available).
  --weights x y [z] [x y [z] ...]
                               Weights for scoring evidence sources: protein alignments, gene predictors, transcript assembly (if provided).
                               Default: 1 2, plus 10 if transcript assembly is included.
  --code CODE                  Genetic code (see NCBI taxonomy utils).
  --utrs UTRS                  PASA SQLite database path for adding UTRs.
  --fungus                     Tune parameters for fungal genomes.
  --protist                    Tune parameters for protist genomes.
  --animal                     Tune parameters for animal genomes.
  --plant                      Tune parameters for plant genomes.
```

### Transcriptome Assembly

Prepare RNA-seq data for input using the `transcriptome_assembly.sh` script. This handles read mapping, assembly, and alignment to produce input files for the main pipeline.

```bash
# Display help
transcriptome_assembly.sh -h
```

#### Usage

```
Usage: transcriptome_assembly.sh [OPTIONS] <ARGS>

        [OPTIONS] and corresponding <ARGS> are:

        Either paired-end:
                [-l] <left reads>
                [-r] <right reads>
        or single-end:
                [-s] <single-end reads>
        [-m] <min intron length> # default 20
        [-M] <max intron length> # default 5000
        [-g] <genome fasta>
        [-p] <phred quality score (33 for MISEQ, 64 for HISEQ)> # default 33
        [-n] <number of CPUs> # default MAX
        [-S] <specify strand-specific assembly, either RF or FR> # default off
        [-A] <switch on read mapping>
        [-E] <switch to extract reads>
        [-T] <switch on Trinity assembly>
        [-e] <switch on StringTie assembly>
        [-P] <switch on PASA alignment>
        [-c] <genetic code according to ncbi table>
        [-h] Display this help message
        [-j] switch on jaccard clipping (for gene-dense organisms and high coverage data)
        [-t] <EndToEnd/Local> # default Local
```

#### Example

```bash
# Assembled paired-end reads with Trinity
transcriptome_assembly.sh -l left_reads.fastq -r right_reads.fastq -g genome.fasta -M 10000 -S RF -A -T -P
```

The pipeline integrates genome mapping, de novo assembly (Trinity), followed by PASA alignment for evidence integration.

## Functional Annotation

The `functional-annotation` directory contains scripts to add functional information to predicted proteins using similarity searches against UniProt-SwissProt and Pfam databases.

### Prerequisites

1. **Databases**: Prepare UniProt-SwissProt and Pfam databases (or use defaults if available).
2. **Dependencies**: Install HMMER, Python 3, and required packages:
   - Python packages: `biopython`, `requests`, `gffutils` (see `requirements.txt`)

   ```bash
   cd functional-annotation
   pip install -r requirements.txt
   ```

3. **Get Databases**: Use `db-fetch.py` to download and format the latest databases:

   ```bash
   python db-fetch.py
   ```

   This will download:
   - `uniprot_sprot.faa`: UniProt-SwissProt protein sequences
   - `Pfam-A.hmm`: Pfam HMM profiles (pressed for hmmscan)

### Running Functional Annotation

The `func-annot` script runs `phmmer` against UniProt and `hmmscan` against Pfam to annotate protein sequences. Results are appended to Fasta headers or GFF3 attributes.

#### Usage

```bash
func-annot --proteins input.faa [--uniprot uniprot_sprot.faa] [--pfam Pfam-A.hmm] [--gff3 input.gff3] [--numcpu N] [--evalue 1e-5]
```

#### Full Arguments

- `--proteins PROTEINS, -p PROTEINS`: Amino acid sequences in Fasta format (required).
- `--uniprot uniprot_sprot.faa`: UniProt-SwissProt database (default: `/share/unsupported/databases/uniprot_sprot/uniprot_sprot.faa`).
- `--pfam Pfam-A.hmm`: Pfam HMM database (default: `/share/unsupported/databases/Pfam/35.0/Pfam-A.hmm`).
- `--gff3 gene_models.gff3`: Optional GFF3 file to annotate with functional information.
- `--numcpu N, -n N`: Number of CPUs (default: all).
- `--evalue Me-N, -e Me-N`: E-value cutoff (default: 1e-1; marginal hits: 1e-3 to 1e-1).

#### Output

- Annotated Fasta: `input.mod.faa` with functional descriptions in headers.
- Optional: `input.mod.gff3` with added `product` and `inference` attributes.

#### Examples

```bash
# Run functional annotation pipeline and append information to fasta headers
func-annot -p input.faa

# Append functional info to Fasta headers from stricter e-values, and update corresponding gff3 feature column with annotations that can be read by table2asn
func-annot -p input.faa --evalue 1e-5 --gff3 input.gff3
```
