# eukan

A eukaryotic genome annotation pipeline for less well-studied organisms

## Installation

Installation is best handled in docker.

First requires a license to run GeneMark-ES/ET/EP+ from topaz.gatech.edu/GeneMark/license_download.cgi

```
docker build -t eukan -f Dockerfile .
```

## Running

Use the helper shell script to run the container.

```
bash annot-docker eukan -h

usage: eukan [-h] --genome genome.fasta --proteins PROTEINS [PROTEINS ...] [--transcriptsFasta transcriptassembly.fasta] [--transcriptsGFF transcriptassembly.gff3] [--rnaseq_hints hints.gff] [--existing_augustus species]
                    [--strand_specific_transcripts] [--numcpu N] [--weights x y [z] [x y [z] ...]] [--code CODE] [--utrs UTRS] [--fungus] [--protist] [--animal] [--plant]

Annotates a genome.

optional arguments:
  -h, --help            show this help message and exit
  --genome genome.fasta, -g genome.fasta
                        REQUIRED. Make sure there are no lower-case letters in the sequences since the pipeline soft-masks the genome by converting upper-case nucleotides in repetitive regions to lower-case
  --proteins PROTEINS [PROTEINS ...], -p PROTEINS [PROTEINS ...]
                        REQUIRED. >=1 protein fasta files, separated by spaces
  --transcriptsFasta transcriptassembly.fasta, -tf transcriptassembly.fasta
                        assembled transcripts fasta
  --transcriptsGFF transcriptassembly.gff3, -tg transcriptassembly.gff3
                        assembled transcripts gff file
  --rnaseq_hints hints.gff, -r hints.gff
                        gff hints generated from RNA-Seq
  --existing_augustus species
                        use existing augustus species parameters
  --strand_specific_transcripts
                        specify that transcripts are alignment oriented
  --numcpu N, -n N      number of CPUs to use, default is MAX
  --weights x y [z] [x y [z] ...]
                        define weights for each evidence source in the following order: protein alignments, gene predictions, transcript assembly (if present), e.g. --weights 2 1 5, or --weights 2 1. Default is 1 2, and 10 if transcript
                        assembly is available
  --code CODE, -C CODE  genetic code to use, as defined at https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG27
  --utrs UTRS           path to PASA sqlite DB to add UTRs
  --fungus              target genome is of a fungus
  --protist             target genome is of a protist
  --animal              target genome is of an animal
  --plant               target genome is of a plant

```

### Data prep

Data needs to be prepared prior to running Eukan. Use the `transcriptome_assembly.sh` script to prepare RNA-Seq read libraries:

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
        [-n] <number of CPUs> # default is MAX
        [-S] <specificy strand-specific assembly, either RF or FR> # default off, i.e. unstranded
        [-A] <switch on read mapping>
        [-E] <switch to extract reads>
        [-T] <switch on Trinity assembly>
        [-e] <switch on StringTie assembly>
        [-P] <switch on PASA alignment>
        [-c] <genetic code according to ncbi table>
        [-h] Display this help message
        [-j] switch on jaccard clipping (for gene-dense orgnanisms and high coverage data)
        [-t] <EndToEnd/Local> # default Local
```

The pipeline starts by mapping reads to the genome, then runs the transcriptome assembly routine, then creates a non-redundant assembly by aligning the resulting transcripts to the genome.
