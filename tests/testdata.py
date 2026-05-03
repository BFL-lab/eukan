"""Test data setup for full pipeline integration tests.

Downloads S. pombe chromosome III, neighbor proteomes, and RNA-seq reads
from NCBI for use as a realistic end-to-end test of the annotation pipeline.

Test organism: Schizosaccharomyces pombe (fission yeast)
- Genome: chromosome III (NC_003424.3, ~2.5 Mbp)
- Proteins: 10 close neighbor proteomes from NCBI
- RNA-seq: 5 SRA runs (paired-end)

Requires:
- NCBI datasets CLI (https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/)
- SRA Toolkit (https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit)
- Internet access
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import gffutils

from eukan.infra.logging import get_logger

log = get_logger(__name__)

DATA_DIR = Path("tests/data")

GENOME_ACCESSION_FILE = DATA_DIR / "test-genome-accession.txt"
GENOME_ASSEMBLY_ACCESSION_FILE = DATA_DIR / "test-genome-assembly-accession.txt"
NEIGHBOR_ACCESSION_FILE = DATA_DIR / "test-close-neighbor-accessions.txt"
SRA_ACCESSION_FILE = DATA_DIR / "test-SRA-accessions.txt"


# ---------------------------------------------------------------------------
# Prerequisite checks
# ---------------------------------------------------------------------------


def _compress_fastq_parallel(files: list[Path]) -> list[Path]:
    """Compress FASTQ files in parallel using pigz (fast) or gzip (fallback).

    Each file is compressed as a separate subprocess for maximum throughput.
    """
    if not files:
        return []

    compressor = shutil.which("pigz") or shutil.which("gzip")
    if not compressor:
        raise RuntimeError("Neither pigz nor gzip found on PATH")

    log.info("  Compressing %d FASTQ files with %s...", len(files), Path(compressor).name)

    procs = [
        subprocess.Popen([compressor, "--force", str(f)])
        for f in files
    ]
    for p in procs:
        p.wait()

    return [f.with_suffix(".fastq.gz") for f in files]


def _require_tool(name: str) -> str:
    """Check that a tool is on PATH, return its path or raise."""
    path = shutil.which(name)
    if not path:
        raise RuntimeError(
            f"`{name}` not found on PATH. "
            f"Install it from: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/"
            if name == "datasets" else
            "Install it from: https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit"
        )
    return path


def _read_accessions(path: Path) -> list[str]:
    """Read accession list from a file, one per line."""
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


# ---------------------------------------------------------------------------
# Download steps
# ---------------------------------------------------------------------------


def download_genome(output_dir: Path) -> Path:
    """Download the target genome (S. pombe chromosome III) from NCBI."""
    output = output_dir / "genome.fasta"
    if output.exists():
        log.info("Genome already exists: %s", output)
        return output

    accession = _read_accessions(GENOME_ACCESSION_FILE)[0]
    url = (
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=nuccore&id={accession}&rettype=fasta&retmode=text"
    )

    log.info("Downloading genome %s...", accession)
    result = subprocess.run(
        ["curl", "-sLf", "-o", str(output), url],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"Failed to download genome: {result.stderr}")

    # Verify it's a valid FASTA
    text = output.read_text()
    if not text.startswith(">"):
        output.unlink()
        raise RuntimeError(f"Downloaded file is not FASTA: {text[:200]}")

    size = output.stat().st_size
    log.info("Genome downloaded: %s (%s bytes)", output, f"{size:,}")
    return output


def _build_reference_gff3(
    gff_files: list[Path], output: Path, *, target_seqid: str | None = None,
) -> None:
    """Build a reference gene set from NCBI GFF3 annotation files.

    Extracts CDS features and reconstructs a minimal gene > mRNA > CDS
    hierarchy where mRNA and gene bounds span from the initial to terminal
    CDS coordinate.

    If *target_seqid* is given, only genes on that sequence are included.
    """
    with open(output, "w") as out:
        out.write("##gff-version 3\n")

        for gff_file in gff_files:
            try:
                db = gffutils.create_db(
                    str(gff_file), ":memory:", merge_strategy="merge",
                )
            except Exception as e:
                log.warning("Skipping %s: %s", gff_file.name, e)
                continue

            for gene in db.features_of_type("gene", order_by=("seqid", "start")):
                if target_seqid and gene.seqid != target_seqid:
                    continue

                mrnas_with_cds: list[tuple[gffutils.Feature, list[gffutils.Feature]]] = []
                for mrna in db.children(gene, featuretype="mRNA", order_by="start"):
                    cds_list = list(
                        db.children(mrna, featuretype="CDS", order_by="start")
                    )
                    if cds_list:
                        mrnas_with_cds.append((mrna, cds_list))

                if not mrnas_with_cds:
                    continue

                # Gene bounds from all child CDS across all isoforms
                all_cds = [c for _, cl in mrnas_with_cds for c in cl]
                gene_start = min(c.start for c in all_cds)
                gene_end = max(c.end for c in all_cds)

                out.write(
                    f"{gene.seqid}\treference\tgene\t{gene_start}\t{gene_end}\t.\t"
                    f"{gene.strand}\t.\tID={gene.id}\n"
                )

                for mrna, cds_list in mrnas_with_cds:
                    mrna_start = cds_list[0].start
                    mrna_end = cds_list[-1].end

                    out.write(
                        f"{mrna.seqid}\treference\tmRNA\t{mrna_start}\t{mrna_end}"
                        f"\t.\t{mrna.strand}\t.\tID={mrna.id};Parent={gene.id}\n"
                    )

                    for cds in cds_list:
                        out.write(
                            f"{cds.seqid}\treference\tCDS\t{cds.start}\t{cds.end}"
                            f"\t.\t{cds.strand}\t{cds.frame}\t"
                            f"ID={cds.id};Parent={mrna.id}\n"
                        )

    n_genes = sum(1 for line in output.open() if "\tgene\t" in line)
    log.info("Reference GFF3: %s (%d genes)", output, n_genes)


def download_proteomes(output_dir: Path) -> Path:
    """Download neighbor proteomes using NCBI datasets CLI."""
    proteins_output = output_dir / "proteins.faa"
    if proteins_output.exists():
        log.info("Proteomes already exist: %s", proteins_output)
        return proteins_output

    _require_tool("datasets")
    _require_tool("unzip")

    accessions = _read_accessions(NEIGHBOR_ACCESSION_FILE)
    zip_path = output_dir / "neighbor_proteomes.zip"

    log.info("Downloading %d neighbor proteomes...", len(accessions))
    result = subprocess.run(
        [
            "datasets", "download", "genome", "accession",
            "--filename", str(zip_path),
            "--include", "protein",
            *accessions,
        ],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"datasets download failed: {result.stderr}")

    # Extract archive
    extract_dir = output_dir / "_proteome_extract"
    extract_dir.mkdir(exist_ok=True)

    subprocess.run(
        ["unzip", "-o", str(zip_path), "-d", str(extract_dir)],
        capture_output=True, text=True, check=True,
    )

    # Concatenate protein FASTAs
    faa_files = sorted(extract_dir.rglob("*.faa"))
    if not faa_files:
        raise RuntimeError("No .faa files found in downloaded archive")

    with open(proteins_output, "w") as out:
        for faa in faa_files:
            out.write(faa.read_text())

    size = proteins_output.stat().st_size
    n_seqs = sum(1 for line in proteins_output.open() if line.startswith(">"))
    log.info("Proteomes: %s (%s bytes, %d sequences)", proteins_output, f"{size:,}", n_seqs)

    # Cleanup
    shutil.rmtree(extract_dir)
    zip_path.unlink()

    return proteins_output


def download_reference(output_dir: Path) -> Path:
    """Download S. pombe's own GFF3 annotation and build a reference gene set.

    Uses the assembly accession to fetch the GFF3 via NCBI datasets, then
    filters to the target chromosome and extracts a gene > mRNA > CDS hierarchy.
    """
    reference_output = output_dir / "reference.gff3"
    if reference_output.exists():
        log.info("Reference already exists: %s", reference_output)
        return reference_output

    _require_tool("datasets")
    _require_tool("unzip")

    assembly_accession = _read_accessions(GENOME_ASSEMBLY_ACCESSION_FILE)[0]
    target_seqid = _read_accessions(GENOME_ACCESSION_FILE)[0]
    zip_path = output_dir / "reference_gff3.zip"

    log.info("Downloading GFF3 annotation for %s...", assembly_accession)
    result = subprocess.run(
        [
            "datasets", "download", "genome", "accession",
            assembly_accession,
            "--filename", str(zip_path),
            "--include", "gff3",
        ],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"datasets download failed: {result.stderr}")

    extract_dir = output_dir / "_reference_extract"
    extract_dir.mkdir(exist_ok=True)

    subprocess.run(
        ["unzip", "-o", str(zip_path), "-d", str(extract_dir)],
        capture_output=True, text=True, check=True,
    )

    gff_files = sorted(extract_dir.rglob("*.gff"))
    if not gff_files:
        raise RuntimeError("No GFF3 files found in downloaded archive")

    _build_reference_gff3(gff_files, reference_output, target_seqid=target_seqid)

    # Cleanup
    shutil.rmtree(extract_dir)
    zip_path.unlink()

    return reference_output


def download_sra_reads(output_dir: Path) -> list[Path]:
    """Download RNA-seq reads using SRA Toolkit (prefetch + fasterq-dump).

    Produces gzipped FASTQ files. Continues on individual run failures.
    """
    _require_tool("prefetch")
    _require_tool("fasterq-dump")

    accessions = _read_accessions(SRA_ACCESSION_FILE)
    fastq_files: list[Path] = []
    failed: list[str] = []

    for acc in accessions:
        # Check if already downloaded (gz or uncompressed)
        existing = list(output_dir.glob(f"{acc}*.fastq.gz"))
        if existing:
            log.info("SRA %s already downloaded: %s", acc, [f.name for f in existing])
            fastq_files.extend(existing)
            continue

        log.info("Downloading SRA run %s...", acc)
        try:
            # prefetch to local cache
            result = subprocess.run(
                ["prefetch", acc, "-O", str(output_dir)],
                capture_output=True, text=True,
            )
            if result.returncode != 0:
                log.warning("prefetch %s failed: %s", acc, result.stderr.strip())

            # Convert to FASTQ (fasterq-dump can work from prefetched or remote)
            # Use the prefetched .sra file if it exists, otherwise let it fetch
            sra_file = output_dir / acc / f"{acc}.sra"
            source = str(sra_file) if sra_file.exists() else acc

            subprocess.run(
                [
                    "fasterq-dump", source,
                    "--outdir", str(output_dir),
                    "--split-3",
                    "--threads", "4",
                ],
                capture_output=True, text=True, check=True,
            )

            # Clean up SRA cache before compression
            sra_cache = output_dir / acc
            if sra_cache.is_dir():
                shutil.rmtree(sra_cache)

            # Compress FASTQ files in parallel (one gzip per file)
            new_files = list(output_dir.glob(f"{acc}*.fastq"))
            if new_files:
                _compress_fastq_parallel(new_files)

            gz_files = list(output_dir.glob(f"{acc}*.fastq.gz"))
            fastq_files.extend(gz_files)
            log.info("  -> %s", [f.name for f in gz_files])

        except subprocess.CalledProcessError as e:
            log.warning("Failed to download %s: %s", acc, e.stderr or str(e))
            failed.append(acc)
            continue

    if failed:
        log.warning("%d SRA runs failed: %s", len(failed), failed)

    log.info("Total FASTQ files: %d", len(fastq_files))
    return fastq_files


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------


def setup_test_data(output_dir: Path) -> dict[str, Path | list[Path]]:
    """Download all test data. Returns paths to downloaded files."""
    output_dir.mkdir(parents=True, exist_ok=True)

    results: dict[str, Path | list[Path]] = {}

    log.info("Setting up test data in %s", output_dir)
    log.info("=" * 60)

    log.info("\n[1/4] Genome")
    results["genome"] = download_genome(output_dir)

    log.info("\n[2/4] Neighbor proteomes")
    results["proteins"] = download_proteomes(output_dir)

    log.info("\n[3/4] Reference gene set")
    results["reference"] = download_reference(output_dir)

    log.info("\n[4/4] RNA-seq reads (SRA)")
    results["reads"] = download_sra_reads(output_dir)

    log.info("\n" + "=" * 60)
    log.info("Test data setup complete.")
    log.info("  Genome:    %s", results["genome"])
    log.info("  Proteins:  %s", results["proteins"])
    log.info("  Reference: %s", results["reference"])
    log.info("  Reads:     %d FASTQ files", len(results["reads"]))

    return results


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def validate_test_data(output_dir: Path) -> list[tuple[str, str, bool]]:
    """Check that test data is present and looks valid.

    Returns list of (name, message, ok) tuples.
    """
    results: list[tuple[str, str, bool]] = []

    # Genome
    genome = output_dir / "genome.fasta"
    if genome.exists():
        size = genome.stat().st_size
        if size > 1_000_000:  # S. pombe chrIII is ~2.5 Mbp
            results.append(("genome", f"genome.fasta ({size:,} bytes)", True))
        else:
            results.append(("genome", f"genome.fasta too small ({size:,} bytes)", False))
    else:
        results.append(("genome", "genome.fasta not found", False))

    # Proteins
    proteins = output_dir / "proteins.faa"
    if proteins.exists():
        n_seqs = sum(1 for line in proteins.open() if line.startswith(">"))
        if n_seqs > 1000:
            results.append(("proteins", f"proteins.faa ({n_seqs:,} sequences)", True))
        else:
            results.append(("proteins", f"proteins.faa only {n_seqs} sequences (expected >1000)", False))
    else:
        results.append(("proteins", "proteins.faa not found", False))

    # Reference GFF3
    reference = output_dir / "reference.gff3"
    if reference.exists():
        n_genes = sum(1 for line in reference.open() if "\tgene\t" in line)
        if n_genes > 100:
            results.append(("reference", f"reference.gff3 ({n_genes:,} genes)", True))
        else:
            results.append(("reference", f"reference.gff3 only {n_genes} genes (expected >100)", False))
    else:
        results.append(("reference", "reference.gff3 not found", False))

    # Reads — expect paired-end _1 and _2 files per accession
    accessions = _read_accessions(SRA_ACCESSION_FILE)
    fwd_count = len(list(output_dir.glob("SRR*_1.fastq.gz")))
    rev_count = len(list(output_dir.glob("SRR*_2.fastq.gz")))

    if fwd_count >= len(accessions) and rev_count >= len(accessions):
        results.append(("reads", f"{fwd_count} paired-end read pairs (.fastq.gz)", True))
    elif fwd_count > 0 or rev_count > 0:
        results.append(("reads", f"{fwd_count} forward + {rev_count} reverse files (expected {len(accessions)} pairs)", False))
    else:
        # Check for old _1/_3 naming from --split-files
        old_files = len(list(output_dir.glob("SRR*_3.fastq.gz")))
        if old_files > 0:
            results.append(("reads", "found _1/_3 files (need _1/_2 pairs — delete and re-download)", False))
        else:
            results.append(("reads", "no FASTQ files found (run: eukan dev setup-test-data)", False))

    return results
