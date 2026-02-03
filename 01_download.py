import argparse
import subprocess
import os
from pathlib import Path
from multiprocessing import Pool

# Usage: python3 01_download.py -o Mus/Homo -e pe -m

def run_command(command, verbose=True):
    """Run shell command and handle errors."""
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if verbose:
        print(result.stdout)
    if result.returncode != 0:
        print(f"Error: {result.stderr}")

# Argument Parsing
parser = argparse.ArgumentParser(description="Download SRA files, convert to FASTQ, and run initial FastQC.")
parser.add_argument('-nv', '--verbose', action='store_false', help='Disable verbose output')
parser.add_argument('-m', '--memory', action='store_true', help='Delete SRA files after conversion to save memory')
parser.add_argument('-o', '--organism', type=str, required=True, choices=['Mus', 'Homo'], help='Organism name')
parser.add_argument('-e', '--end', type=str, default='se', choices=['se', 'pe'], help='Single-end (se) or paired-end (pe) sequencing')

args = parser.parse_args()

# Set paths
genome_paths = {
    'Mus': ('temporary/bowtie2/refGenomes/mm10.fasta', 'temporary/feature_counts/data/mm10_RefSeq_exon.txt'),
    'Homo': ('temporary/bowtie2/refGenomes/GRCh38_noalt_as', 'temporary/feature_counts/data/hg38_RefSeq_exon.txt')
}

refGenomePath, annotationPath = genome_paths[args.organism]

# Read Accession List
accession_list = [line.strip() for line in Path('input/SRRAccList.txt').read_text().splitlines()]

# Download SRA Data
if args.verbose:
    print("Downloading SRA data...")

run_command('prefetch --option-file input/SRRAccList.txt', args.verbose)

# Convert SRA to FASTQ
def convert_sra(acc):
    """Convert SRA to FASTQ using fasterq-dump."""
    output_dir = Path('temporary/fastq_dump')
    output_dir.mkdir(parents=True, exist_ok=True)

    command = f"fasterq-dump {acc} -O {output_dir}"
    run_command(command, args.verbose)

    expected_files = [f"{acc}.fastq"] if args.end == 'se' else [f"{acc}_1.fastq", f"{acc}_2.fastq"]
    if all((output_dir / f).exists() for f in expected_files):
        if args.verbose:
            print(f"Conversion complete for {acc}")
    else:
        print(f"Error: Conversion failed for {acc}")

# Use multiprocessing for parallel processing
with Pool(processes=4) as pool:
    pool.map(convert_sra, accession_list)

# Delete SRA Files if Low Memory Mode is Enabled
if args.memory:
    sra_dir = Path('temporary/sratoolkit/sra')
    for file in sra_dir.glob("*.sra"):
        file.unlink()
    if args.verbose:
        print("Deleted SRA files to save memory.")

# Run FastQC
fastqc_output = Path('temporary/fastqc/download_fastqc')
fastqc_output.mkdir(parents=True, exist_ok=True)

fastq_files = list(Path('temporary/fastq_dump').glob("*.fastq"))
for fastq_file in fastq_files:
    if args.verbose:
        print(f"Running FastQC on {fastq_file.name}...")
    run_command(f"fastqc -o {fastqc_output} {fastq_file}", args.verbose)

if args.verbose:
    print("FastQC analysis completed.")

