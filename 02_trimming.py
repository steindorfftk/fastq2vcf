import os
import argparse
import subprocess
import time
from pathlib import Path

def run_command(command, verbose):
    if verbose:
        print(f"Running: {command}")
    subprocess.run(command, shell=True, check=True)

def main():
    parser = argparse.ArgumentParser(description="Trim bases from sequencing reads using Trimmomatic.")
    parser.add_argument("-n", "--num_bases", type=int, required=True, help="Number of bases to remove.")
    parser.add_argument("-p", "--position", choices=["start", "end"], required=True, help="Remove bases from the start or end of the reads.")
    parser.add_argument("-e", "--ending", choices=["se", "pe"], required=True, help="Specify if the reads are single or paired-ended.")
    parser.add_argument("-nv", "--verbose", action="store_false", help="Disable verbose output.")    
    parser.add_argument("-max", "--max_length", type=int, help="Maximum length to crop.")  
    parser.add_argument("-th", "--threads", help="Number of threads to use (default: 6)", type=int, default=6)
    parser.add_argument("-min", "--minlen", type=int, help="Minimum length of the reads after trimming (optional).")
    args = parser.parse_args()

    trimmomatic_path = "temporary/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar"
    
    # Update to use TAILCROP for the end position.
    trimmomatic_option = "HEADCROP" if args.position == "start" else "CROP"

    # Define crop target if position is 'end'
    if args.position == 'end':
         args.num_bases = args.max_length - args.num_bases

    srr_list = []
    with open('input/SRRAccList.txt', 'r') as text:
        for line in text:
            srr_list.append(line.strip())

    output_dir = "temporary/trimmomatic/output"
    os.makedirs(output_dir, exist_ok=True)
    
    # Create the FastQC output directory
    fastqc_output = "temporary/fastqc/post_trimming_fastqc"
    os.makedirs(fastqc_output, exist_ok=True)

    for srr in srr_list:
        input_prefix = f'temporary/fastq_dump/{srr}'
        output_prefix = f'{output_dir}/{srr}'

        # Start timer for each SRR
        start_time = time.time()

        # Construct the base Trimmomatic command
        trimmomatic_command = (
            f'java -jar {trimmomatic_path} {args.ending.upper()} '
            f'{input_prefix}_1.fastq {input_prefix}_2.fastq '
            f'{output_prefix}_1.fastq {output_prefix}_1_unpaired.fastq '
            f'{output_prefix}_2.fastq {output_prefix}_2_unpaired.fastq '
            f'{trimmomatic_option}:{args.num_bases} '
            f'-threads {args.threads}'
        ) if args.ending == "pe" else (
            f'java -jar {trimmomatic_path} SE '
            f'{input_prefix}.fastq '
            f'{output_prefix}.fastq '
            f'{trimmomatic_option}:{args.num_bases} '
            f'-threads {args.threads}'
        )

        # Add MINLEN argument if it's provided
        if args.minlen:
            trimmomatic_command += f' MINLEN:{args.minlen}'

        # Execute the Trimmomatic command
        run_command(trimmomatic_command, args.verbose)

        # End timer for SRR
        end_time = time.time()
        elapsed_time = end_time - start_time

        if args.verbose:
            print(f"Time taken for {srr}: {elapsed_time:.2f} seconds")

        # Now run FastQC on the trimmed fastq files
        fastq_files = list(Path(output_dir).glob(f"{srr}*.fastq"))
        for fastq_file in fastq_files:
            if fastq_file.stat().st_size == 0:  # Check if the file is empty (0 bytes)
                if args.verbose:
                    print(f"File {fastq_file.name} is empty. Deleting it.")
                fastq_file.unlink()  # Delete the empty file
            else:
                if args.verbose:
                    print(f"Running FastQC on {fastq_file.name}...")
                # Use the same threads for FastQC as for Trimmomatic
                run_command(f"fastqc -t {args.threads} -o {fastqc_output} {fastq_file}", args.verbose)

    if args.verbose:
        print("FastQC analysis completed.")

if __name__ == "__main__":
    main()

