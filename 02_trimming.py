import os
import argparse
import subprocess
import time
from pathlib import Path

def run_command(command, verbose):
    if verbose:
        print(f"Running: {command}")
    # Using shell=True to allow the conda wrapper to be found in the environment path
    subprocess.run(command, shell=True, check=True)

def main():
    parser = argparse.ArgumentParser(description="Trim bases from sequencing reads using Trimmomatic.")
    parser.add_argument("-n", "--num_bases", type=int, required=True, help="Number of bases to remove.")
    parser.add_argument("-p", "--position", choices=["start", "end"], required=True, help="Remove bases from the start or end of the reads.")
    parser.add_argument("-e", "--ending", choices=["se", "pe"], required=True, help="Specify if the reads are single or paired-ended.")
    parser.add_argument("-nv", "--verbose", action="store_true", help="Enable verbose output.")    
    parser.add_argument("-max", "--max_length", type=int, help="Maximum length to crop.")  
    parser.add_argument("-th", "--threads", help="Number of threads to use (default: 6)", type=int, default=6)
    parser.add_argument("-min", "--minlen", type=int, help="Minimum length of the reads after trimming (optional).")
    args = parser.parse_args()

    # In Conda, 'trimmomatic' is a direct command. No .jar path needed.
    trimmomatic_option = "HEADCROP" if args.position == "start" else "CROP"

    # Define crop target if position is 'end'
    if args.position == 'end':
        if args.max_length is None:
            print("Error: --max_length is required when position is 'end'")
            return
        # CROP:X keeps the first X bases. To remove N bases from end, we keep (Total - N)
        args.num_bases = args.max_length - args.num_bases

    srr_list = []
    # Ensure this path exists or matches your GitHub architecture
    input_acc_list = 'input/SRRAccList.txt'
    if not os.path.exists(input_acc_list):
        print(f"Error: {input_acc_list} not found.")
        return

    with open(input_acc_list, 'r') as text:
        for line in text:
            if line.strip():
                srr_list.append(line.strip())

    output_dir = "temporary/trimmomatic/output"
    os.makedirs(output_dir, exist_ok=True)
    
    fastqc_output = "temporary/fastqc/post_trimming_fastqc"
    os.makedirs(fastqc_output, exist_ok=True)

    for srr in srr_list:
        input_prefix = f'temporary/fastq_dump/{srr}'
        output_prefix = f'{output_dir}/{srr}'

        start_time = time.time()

        # Conda version uses lowercase 'trimmomatic' and handles java internally
        if args.ending == "pe":
            trimmomatic_command = (
                f'trimmomatic PE -threads {args.threads} '
                f'{input_prefix}_1.fastq {input_prefix}_2.fastq '
                f'{output_prefix}_1.fastq {output_prefix}_1_unpaired.fastq '
                f'{output_prefix}_2.fastq {output_prefix}_2_unpaired.fastq '
                f'{trimmomatic_option}:{args.num_bases}'
            )
        else:
            trimmomatic_command = (
                f'trimmomatic SE -threads {args.threads} '
                f'{input_prefix}.fastq '
                f'{output_prefix}.fastq '
                f'{trimmomatic_option}:{args.num_bases}'
            )

        if args.minlen:
            trimmomatic_command += f' MINLEN:{args.minlen}'

        try:
            run_command(trimmomatic_command, args.verbose)
        except subprocess.CalledProcessError as e:
            print(f"Error processing {srr}: {e}")
            continue

        elapsed_time = time.time() - start_time
        if args.verbose:
            print(f"Time taken for {srr}: {elapsed_time:.2f} seconds")

        # Post-trimming FastQC logic
        fastq_files = list(Path(output_dir).glob(f"{srr}*.fastq"))
        for fastq_file in fastq_files:
            if fastq_file.stat().st_size == 0:
                if args.verbose:
                    print(f"File {fastq_file.name} is empty. Deleting.")
                fastq_file.unlink()
            else:
                if args.verbose:
                    print(f"Running FastQC on {fastq_file.name}...")
                run_command(f"fastqc -t {args.threads} -o {fastqc_output} {fastq_file}", args.verbose)

    if args.verbose:
        print("Pipeline step completed.")

if __name__ == "__main__":
    main()
