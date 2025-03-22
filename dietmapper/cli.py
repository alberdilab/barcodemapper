import argparse
import os
import sys
import subprocess
import yaml
import re
import json
import requests
import pandas as pd
from pathlib import Path
import shutil
from datetime import datetime
from collections import defaultdict

#####
# dietmapper installation path
#####

PACKAGE_DIR = Path(__file__).parent

#####
# Function definitions
#####

def unlock_snakemake(tmp_dir, config_path, profile):
    unlock_command = [
        "/bin/bash", "-c",  # Ensures the module system works properly
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {tmp_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {config_path} "
        f"--config package_dir={PACKAGE_DIR} "
        f"--unlock"
    ]

    subprocess.run(unlock_command, shell=False, check=True)
    print(f"The output directory {output_dir} has been succesfully unlocked.")

def run_snakemake(tmp_dir, config_path, profile):
    snakemake_command = [
        "/bin/bash", "-c",
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {tmp_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {config_path} "
        f"--config package_dir={PACKAGE_DIR}"
    ]
    subprocess.run(snakemake_command, shell=False, check=True)

def inputdir_to_samples(input_dir, output_dir):
    # Define the directory containing the raw reads
    READS_DIR = Path(input_dir).resolve()

    # Initialize dictionaries
    SAMPLE_TO_READS1 = defaultdict(list)
    SAMPLE_TO_READS2 = defaultdict(list)

    # Regular expression to capture sample names
    pattern = re.compile(r"^(.*)_\d\.fq\.gz$")  # Captures everything before "_1.fq.gz" or "_2.fq.gz"

    # Scan the directory
    for filename in os.listdir(READS_DIR):
        if filename.endswith(".fq.gz"):
            full_path = os.path.join(READS_DIR, filename)

            # Extract sample name using regex
            match = pattern.match(filename)
            if match:
                sample_name = match.group(1)  # Everything before _1.fq.gz or _2.fq.gz

                # Sort into forward and reverse reads
                if "_1.fq.gz" in filename:
                    SAMPLE_TO_READS1[sample_name].append(full_path)
                elif "_2.fq.gz" in filename:
                    SAMPLE_TO_READS2[sample_name].append(full_path)

    # Convert defaultdict to standard dict (optional)
    SAMPLE_TO_READS1 = dict(SAMPLE_TO_READS1)
    SAMPLE_TO_READS2 = dict(SAMPLE_TO_READS2)

    os.makedirs(f"{output_dir}/data", exist_ok=True)
    with open(f"{output_dir}/data/sample_to_reads1.json", "w") as f:
        json.dump(SAMPLE_TO_READS1, f, indent=4)

    with open(f"{output_dir}/data/sample_to_reads2.json", "w") as f:
        json.dump(SAMPLE_TO_READS2, f, indent=4)

def inputlist_to_samples(read1, read2, output_dir):
    forward_files = read1.split(",")
    reverse_files = read2.split(",")

    # Initialize dictionaries
    SAMPLE_TO_READS1 = defaultdict(list)
    SAMPLE_TO_READS2 = defaultdict(list)

    # Regular expression to capture sample names
    pattern = re.compile(r"^(.*)_\d\.fq\.gz$")  # Captures everything before "_1.fq.gz" or "_2.fq.gz"

    # Process forward reads
    for file in forward_files:
        filename = os.path.basename(file)
        match = pattern.match(filename)
        if match:
            sample_name = match.group(1)
            SAMPLE_TO_READS1[sample_name].append(file)

    # Process reverse reads
    for file in reverse_files:
        filename = os.path.basename(file)
        match = pattern.match(filename)
        if match:
            sample_name = match.group(1)
            SAMPLE_TO_READS2[sample_name].append(file)

    # Convert defaultdict to standard dict
    SAMPLE_TO_READS1 = dict(SAMPLE_TO_READS1)
    SAMPLE_TO_READS2 = dict(SAMPLE_TO_READS2)

    # Ensure output directory exists
    os.makedirs(f"{output_dir}/data", exist_ok=True)

    # Save JSON files
    with open(f"{output_dir}/data/sample_to_reads1.json", "w") as f:
        json.dump(SAMPLE_TO_READS1, f, indent=4)

    with open(f"{output_dir}/data/sample_to_reads2.json", "w") as f:
        json.dump(SAMPLE_TO_READS2, f, indent=4)

#####
# dietmapper execution
#####

def main():
    parser = argparse.ArgumentParser(description="DietMapper: dietary profiling from metagenomic data.")
    parser.add_argument("-i", "--input", type=str, required=False, help="Path to the input directory containing the sequencing reads.")
    parser.add_argument("-1", "--read1", type=str, required=False, help="Comma-separated list of forward reads.")
    parser.add_argument("-2", "--read2", type=str, required=False, help="Comma-separated list of reverse reads.")
    parser.add_argument("-o", "--output", type=str, required=False, help="Output taxonomy file.")
    parser.add_argument("-p", "--plot", action="store_true", required=False, help="Output sunburst plot html.")
    parser.add_argument("-d", "--database", type=str, required=True, help="Combined DietMapper database (fasta).")
    parser.add_argument("-b", "--bold", type=str, required=False, help="Bold database (fasta).")
    parser.add_argument("-u", "--unite", type=str, required=False, help="Unite database (fasta).")
    parser.add_argument("-x", "--bold_retain", type=str, required=False, default="k__Animalia", help="Comma-separated list of taxa to consider in the BOLD database (e.g. 'o__Coleoptera,o__Lepidoptera')")
    parser.add_argument("-y", "--unite_retain", type=str, required=False, default="k__Viridiplantae,p__Basidiomycota", help="Comma-separated list of taxa to consider in the UNITE database (e.g. 'k__Fungi,k__Viridioplantae')")
    parser.add_argument("-t", "--tmpdir", type=str, required=False, help="Directory where the temporary files are stored")
    parser.add_argument("-s", "--slurm", action="store_true", required=False, help="Whether to use slurm")
    parser.add_argument("--unlock", action="store_true", required=False, help="Whether to unlock the directory")
    parser.add_argument("--build", action="store_true", required=False, help="Only build the database")

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    #####
    # select profile
    #####

    if args.slurm:
        profile = 'slurm'
    else:
        profile = 'local'

    #####
    # input check
    #####

    if args.input and (args.read1 or args.read2):
        print(f"    Please, provide either an input directory or a list of input files.")
        return

    if args.read1 and not args.read2:
        print(f"    Please, provide both forward and reverse sequencing read files.")
        return

    if args.read2 and not args.read1:
        print(f"    Please, provide both forward and reverse sequencing read files.")
        return

    if not args.build and not (args.input or (args.read1 and args.read2)):
        print(f"    No input data was specified. Use either the -i, or the -1 and -2 arguments.")
        return

    if args.database and not (args.bold and args.unite):
        dietmapper_db_path = Path(args.database)
        if not dietmapper_db_path.exists():
            print(f"Error: Database file {dietmapper_db_path} does not exist.")
            sys.exit(1)
        else:
            if dietmapper_db_path.suffix.lower() not in [".fa", ".fasta"]:
                print(f"Error: Database file {dietmapper_db_path} does not have a correct (.fa or .fasta) extension.")
                sys.exit(1)

    if args.bold:
        bold_db_path = Path(args.bold)
        if not bold_db_path.exists():
            print(f"Error: BOLD database file {bold_db_path} does not exist.")
            sys.exit(1)

    if args.unite:
        unite_db_path = Path(args.unite)
        if not unite_db_path.exists():
            print(f"Error: UNITE database file {unite_db_path} does not exist.")
            sys.exit(1)

    if args.database and args.bold and args.unite:
        dietmapper_db_path = Path(args.database)
        if dietmapper_db_path.exists():
            print(f"Warning: The DietMapper database {dietmapper_db_path} already exists.")
            print(f"If you want to re-build the DietMapper database first remove this one.")

    #####
    # depenency check
    #####

    def check_tool_installed(tool_name):
        if shutil.which(tool_name) is None:
            print(f"    Error: {tool_name} is not installed or not in your PATH.", file=sys.stderr)
            return False
        return True

    tools = ["bowtie2", "samtools", "fastp"]
    missing = [tool for tool in tools if not check_tool_installed(tool)]

    if missing:
        sys.exit(f"Make sure {', '.join(missing)} are installed before running DietMapper. Optionally, create the DietMapper conda environment containing all dependencies following the explanations in the DietMapper documentation.")

    #####
    # tmp directory
    #####

    if not args.tmpdir:
        tmp_dir = "dietmapper_" + datetime.now().strftime("%Y%m%d%H%M")
    else:
        tmp_dir = args.tmpdir

    os.makedirs(Path(tmp_dir).resolve(), exist_ok=True)

    #####
    # database build
    #####

    if args.build:
        buildonly="yes"

        if not args.bold:
            print(f"Error: You did not specify a BOLD database.")
            sys.exit(1)
        else:
            bold_db_path = Path(args.bold)
            if not bold_db_path.exists():
                print(f"Error: BOLD database file {bold_db_path} does not exist.")
                sys.exit(1)

        if not args.unite:
            print(f"Error: You did not specify a UNITE database.")
        else:
            unite_db_path = Path(args.unite)
            if not unite_db_path.exists():
                print(f"Error: UNITE database file {unite_db_path} does not exist.")
                sys.exit(1)

        if not args.output:
            output_file="no"

    else:
        buildonly="no"
        output_file = str(Path(args.output).resolve())

    #####
    # plot
    #####

    if args.plot:
        output_path=Path(args.output).resolve()
        sunburst_file=str(output_path.with_name(output_path.stem + '.html'))
    else:
        sunburst_file="no"

    #####
    # config
    #####

    config_dir = Path(tmp_dir) / "data"
    config_dir.mkdir(parents=True, exist_ok=True)
    config_path = config_dir / "config.yaml"
    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    if args.database and not (args.bold and args.unite):
        config_data = {
            "run": current_time,
            "dietmapper_db": str(Path(args.database).resolve()),
            "tmp_dir": str(Path(tmp_dir).resolve()),
            "output_file": output_file,
            "build_only": buildonly,
            "sunburst_file": sunburst_file
        }
    else:
        config_data = {
            "run": current_time,
            "dietmapper_db": str(Path(args.database).resolve()),
            "bold_db": str(Path(args.bold).resolve()),
            "unite_db": str(Path(args.unite).resolve()),
            "bold_retain": args.bold_retain,
            "unite_retain": args.unite_retain,
            "tmp_dir": str(Path(tmp_dir).resolve()),
            "output_file": output_file,
            "build_only": buildonly,
            "sunburst_file": sunburst_file
        }

    with open(config_path, "w") as f:
        yaml.dump(config_data, f)

    #####
    # run
    #####

    if args.read1 and args.read2:
        inputlist_to_samples(args.read1, args.read2, Path(tmp_dir).resolve())

    if args.input:
        inputdir_to_samples(args.input, Path(tmp_dir).resolve())

    if args.unlock:
        unlock_snakemake(Path(tmp_dir).resolve(), config_path, profile)
    else:
        run_snakemake(Path(tmp_dir).resolve(), config_path, profile)
