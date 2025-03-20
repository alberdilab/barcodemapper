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
# dietscan installation path
#####

PACKAGE_DIR = Path(__file__).parent
CONFIG_PATH = PACKAGE_DIR / "workflow" / "config.yaml"

def load_config():
    """Load fixed variables from config.yaml."""
    if CONFIG_PATH.exists():
        with open(CONFIG_PATH, "r") as f:
            return yaml.safe_load(f)
    return {}

config_vars = load_config()

#####
# Function definitions
#####

def unlock_snakemake(tmpdir, profile):
    unlock_command = [
        "/bin/bash", "-c",  # Ensures the module system works properly
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {output_dir} "
        f"--configfile {CONFIG_PATH} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--unlock"
    ]

    subprocess.run(unlock_command, shell=False, check=True)
    print(f"The output directory {output_dir} has been succesfully unlocked.")

def run_snakemake_origin(tmp_dir, outputfile, bold_db, unite_db, bold_retain, unite_retain, profile):
    snakemake_command = [
        "/bin/bash", "-c",
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {tmp_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} bold_db={bold_db} unite_db={unite_db} bold_retain={bold_retain} unite_retain={unite_retain} tmp_dir={tmp_dir} output_file={output_file}"
    ]
    subprocess.run(snakemake_command, shell=False, check=True)

def run_snakemake_database(tmp_dir, output_file, profile):
    snakemake_command = [
        "/bin/bash", "-c",
        "snakemake "
        f"-s {PACKAGE_DIR / 'workflow' / 'Snakefile'} "
        f"--directory {tmp_dir} "
        f"--workflow-profile {PACKAGE_DIR / 'profile' / profile} "
        f"--configfile {CONFIG_PATH} "
        f"--config package_dir={PACKAGE_DIR} tmp_dir={tmp_dir} output_file={output_file}"
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
# dietscan execution
#####

def main():
    parser = argparse.ArgumentParser(description="Dietary profiling from metagenomic data.")
    parser.add_argument("-i", "--input", type=str, required=False, help="Path to the input directory containing the sequencing reads.")
    parser.add_argument("-1", "--read1", type=str, required=False, help="Comma-separated list of forward reads.")
    parser.add_argument("-2", "--read2", type=str, required=False, help="Comma-separated list of reverse reads.")
    parser.add_argument("-o", "--output", type=str, required=False, help="Cryosection identifier (required).")
    parser.add_argument("-d", "--database", type=str, required=False, help="Combined database (fasta).")
    parser.add_argument("-b", "--bold", type=str, required=False, help="Bold database (fasta).")
    parser.add_argument("-u", "--unite", type=str, required=False, help="Unite database (fasta).")
    parser.add_argument("-x", "--bold_retain", type=str, required=False, default="k__Animalia", help="Comma-separated list of taxa to consider in the BOLD database (e.g. 'o__Coleoptera,o__Lepidoptera')")
    parser.add_argument("-y", "--unite_retain", type=str, required=False, default="k__Viridiplantae,p__Basidiomycota", help="Comma-separated list of taxa to consider in the UNITE database (e.g. 'k__Fungi,k__Viridioplantae')")
    parser.add_argument("-t", "--tmpdir", type=str, required=False, help="Directory where the temporary files are stored")
    parser.add_argument("-s", "--slurm", action="store_true", required=False, help="Whether to use slurm")

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputdir = args.input
    read1 = args.read1
    read2 = args.read2
    outputfile = args.output
    database = args.database
    bold_db = args.bold
    unite_db = args.unite
    bold_retain = args.bold_retain
    unite_retain = args.unite_retain
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

    if args.database and (args.bold or args.unite):
        print(f"    Please, provide either the combined database or the original bold and unite databases.")
        return

    #####
    # tmp directory
    #####

    tmp_folder= "dietscan" + datetime.now().strftime("%Y%m%d%H%M")
    if not args.tmpdir:
        tmp_dir = tmp_folder
    else:
        tmp_dir = os.path.join(tmpdir, dir_name)

    os.makedirs(tmp_dir, exist_ok=True)

    #####
    # run
    #####

    if args.read1 and args.read2:
        inputlist_to_samples(read1, read2, tmp_dir)

    if args.input:
        inputdir_to_samples(inputdir, tmp_dir)

    if args.database:
        dest_file = os.path.join(tmp_dir, "database/dietscan_db.fa")
        os.makedirs(os.path.join(tmp_dir, "database"), exist_ok=True)
        shutil.copy(database, dest_file)
        run_snakemake_database(tmp_dir, outputfile, profile)
    else:
        run_snakemake_origin(tmp_dir, outputfile, bold_db, unite_db, bold_retain, unite_retain, profile)
