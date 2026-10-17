#!/usr/bin/env bash
#SBATCH --partition=filetransfer
#SBATCH --qos=filetransfer
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=1-00:00:00

set -euo pipefail
: "${BENCH_ROOT:?Export the benchmark root before submission}"
: "${STUDY:?Export the study accession before submission}"
: "${RUN:?Export the run accession before submission}"

inputs="$BENCH_ROOT/projects/01-mock-communities/inputs"
python3 "$BENCH_ROOT/plans/download_ena_run.py" \
  "$inputs/ena-$STUDY.tsv" "$RUN" "$inputs"
