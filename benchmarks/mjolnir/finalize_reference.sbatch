#!/usr/bin/env bash
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00

set -euo pipefail
: "${BENCH_ROOT:?Export the benchmark root before submission}"
: "${REFERENCE:?Export the reference basename before submission}"

python3 "$BENCH_ROOT/plans/finalize_reference.py" \
  "$BENCH_ROOT/references/$REFERENCE.fa.gz"
