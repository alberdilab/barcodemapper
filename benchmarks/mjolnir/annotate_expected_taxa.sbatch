#!/usr/bin/env bash
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00

set -euo pipefail
: "${BENCH_ROOT:?Export the benchmark root before submission}"

python3 "$BENCH_ROOT/plans/annotate_expected_taxa.py" "$BENCH_ROOT"
