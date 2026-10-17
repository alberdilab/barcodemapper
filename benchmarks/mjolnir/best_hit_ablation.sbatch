#!/usr/bin/env bash
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00

set -euo pipefail
: "${BENCH_ROOT:?Export the benchmark root before submission}"
: "${BENCH_ENV:?Export the Project 03 Conda prefix before submission}"
: "${RUN:?Export the pilot ENA run accession before submission}"

"$BENCH_ENV/bin/python" "$BENCH_ROOT/plans/best_hit_ablation.py" "$BENCH_ROOT" "$RUN"
