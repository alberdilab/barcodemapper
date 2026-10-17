#!/usr/bin/env bash
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00

set -euo pipefail
: "${BENCH_ROOT:?Export the benchmark root before submission}"

inputs="$BENCH_ROOT/projects/01-mock-communities/inputs"
python3 "$BENCH_ROOT/plans/make_pilot.py" "$inputs" SRR14351851 --pairs 100000
