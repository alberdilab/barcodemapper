#!/usr/bin/env bash
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00

set -euo pipefail
: "${BENCH_ROOT:?}"
: "${BENCH_ENV:?}"
: "${RUN:?}"

"$BENCH_ENV/bin/python" "$BENCH_ROOT/tasks/1.2.3-pilots/best_hit_ablation_1.2.3.py" "$BENCH_ROOT" "$RUN"
