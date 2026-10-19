#!/usr/bin/env bash
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
# Env: BENCH_ROOT. Converts and writes manifests for all 6 (setting, run) combinations.
set -euo pipefail
: "${BENCH_ROOT:?}"
export PROJECT_ENV="$BENCH_ROOT/projects/03-method-comparison/software/env-1.2.3"
exec bash -ic 'useconda; conda activate "$PROJECT_ENV"; exec bash "$BENCH_ROOT/tasks/1.2.3-kraken2/finish_kraken2.sh"'
