#!/usr/bin/env bash
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=04:00:00
#SBATCH --job-name=cmp-nmis-test
# Tiny test: 300 reads of ERR6001161 R1, vsearch with --n_mismatch (same settings as the frozen
# config), then the exact acceptance filter and the N-run diagnostic. Writes only under work/.../test-nmismatch.
set -euo pipefail
: "${BENCH_ROOT:?}"
export PROJECT_ENV="$BENCH_ROOT/projects/03-method-comparison/software/env-1.2.3"
exec bash -ic 'useconda; conda activate "$PROJECT_ENV"; exec bash "$BENCH_ROOT/tasks/1.2.3-comparators/test_nmismatch.sh"'
