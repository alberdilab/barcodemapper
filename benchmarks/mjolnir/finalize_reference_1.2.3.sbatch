#!/usr/bin/env bash
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
# Freeze a 1.2.3-built reference in references-1.2.3/ (same finalize_reference.py as 1.1.2).
set -euo pipefail
: "${BENCH_ROOT:?}"
: "${REFERENCE:?}"
python3 "$BENCH_ROOT/plans/finalize_reference.py" \
  "$BENCH_ROOT/references-1.2.3/$REFERENCE.fa.gz"
