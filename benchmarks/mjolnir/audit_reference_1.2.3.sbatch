#!/usr/bin/env bash
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
# Audit a frozen 1.2.3 reference (same audit_reference_taxa.py as 1.1.2).
set -euo pipefail
: "${BENCH_ROOT:?}"
: "${REFERENCE:?}"
python3 "$BENCH_ROOT/plans/audit_reference_taxa.py" \
  "$BENCH_ROOT/references-1.2.3/$REFERENCE.fa.gz" \
  "$BENCH_ROOT/references-1.2.3/audits/$REFERENCE"
