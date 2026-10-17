#!/usr/bin/env bash
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00

set -euo pipefail
: "${BENCH_ROOT:?Export the benchmark root before submission}"
: "${REFERENCE:?Export the reference basename before submission}"

python3 "$BENCH_ROOT/plans/audit_reference_taxa.py" \
  "$BENCH_ROOT/references/$REFERENCE.fa.gz" \
  "$BENCH_ROOT/references/audits/$REFERENCE"
