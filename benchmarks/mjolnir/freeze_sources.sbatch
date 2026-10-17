#!/usr/bin/env bash
#SBATCH --partition=cpuqueue
#SBATCH --qos=normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00

set -euo pipefail
: "${BENCH_ROOT:?Export the benchmark root before submission}"
: "${BENCH_ENV:?Export an absolute Conda prefix before submission}"

export BENCH_ROOT BENCH_ENV
bash -ic '
  set -euo pipefail
  useconda
  conda activate "$BENCH_ENV"
  python "$BENCH_ROOT/plans/freeze_sources.py" "$BENCH_ROOT/references/sources"
  cd "$BENCH_ROOT/references/sources"
  sha256sum -c SHA256SUMS
  chmod a-w bold.fasta unite.fasta source-releases.json SHA256SUMS
'
