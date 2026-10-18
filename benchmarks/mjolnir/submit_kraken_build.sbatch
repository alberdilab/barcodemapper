#!/usr/bin/env bash
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00

set -euo pipefail
: "${BENCH_ROOT:?Export the benchmark root before submission}"
: "${BENCH_ENV:?Export the Project 03 Conda prefix before submission}"
: "${REFERENCE_FASTA:?Export the frozen BarcodeMapper .fa.gz path}"
: "${KRAKEN_DB:?Export the Kraken2 database output directory}"
: "${K:?Export the Kraken2 k-mer length}"
: "${MINIMIZER_LEN:?Export the Kraken2 minimizer length}"
: "${MAX_DB_SIZE:?Export the Kraken2 maximum database size in bytes}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-16}}"

python_cmd=("$BENCH_ENV/bin/python" "$BENCH_ROOT/plans/build_kraken_inputs.py" "$REFERENCE_FASTA" "$KRAKEN_DB")
add_cmd=("$BENCH_ENV/bin/kraken2-build" --add-to-library "$KRAKEN_DB/library.fa" --db "$KRAKEN_DB")
build_cmd=("$BENCH_ENV/bin/kraken2-build" --build --db "$KRAKEN_DB" \
  --threads "$THREADS" --kmer-len "$K" --minimizer-len "$MINIMIZER_LEN" \
  --max-db-size "$MAX_DB_SIZE")

mkdir -p "$KRAKEN_DB"
commands="$KRAKEN_DB/build-commands.sh"
if [[ -e "$commands" ]]; then
  echo "Refusing to overwrite $commands" >&2
  exit 2
fi
{
  printf '#!/usr/bin/env bash\nset -euo pipefail\n'
  printf '%q ' "${python_cmd[@]}"; printf '\n'
  printf '%q ' "${add_cmd[@]}"; printf '\n'
  printf '%q ' "${build_cmd[@]}"; printf '\n'
} > "$commands"

"${python_cmd[@]}"
"${add_cmd[@]}"
"${build_cmd[@]}"
