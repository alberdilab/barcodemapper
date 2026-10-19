#!/usr/bin/env bash
#SBATCH --partition=cpuqueue
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=05-taxid-paths-%j.out
#SBATCH --error=05-taxid-paths-%j.err

set -euo pipefail

BENCH_ROOT="${BENCH_ROOT:-/projects/alberdilab/scratch/jpl786/projects/barcodemapper-benchmarks}"
BENCH_ENV="${BENCH_ENV:-$BENCH_ROOT/projects/03-method-comparison/software/env-1.2.3}"
PLAN_DIR="${PLAN_DIR:-$BENCH_ROOT/plans}"
PROJECT_DIR="${PROJECT_DIR:-$BENCH_ROOT/projects/05-default-databases}"
OUT_DIR="${OUT_DIR:-$PROJECT_DIR/taxonomy}"
NCBI_TAXDUMP="${NCBI_TAXDUMP:-/maps/datasets/globe_databases/ncbi/taxonomy/20241001}"
KRAKEN_KINGDOM_TAXDUMP="${KRAKEN_KINGDOM_TAXDUMP:-/maps/datasets/globe_databases/ncbi/taxonomy/20241001}"
KRAKEN_DB="${KRAKEN_DB:-/maps/datasets/globe_databases/kraken2/kraken2_nt/20240530}"
KTAXONOMY="${KTAXONOMY:-$KRAKEN_DB/ktaxonomy.tsv}"
NAME_MAP="${NAME_MAP:-$PLAN_DIR/metadata/05-default-databases/name-map.tsv}"
PYTHON="$BENCH_ENV/bin/python"
BUILDER="$PLAN_DIR/build_ncbi_path_map.py"

for input in "$PYTHON" "$BUILDER" "$NAME_MAP" "$KTAXONOMY" \
             "$NCBI_TAXDUMP/nodes.dmp" "$NCBI_TAXDUMP/names.dmp" \
             "$NCBI_TAXDUMP/merged.dmp" "$NCBI_TAXDUMP/delnodes.dmp" \
             "$KRAKEN_KINGDOM_TAXDUMP/nodes.dmp"; do
    [[ -f "$input" ]] || { echo "Missing input: $input" >&2; exit 2; }
done
[[ ! -e "$OUT_DIR" && ! -L "$OUT_DIR" ]] || { echo "Refusing to overwrite $OUT_DIR" >&2; exit 2; }
mkdir -p "$PROJECT_DIR" "$(dirname "$OUT_DIR")"
stage=$(mktemp -d "$(dirname "$OUT_DIR")/.taxonomy.XXXXXXXX")
trap 'rm -rf "$stage"' EXIT

"$PYTHON" "$BUILDER" \
    --taxdump "$NCBI_TAXDUMP" --name-map "$NAME_MAP" \
    --output "$stage/blast-taxid-paths.tsv" \
    --audit "$stage/blast-taxid-paths.audit.json" \
    --report "$stage/blast-taxid-paths.unresolved.tsv"
"$PYTHON" "$BUILDER" \
    --ktaxonomy "$KTAXONOMY" --kingdom-from-taxdump "$KRAKEN_KINGDOM_TAXDUMP" \
    --name-map "$NAME_MAP" \
    --output "$stage/kraken-taxid-paths.tsv" \
    --audit "$stage/kraken-taxid-paths.audit.json" \
    --report "$stage/kraken-taxid-paths.unresolved.tsv"

[[ ! -e "$OUT_DIR" && ! -L "$OUT_DIR" ]] || { echo "Refusing to overwrite $OUT_DIR" >&2; exit 2; }
mv -T -n "$stage" "$OUT_DIR"
[[ ! -e "$stage" ]] || { echo "Refusing to overwrite $OUT_DIR" >&2; exit 2; }
trap - EXIT
