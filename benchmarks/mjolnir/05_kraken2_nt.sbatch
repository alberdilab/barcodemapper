#!/usr/bin/env bash
#SBATCH --partition=cpuqueue
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=1150G
#SBATCH --time=7-00:00:00
#SBATCH --output=05-kraken2-nt-%j.out
#SBATCH --error=05-kraken2-nt-%j.err
# A 1.1 TB request will queue until a 2 TB node has enough free memory.
# To prefer a particular eligible node, pass --nodelist to sbatch at submission.

set -euo pipefail
: "${RUN:?Export the pilot run accession before submission}"

BENCH_ROOT="${BENCH_ROOT:-/projects/alberdilab/scratch/jpl786/projects/barcodemapper-benchmarks}"
BENCH_ENV="${BENCH_ENV:-$BENCH_ROOT/projects/03-method-comparison/software/env-1.2.3}"
PLAN_DIR="${PLAN_DIR:-$BENCH_ROOT/plans}"
PROJECT_DIR="${PROJECT_DIR:-$BENCH_ROOT/projects/05-default-databases}"
OUT_DIR="${OUT_DIR:-$PROJECT_DIR/kraken2/$RUN}"
FASTP_ROOT="${FASTP_ROOT:-$BENCH_ROOT/projects/01-mock-communities/work}"
R1="${R1:-$FASTP_ROOT/$RUN-pilot-standard/fastp/${RUN}_1.fq.gz}"
R2="${R2:-$FASTP_ROOT/$RUN-pilot-standard/fastp/${RUN}_2.fq.gz}"
KRAKEN_DB="${KRAKEN_DB:-/maps/datasets/globe_databases/kraken2/kraken2_nt/20240530}"
NT_MD5="${NT_MD5:-$KRAKEN_DB/nt.md5}"
TAXID_PATHS="${TAXID_PATHS:-$PROJECT_DIR/taxonomy/kraken-taxid-paths.tsv}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-32}}"
PRELOAD="${PRELOAD:-1}"
MEM="${MEM:-1150G}" # Submission override: sbatch --mem="$MEM"; recorded below.
PYTHON="$BENCH_ENV/bin/python"
KRAKEN2="$BENCH_ENV/bin/kraken2"
export PATH="$BENCH_ENV/bin:$PATH"

[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "THREADS must be a positive integer" >&2; exit 2; }
if [[ -n "${SLURM_CPUS_PER_TASK:-}" && "$THREADS" != "$SLURM_CPUS_PER_TASK" ]]; then
    echo "THREADS must equal SLURM_CPUS_PER_TASK for time_classifier.sh" >&2
    exit 2
fi
[[ "$PRELOAD" == 0 || "$PRELOAD" == 1 ]] || { echo "PRELOAD must be 0 or 1" >&2; exit 2; }
for input in "$PYTHON" "$KRAKEN2" "$PLAN_DIR/time_classifier.sh" \
             "$PLAN_DIR/convert_calls.py" "$TAXID_PATHS" "$R1" "$R2" \
             "$KRAKEN_DB/hash.k2d" "$KRAKEN_DB/opts.k2d" \
             "$KRAKEN_DB/taxo.k2d" "$KRAKEN_DB/ktaxonomy.tsv" "$NT_MD5"; do
    [[ -f "$input" ]] || { echo "Missing input: $input" >&2; exit 2; }
done
[[ ! -e "$OUT_DIR" && ! -L "$OUT_DIR" ]] || { echo "Refusing to overwrite $OUT_DIR" >&2; exit 2; }
mkdir -p "$PROJECT_DIR" "$(dirname "$OUT_DIR")"
stage=$(mktemp -d "$(dirname "$OUT_DIR")/.kraken2-${RUN}.XXXXXXXX")
preserve_failed_stage() {
    local status=$? failed
    trap - EXIT
    set +e
    if [[ -d "$stage" ]]; then
        failed="${OUT_DIR}.failed.${SLURM_JOB_ID:-manual-$$}"
        if [[ ! -e "$failed" && ! -L "$failed" ]] && \
                mv -T -n "$stage" "$failed" && [[ ! -e "$stage" ]]; then
            echo "Job failed (exit $status); preserved outputs at $failed" >&2
        else
            echo "Job failed (exit $status); could not move stage, outputs remain at $stage" >&2
        fi
    fi
    exit "$status"
}
trap preserve_failed_stage EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

{
    printf 'database=%s\nnt_md5=%s\n' "$KRAKEN_DB" "$NT_MD5"
    printf 'sha256 (hash.k2d deliberately excluded):\n'
    sha256sum "$KRAKEN_DB/opts.k2d" "$KRAKEN_DB/taxo.k2d" \
        "$KRAKEN_DB/ktaxonomy.tsv" "$NT_MD5"
    printf 'directory listing:\n'
    LC_ALL=C ls -ld "$KRAKEN_DB"
    LC_ALL=C ls -l "$KRAKEN_DB"
} > "$stage/input-manifest.txt"
sha256sum "$R1" "$R2" "$TAXID_PATHS" > "$stage/input-sha256.txt"

# This NFS pre-read is best effort: loading the hash into Kraken2 may evict
# cached pages. Classifier-only time excludes the pre-read itself but can still
# include a network read of hash.k2d.
preload_seconds=null
if [[ "$PRELOAD" == 1 ]]; then
    [[ -x /usr/bin/time ]] || { echo "GNU /usr/bin/time is required for PRELOAD=1" >&2; exit 2; }
    /usr/bin/time -f '%e' -o "$stage/preload-seconds.txt" \
        cat "$KRAKEN_DB/hash.k2d" > /dev/null
    preload_seconds=$(cat "$stage/preload-seconds.txt")
    rm "$stage/preload-seconds.txt"
fi
"$PYTHON" - "$stage/preload.json" "$PRELOAD" "$preload_seconds" \
    "$KRAKEN_DB/hash.k2d" "$MEM" <<'PY'
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

output, preload, seconds, hashed_file, requested_mem = sys.argv[1:]
Path(output).write_text(json.dumps({
    "recorded_utc": datetime.now(timezone.utc).isoformat(),
    "preload": "best_effort" if preload == "1" else "disabled",
    "preload_enabled": preload == "1",
    "preload_seconds": float(seconds) if preload == "1" else None,
    "preload_file": hashed_file if preload == "1" else None,
    "classifier_only_excludes_preload": True,
    "requested_mem": requested_mem,
    "slurm_mem_per_node_mb": os.environ.get("SLURM_MEM_PER_NODE"),
}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
PY

# Kraken2 loads the database into memory. Do not add --memory-mapping to this
# timing run: an uncached NFS first read may dominate the measurement.
"$PLAN_DIR/time_classifier.sh" --label classifier-only --cpus "$THREADS" \
    --input "$R1" --input "$R2" --index "$KRAKEN_DB" \
    --output "$stage/classifier-timing.json" -- \
    "$KRAKEN2" --db "$KRAKEN_DB" --paired --threads "$THREADS" \
    --report "$stage/kraken-report.tsv" --output "$stage/kraken-output.tsv" \
    "$R1" "$R2"
"$PYTHON" "$PLAN_DIR/convert_calls.py" kraken2 \
    --read-ids "$R1" --input "$stage/kraken-output.tsv" \
    --taxid-paths "$TAXID_PATHS" --output "$stage/calls.tsv"

[[ ! -e "$OUT_DIR" && ! -L "$OUT_DIR" ]] || { echo "Refusing to overwrite $OUT_DIR" >&2; exit 2; }
mv -T -n "$stage" "$OUT_DIR"
[[ ! -e "$stage" ]] || { echo "Refusing to overwrite $OUT_DIR" >&2; exit 2; }
trap - EXIT
