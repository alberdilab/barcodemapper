#!/usr/bin/env bash
# Extend the 1M-pair hit-cap comparison (see submit_max_hits_comparison.sh):
#  - caps 75, 100 and 150 on the existing animal (ERR6001161) and fungal (SRR14351851) inputs;
#  - four more samples (three other animal mock communities, one more fungal sample) with the
#    all-hits baseline (a) and caps 50 and 200. The input-cutting job runs first for these.
# Usage: submit_hit_cap_extension.sh
# Set DRY_RUN=1 to echo sbatch commands with fake job IDs and write nothing.
set -euo pipefail

BENCH_ROOT=/projects/alberdilab/scratch/jpl786/projects/barcodemapper-benchmarks
DRY_RUN=${DRY_RUN:-0}
PAIRS=${PAIRS:-1000000}
NEW_RUNS="ERR6001164 ERR6001167 ERR6001170 SRR14351849"
PLAN="$BENCH_ROOT/plans"
REFERENCE_DIR="$BENCH_ROOT/references"
PROJECT04="$BENCH_ROOT/projects/04-scaling"
PROJECT05="$BENCH_ROOT/projects/05-hit-cap"
OVERLAY="$PROJECT04/overlay-maxhits"
LOG_ROOT="$BENCH_ROOT/tasks/hit-cap-extension"
JOB_MAP="$LOG_ROOT/job-map.tsv"

if ! grep -q max_hits "$OVERLAY/barcodemapper/cli.py"; then
    echo "Overlay with -k support missing: $OVERLAY" >&2
    exit 2
fi

FAKE_FILE=$(mktemp)
echo 9000000 > "$FAKE_FILE"
sbatch() {
    if [[ "$DRY_RUN" == 1 ]]; then
        local FAKE_ID=$(( $(cat "$FAKE_FILE") + 1 ))
        echo "$FAKE_ID" > "$FAKE_FILE"
        echo "sbatch $*" >&2
        echo "$FAKE_ID"
    else
        command sbatch "$@"
    fi
}

if [[ "$DRY_RUN" == 1 ]]; then
    JOB_MAP=/dev/null
fi
if [[ "$DRY_RUN" != 1 && -e "$JOB_MAP" ]]; then
    echo "Extension was already submitted: $JOB_MAP" >&2
    exit 2
fi
if [[ "$DRY_RUN" != 1 ]]; then
    mkdir -p "$LOG_ROOT"
    for run in $NEW_RUNS; do mkdir -p "$PROJECT05/results/$run"; done
fi
printf 'label\tmax_hits\tjob_id\n' > "$JOB_MAP"

INPUTS_JOB=$(sbatch --parsable \
    --export="ALL,BENCH_ROOT=$BENCH_ROOT,PAIRS=$PAIRS,RUNS=$NEW_RUNS" \
    --job-name=bm-cap-inputs-ext \
    --output="$LOG_ROOT/inputs-%j.out" \
    "$PLAN/make_hit_cap_inputs.sbatch")
printf 'inputs\t%s\t%s\n' "$PAIRS" "$INPUTS_JOB" | tee -a "$JOB_MAP"

# submit LABEL RUN REFERENCE K [DEPENDENCY_JOB]; K "a" runs without -k
submit() {
    local label="$1" run="$2" reference="$3" k="$4" dependency="${5:-}"
    local tag="k$k" cap=(-k "$k") dep=()
    if [[ "$k" == a ]]; then tag=all; cap=(); fi
    if [[ -n "$dependency" ]]; then dep=(--dependency="afterok:$dependency"); fi
    local job
    job=$(sbatch --parsable "${dep[@]}" \
        --export="ALL,BENCH_ENV=$PROJECT04/software/env,PYTHONPATH=$OVERLAY" \
        --job-name="bm-cap-$label-$tag" \
        --output="$LOG_ROOT/controller-$label-$tag-%j.out" \
        "$PLAN/controller.sbatch" \
        barcodemapper --slurm \
        -1 "$PROJECT05/inputs/${run}_1.fastq.gz" \
        -2 "$PROJECT05/inputs/${run}_2.fastq.gz" \
        -d "$REFERENCE_DIR/$reference.fa.gz" \
        -o "$PROJECT05/results/$label/$tag.tsv" \
        -t "$PROJECT05/work/$label/$tag" \
        --standard --adapter-r1 auto --adapter-r2 auto "${cap[@]}")
    printf '%s\t%s\t%s\n' "$label" "$k" "$job" | tee -a "$JOB_MAP"
}

# Step 2: intermediate caps on the existing inputs (already in place).
for k in 75 100 150; do
    submit animal ERR6001161 default "$k"
    submit fungal SRR14351851 all-fungi "$k"
done

# Step 1: more samples, after their inputs are cut.
for run in $NEW_RUNS; do
    reference=default
    if [[ "$run" == SRR* ]]; then reference=all-fungi; fi
    for k in a 50 200; do
        submit "$run" "$run" "$reference" "$k" "$INPUTS_JOB"
    done
done
