#!/usr/bin/env bash
# Compare Bowtie2 hit caps (-k N) against an all-hits (-a) baseline on the first
# 1M read pairs of the animal and fungal raw runs, using the patched BarcodeMapper
# overlay (needs -k). One controller job per sample and setting; "a" is the all-hits
# baseline and runs without -k. The input-cutting job runs first.
# Note: the pipeline's mapping resources scale with input size only, so the -a baseline
# at 1M may need Snakemake's automatic retries (memory and time double each attempt).
# Usage: submit_max_hits_comparison.sh
# Set DRY_RUN=1 to echo sbatch commands with fake job IDs and write nothing.
# Set KS="5 50" to override the caps (default: a 1 5 10 50 200).
set -euo pipefail

BENCH_ROOT=/projects/alberdilab/scratch/jpl786/projects/barcodemapper-benchmarks
DRY_RUN=${DRY_RUN:-0}
KS=${KS:-"a 1 5 10 50 200"}
PLAN="$BENCH_ROOT/plans"
REFERENCE_DIR="$BENCH_ROOT/references"
PROJECT01="$BENCH_ROOT/projects/01-mock-communities"
PROJECT04="$BENCH_ROOT/projects/04-scaling"
PROJECT05="$BENCH_ROOT/projects/05-hit-cap"
PAIRS=${PAIRS:-1000000}
OVERLAY="$PROJECT04/overlay-maxhits"
LOG_ROOT="$BENCH_ROOT/tasks/hit-cap"
JOB_MAP="$LOG_ROOT/job-map.tsv"

if [[ ! -f "$OVERLAY/barcodemapper/cli.py" ]] || ! grep -q max_hits "$OVERLAY/barcodemapper/cli.py"; then
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
    echo "Hit-cap comparison was already submitted: $JOB_MAP" >&2
    exit 2
fi
if [[ "$DRY_RUN" != 1 ]]; then mkdir -p "$LOG_ROOT" "$PROJECT05/results/animal" "$PROJECT05/results/fungal"; fi
printf 'sample\tmax_hits\tjob_id\n' > "$JOB_MAP"

INPUTS_JOB=$(sbatch --parsable \
    --export="ALL,BENCH_ROOT=$BENCH_ROOT,PAIRS=$PAIRS" \
    --output="$LOG_ROOT/inputs-%j.out" \
    "$PLAN/make_hit_cap_inputs.sbatch")
printf 'inputs\t%s\t%s\n' "$PAIRS" "$INPUTS_JOB" | tee -a "$JOB_MAP"

submit() {
    local label="$1" run="$2" reference="$3" k="$4"
    local tag="k$k" cap=(-k "$k")
    if [[ "$k" == a ]]; then tag=all; cap=(); fi
    local job
    job=$(sbatch --parsable --dependency="afterok:$INPUTS_JOB" \
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

for k in $KS; do
    submit animal ERR6001161 default "$k"
    submit fungal SRR14351851 all-fungi "$k"
done
