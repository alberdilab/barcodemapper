#!/usr/bin/env bash
# BarcodeMapper 1.1.2 resubmission of the reference builds and their dependent jobs.
# Builds have no dependency (source cleaning job 45898609 already completed).
# Usage: submit_clean_source_recovery_1.1.2.sh [ANIMAL_INPUT_JOB_ID]
# Set DRY_RUN=1 to echo sbatch commands with fake job IDs and write nothing.
set -euo pipefail

BENCH_ROOT=/projects/alberdilab/scratch/jpl786/projects/barcodemapper-benchmarks
ANIMAL_INPUT_JOB_ID=${1:-45898272}
DRY_RUN=${DRY_RUN:-0}
PLAN="$BENCH_ROOT/plans"
REFERENCE_DIR="$BENCH_ROOT/references"
PROJECT01="$BENCH_ROOT/projects/01-mock-communities"
PROJECT04="$BENCH_ROOT/projects/04-scaling"
LOG_ROOT="$BENCH_ROOT/tasks/clean-source-recovery"
JOB_MAP="$LOG_ROOT/job-map.tsv"

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
    echo "Recovery jobs were already submitted: $JOB_MAP" >&2
    exit 2
fi
if [[ "$DRY_RUN" != 1 ]]; then mkdir -p "$LOG_ROOT" "$REFERENCE_DIR/audits"; fi
printf 'stage\tjob_id\tdependency\n' > "$JOB_MAP"

record() {
    printf '%s\t%s\t%s\n' "$1" "$2" "$3" | tee -a "$JOB_MAP"
}

submit_build() {
    local name="$1"
    local extra=()
    if [[ "$name" == all-fungi ]]; then
        extra=(-x k__Animalia -y k__Fungi)
    fi
    local job
    job=$(sbatch --parsable \
        --export="ALL,BENCH_ENV=$PROJECT01/software/env" \
        --job-name="bm-build-$name-clean" \
        --output="$LOG_ROOT/build-$name-%j.out" \
        "$PLAN/controller.sbatch" \
        barcodemapper --build --slurm \
        -b "$REFERENCE_DIR/cleaned-sources/bold.fasta" \
        -u "$REFERENCE_DIR/cleaned-sources/unite.fasta" \
        ${extra[@]+"${extra[@]}"} \
        -d "$REFERENCE_DIR/$name.fa.gz" \
        -t "$REFERENCE_DIR/work-$name-clean")
    record "build-$name" "$job" "none" >&2
    printf '%s' "$job"
}

submit_finalize_reference() {
    local name="$1" dependency="$2"
    local job
    job=$(sbatch --parsable --dependency="afterok:$dependency" \
        --export="ALL,BENCH_ROOT=$BENCH_ROOT,REFERENCE=$name" \
        --job-name="bm-freeze-$name-clean" \
        --output="$LOG_ROOT/freeze-$name-%j.out" \
        "$PLAN/finalize_reference.sbatch")
    record "freeze-$name" "$job" "$dependency" >&2
    printf '%s' "$job"
}

submit_audit() {
    local name="$1" dependency="$2"
    local job
    job=$(sbatch --parsable --dependency="afterok:$dependency" \
        --export="ALL,BENCH_ROOT=$BENCH_ROOT,REFERENCE=$name" \
        --job-name="bm-audit-$name-clean" \
        --output="$LOG_ROOT/audit-$name-%j.out" \
        "$PLAN/audit_reference.sbatch")
    record "audit-$name" "$job" "$dependency" >&2
    printf '%s' "$job"
}

submit_pilot() {
    local run="$1" name="$2" dependency="$3" label="$4"
    local task="pilot-$label-clean"
    local job
    job=$(sbatch --parsable --dependency="afterok:$dependency" \
        --export="ALL,BENCH_ENV=$PROJECT01/software/env" \
        --job-name="bm-pilot-$label-clean" \
        --output="$BENCH_ROOT/tasks/$task/controller-%j.out" \
        "$PLAN/controller.sbatch" \
        barcodemapper --slurm \
        -1 "$PROJECT01/inputs/pilot/${run}_1.fastq.gz" \
        -2 "$PROJECT01/inputs/pilot/${run}_2.fastq.gz" \
        -d "$REFERENCE_DIR/$name.fa.gz" \
        -o "$PROJECT01/results/$run-pilot-standard.tsv" \
        -t "$PROJECT01/work/$run-pilot-standard" \
        --standard --adapter-r1 auto --adapter-r2 auto)
    record "pilot-$label" "$job" "$dependency" >&2
    local finalized
    finalized=$(sbatch --parsable --dependency="afterok:$job" \
        --export="ALL,BENCH_ROOT=$BENCH_ROOT,BENCH_ENV=$PROJECT01/software/env,RUN=$run,REFERENCE=$name,CONTROLLER_JOBID=$job,CONTROLLER_TASK=$task" \
        --job-name="bm-verify-$label-clean" \
        --output="$LOG_ROOT/verify-$label-%j.out" \
        "$PLAN/finalize_pilot.sbatch")
    record "verify-$label" "$finalized" "$job" >&2
    printf '%s' "$finalized"
}

if [[ "$DRY_RUN" != 1 ]]; then
    mkdir -p "$BENCH_ROOT/tasks/pilot-fungal-clean" \
        "$BENCH_ROOT/tasks/pilot-animal-clean" \
        "$BENCH_ROOT/tasks/pilot-blank-clean"
fi

DEFAULT_BUILD=$(submit_build default)
FUNGI_BUILD=$(submit_build all-fungi)
DEFAULT_FREEZE=$(submit_finalize_reference default "$DEFAULT_BUILD")
FUNGI_FREEZE=$(submit_finalize_reference all-fungi "$FUNGI_BUILD")
AUDIT_DEFAULT=$(submit_audit default "$DEFAULT_FREEZE")
AUDIT_FUNGI=$(submit_audit all-fungi "$FUNGI_FREEZE")

FUNGI_VERIFY=$(submit_pilot SRR14351851 all-fungi "$FUNGI_FREEZE" fungal)
ANIMAL_VERIFY=$(submit_pilot ERR6001161 default "$DEFAULT_FREEZE:$ANIMAL_INPUT_JOB_ID" animal)
submit_pilot ERR6001173 default "$DEFAULT_FREEZE" blank > /dev/null

EXPECTED=$(sbatch --parsable --dependency="afterok:$AUDIT_DEFAULT:$AUDIT_FUNGI" \
    --export="ALL,BENCH_ROOT=$BENCH_ROOT" \
    --job-name=bm-expected-taxa \
    --output="$LOG_ROOT/expected-taxa-%j.out" \
    "$PLAN/annotate_expected_taxa.sbatch")
record expected-taxa "$EXPECTED" "$AUDIT_DEFAULT:$AUDIT_FUNGI"

submit_besthit() {
    local label="$1" run="$2" dependency="$3"
    local job
    job=$(sbatch --parsable --dependency="afterok:$dependency" \
        --export="ALL,BENCH_ROOT=$BENCH_ROOT,BENCH_ENV=$BENCH_ROOT/projects/03-method-comparison/software/env,RUN=$run" \
        --job-name="bm-besthit-$label" \
        --output="$LOG_ROOT/besthit-$label-%j.out" \
        "$PLAN/best_hit_ablation.sbatch")
    record "besthit-$label" "$job" "$dependency"
}
submit_besthit fungal SRR14351851 "$FUNGI_VERIFY"
submit_besthit animal ERR6001161 "$ANIMAL_VERIFY"

SCALING=$(sbatch --parsable --dependency="afterok:$ANIMAL_VERIFY" \
    --export="ALL,BENCH_ENV=$PROJECT04/software/env" \
    --job-name=bm-scale-100k-clean \
    --output="$LOG_ROOT/scale-100k-%j.out" \
    "$PLAN/controller.sbatch" \
    barcodemapper --slurm \
    -1 "$PROJECT01/inputs/pilot/ERR6001161_1.fastq.gz" \
    -2 "$PROJECT01/inputs/pilot/ERR6001161_2.fastq.gz" \
    -d "$REFERENCE_DIR/default.fa.gz" \
    -o "$PROJECT04/results/background-heavy/100k/rep1.tsv" \
    -t "$PROJECT04/work/background-heavy/100k/rep1" \
    --standard --adapter-r1 auto --adapter-r2 auto)
record scale-100k "$SCALING" "$ANIMAL_VERIFY"

SAMPLER=$(sbatch --parsable --dependency="after:$SCALING" \
    --export="ALL,BENCH_ROOT=$BENCH_ROOT,CONTROLLER_JOBID=$SCALING" \
    --job-name=bm-sample-scratch-clean \
    --output="$LOG_ROOT/sample-scratch-%j.out" \
    "$PLAN/sample_scratch.sbatch")
record sample-scratch "$SAMPLER" "$SCALING"

echo "Submitted recovery chain; job map: $JOB_MAP"
