#!/usr/bin/env bash
# End-to-end test of streaming extraction (no BAM): run the real pipeline with the
# overlay-stream package and compare its final tables with the existing BAM-based runs.
#  - animal and fungal 1M pairs at -k 200 (results/stream/<sample>/k200.tsv), compare with
#    results/<sample>/k200.tsv from the hit-cap benchmark;
#  - animal 100k pairs with all hits, once with overlay-maxhits (BAM) and once with
#    overlay-stream, to cover the -a path in minutes.
# Usage: submit_stream_test.sh
set -euo pipefail
BENCH_ROOT=/projects/alberdilab/scratch/jpl786/projects/barcodemapper-benchmarks
PLAN="$BENCH_ROOT/plans"
REF="$BENCH_ROOT/references"
P04="$BENCH_ROOT/projects/04-scaling"
P05="$BENCH_ROOT/projects/05-hit-cap"
LOG="$BENCH_ROOT/tasks/stream-test"
OUT="$P05/results/stream"
mkdir -p "$LOG" "$OUT"/{animal,fungal,animal-100k-bam,animal-100k-stream} "$P05/inputs-100k"

CUT=$(sbatch --parsable --job-name=bm-stream-cut --time=00:20:00 --mem=2G --output="$LOG/cut-%j.out" --wrap \
  "set -euo pipefail; for r in 1 2; do zcat $P05/inputs/ERR6001161_\$r.fastq.gz | head -n 400000 | gzip > $P05/inputs-100k/ERR6001161_\$r.fastq.gz; done")

submit() {  # label overlay run reference outdir tag workdir dep k
  local label=$1 overlay=$2 run=$3 reference=$4 inputs=$5 outdir=$6 tag=$7 dep=$8 cap=${9:-}
  local capargs=(); [[ -n "$cap" ]] && capargs=(-k "$cap")
  local depargs=(); [[ -n "$dep" ]] && depargs=(--dependency="afterok:$dep")
  sbatch --parsable "${depargs[@]}" \
    --export="ALL,BENCH_ENV=$P04/software/env,PYTHONPATH=$P04/$overlay" \
    --job-name="bm-stream-$label" --output="$LOG/controller-$label-%j.out" \
    "$PLAN/controller.sbatch" barcodemapper --slurm \
    -1 "$inputs/${run}_1.fastq.gz" -2 "$inputs/${run}_2.fastq.gz" \
    -d "$REF/$reference.fa.gz" -o "$outdir/$tag.tsv" -t "$P05/work/stream-$label" \
    --standard --adapter-r1 auto --adapter-r2 auto "${capargs[@]}"
}
printf 'label\tjob\n' > "$LOG/job-map.tsv"
for spec in "animal-k200 overlay-stream ERR6001161 default $P05/inputs $OUT/animal k200  '' 200" \
            "fungal-k200 overlay-stream SRR14351851 all-fungi $P05/inputs $OUT/fungal k200 '' 200" \
            "animal100k-bam overlay-maxhits ERR6001161 default $P05/inputs-100k $OUT/animal-100k-bam all $CUT ''" \
            "animal100k-stream overlay-stream ERR6001161 default $P05/inputs-100k $OUT/animal-100k-stream all $CUT ''"; do
  eval "set -- $spec"
  job=$(submit "$@")
  printf '%s\t%s\n' "$1" "$job" | tee -a "$LOG/job-map.tsv"
done
