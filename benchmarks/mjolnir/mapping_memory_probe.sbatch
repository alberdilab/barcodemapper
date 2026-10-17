#!/usr/bin/env bash
# Probe Bowtie2 mapping memory, time and disk for ONE case: the first PAIRS pairs of
# the deep ERR6001161 raw run through the patched BarcodeMapper mapping pipeline
# (bowtie2 | samtools view -b, no samtools collate), optionally with -a, recording
# peak RSS per process. Submit one job per case (see PAIRS, LABEL, FLAGS) so cases run
# in parallel and one failure keeps the others. Submit by hand after review.
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=1-12:00:00
#SBATCH --job-name=bm-map-probe
#SBATCH --output=map-probe-%x-%j.out

set -euo pipefail
: "${BENCH_ROOT:?Export the benchmark root before submission}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
: "${PAIRS:?Export PAIRS (read pairs) before submission}"
: "${LABEL:?Export LABEL (with-a or without-a) before submission}"
FLAGS="${FLAGS:-}"  # extra bowtie2 flags, e.g. -a

project="$BENCH_ROOT/projects/04-scaling"
reads="$BENCH_ROOT/projects/01-mock-communities/inputs/raw"
index="$BENCH_ROOT/references/default"
env_bin="$project/software/env/bin"
out="$project/results/mapping-memory-probe/$SLURM_JOB_ID-p$PAIRS-$LABEL"
tmp="$project/work/mapping-memory-probe/$SLURM_JOB_ID-p$PAIRS-$LABEL"
mkdir -p "$out" "$tmp"
export PATH="$env_bin:$PATH"

make_subset() {  # $1 = pairs; head closes the pipe early, so ignore SIGPIPE statuses
    local n="$1" mate
    for mate in 1 2; do
        { zcat "$reads/ERR6001161_$mate.fastq.gz" | head -n $((n * 4)) | gzip -1 \
            > "$tmp/p${n}_$mate.fq.gz"; } || true
    done
    sha256sum "$tmp"/p${n}_*.fq.gz > "$out/p${n}.input.SHA256SUMS"
}

sample_rss() {  # $1 = label; runs until its sentinel file disappears
    local label="$1" sentinel="$2"
    echo "utc,pid,comm,rss_kb" > "$out/$label.rss.csv"
    while [[ -e "$sentinel" ]]; do
        ps -u "$USER" -o pid=,comm=,rss= | awk -v t="$(date -u +%FT%TZ)" \
            '$2 ~ /bowtie2|samtools/ {print t","$1","$2","$3}' >> "$out/$label.rss.csv"
        echo "$(date -u +%FT%TZ),$(du -sb "$tmp" | cut -f1)" >> "$out/$label.tmp-bytes.csv"
        sleep 5
    done
}

run_case() {  # $1 = pairs, $2 = label, remaining = extra bowtie2 flags
    local n="$1" label="p$1-$2"; shift 2
    local sentinel="$tmp/$label.running"
    touch "$sentinel"
    sample_rss "$label" "$sentinel" &
    local sampler=$!
    local start=$SECONDS
    bash -c 'set -o pipefail
                 bowtie2 --time --threads "$1" -x "$2" -1 "$3" -2 "$4" "${@:6}" --seed 1337 \
                   | samtools view -b -@ "$1" -o "$5" -' \
        _ "$THREADS" "$index" "$tmp/p${n}_1.fq.gz" "$tmp/p${n}_2.fq.gz" "$tmp/$label.bam" "$@" \
        2> "$out/$label.bowtie2.log" || echo "exit $?" > "$out/$label.failed"
    echo "elapsed_seconds $((SECONDS - start))" > "$out/$label.time.txt"
    rm -f "$sentinel"; wait "$sampler" || true
    [[ -s "$tmp/$label.bam" ]] || return 0  # failed case: keep its time and RSS records
    stat -c '%s' "$tmp/$label.bam" > "$out/$label.bam-bytes.txt"
    samtools view -c "$tmp/$label.bam" > "$out/$label.alignment-rows.txt"
}

make_subset "$PAIRS"
# shellcheck disable=SC2086  # FLAGS is a deliberate word-split flag list
run_case "$PAIRS" "$LABEL" $FLAGS

echo "Wrote $out"
if compgen -G "$out/*.failed" > /dev/null; then
    echo "Case failed: $(cat "$out"/*.failed)" >&2
    exit 1
fi
