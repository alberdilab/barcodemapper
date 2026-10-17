#!/usr/bin/env bash
# Cut the first PAIRS read pairs of the deep raw runs listed in RUNS for the hit-cap
# comparison (default: the animal ERR6001161 and fungal SRR14351851 runs). Deterministic prefixes, so the 100k pilot inputs are nested inside them.
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=02:00:00
#SBATCH --job-name=bm-cap-inputs

set -euo pipefail
: "${BENCH_ROOT:?Export the benchmark root before submission}"
PAIRS="${PAIRS:-1000000}"
RUNS="${RUNS:-ERR6001161 SRR14351851}"
raw="$BENCH_ROOT/projects/01-mock-communities/inputs/raw"
out="$BENCH_ROOT/projects/05-hit-cap/inputs"
mkdir -p "$out"

for run in $RUNS; do
    for mate in 1 2; do
        target="$out/${run}_${mate}.fastq.gz"
        if [[ -e "$target" ]]; then echo "Refusing to overwrite $target" >&2; exit 2; fi
        # head closes the pipe early, so ignore SIGPIPE statuses from zcat
        { zcat "$raw/${run}_${mate}.fastq.gz" | head -n $((PAIRS * 4)) | gzip -6 > "$target"; } || true
        lines=$(zcat "$target" | wc -l)
        if (( lines != PAIRS * 4 )); then echo "$target has $lines lines, expected $((PAIRS * 4))" >&2; exit 3; fi
    done
    # The two mates must carry matching read IDs in order.
    diff <(zcat "$out/${run}_1.fastq.gz" | awk 'NR%4==1{print $1}') \
         <(zcat "$out/${run}_2.fastq.gz" | awk 'NR%4==1{print $1}') > /dev/null \
        || { echo "Mate read IDs differ for $run" >&2; exit 4; }
done
(cd "$out" && sha256sum ./*.fastq.gz > inputs.SHA256SUMS)
printf '{"pairs": %d, "source": "first pairs of the raw runs", "created_utc": "%s"}\n' \
    "$PAIRS" "$(date -u +%FT%TZ)" > "$out/inputs.json"
echo "Wrote $out"
