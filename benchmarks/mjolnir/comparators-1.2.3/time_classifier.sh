#!/usr/bin/env bash
# Time one classifier command on a fixed set of allowed CPUs. Pass the same
# thread count to the classifier itself, and list the FASTQs and index files
# (or directories) whose logical sizes should appear in the JSON record.
#
# Example:
#   time_classifier.sh --label classifier-only --cpus 8 \
#     --input reads_1.fq.gz --input reads_2.fq.gz --index blast-db/ \
#     --output blast-timing.json -- blastn -num_threads 8 ...
set -euo pipefail

usage() {
    echo "Usage: $0 --label classifier-only|end-to-end --cpus N --input PATH [--input PATH ...] --index PATH [--index PATH ...] --output JSON -- COMMAND [ARG ...]" >&2
    exit 2
}

label=""
cpus=""
output=""
inputs=()
indexes=()
while (($#)); do
    case "$1" in
        --label|--cpus|--input|--index|--output)
            (($# >= 2)) || usage
            case "$1" in
                --label) label="$2" ;;
                --cpus) cpus="$2" ;;
                --input) inputs+=("$2") ;;
                --index) indexes+=("$2") ;;
                --output) output="$2" ;;
            esac
            shift 2 ;;
        --) shift; break ;;
        *) usage ;;
    esac
done
[[ "$label" == "classifier-only" || "$label" == "end-to-end" ]] || usage
[[ "$cpus" =~ ^[1-9][0-9]*$ ]] || usage
[[ -n "$output" && ${#inputs[@]} -gt 0 && ${#indexes[@]} -gt 0 && $# -gt 0 ]] || usage
if [[ -n "${SLURM_CPUS_PER_TASK:-}" && "$SLURM_CPUS_PER_TASK" != "$cpus" ]]; then
    echo "--cpus must match SLURM_CPUS_PER_TASK ($SLURM_CPUS_PER_TASK)" >&2
    exit 2
fi
[[ ! -e "$output" ]] || { echo "Timing record already exists: $output" >&2; exit 2; }
[[ -x /usr/bin/time ]] || { echo "GNU /usr/bin/time is required" >&2; exit 2; }
command -v taskset >/dev/null || { echo "taskset is required" >&2; exit 2; }

scratch=$(mktemp -d)
trap 'rm -rf "$scratch"' EXIT
metadata="$scratch/metadata.json"
timing="$scratch/time.txt"
cpu_list=$(python3 - "$cpus" "$metadata" "${#inputs[@]}" "${inputs[@]}" "${indexes[@]}" <<'PY'
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

count = int(sys.argv[1])
metadata = Path(sys.argv[2])
n_inputs = int(sys.argv[3])
inputs = [Path(value).expanduser().resolve(strict=True) for value in sys.argv[4:4 + n_inputs]]
indexes = [Path(value).expanduser().resolve(strict=True) for value in sys.argv[4 + n_inputs:]]
allowed = sorted(os.sched_getaffinity(0))
if len(allowed) < count:
    raise SystemExit(f"Requested {count} CPUs, but only {len(allowed)} are available")

def size(path):
    if path.is_file():
        return path.stat().st_size
    if not path.is_dir():
        raise ValueError(f"Not a regular file or directory: {path}")
    return sum((Path(root) / name).stat().st_size
               for root, _, names in os.walk(path) for name in names)

input_sizes = [{"path": str(path), "bytes": size(path)} for path in inputs]
index_sizes = [{"path": str(path), "bytes": size(path)} for path in indexes]
metadata.write_text(json.dumps({
    "started_utc": datetime.now(timezone.utc).isoformat(),
    "input_sizes": input_sizes,
    "index_sizes": index_sizes,
    "input_bytes": sum(item["bytes"] for item in input_sizes),
    "index_bytes": sum(item["bytes"] for item in index_sizes),
    "cpu_affinity": allowed[:count],
}), encoding="utf-8")
print(",".join(map(str, allowed[:count])))
PY
)

set +e
taskset -c "$cpu_list" /usr/bin/time -v -o "$timing" -- "$@"
status=$?
set -e

mkdir -p "$(dirname "$output")"
python3 - "$metadata" "$timing" "$label" "$cpus" "$status" "$output" <<'PY'
import json
import os
import platform
import socket
import sys
from datetime import datetime, timezone
from pathlib import Path

metadata, timing, label, cpus, status, output = sys.argv[1:]
record = json.loads(Path(metadata).read_text(encoding="utf-8"))
values = {}
for line in Path(timing).read_text(encoding="utf-8").splitlines():
    key, separator, value = line.strip().partition(": ")
    if separator:
        values[key] = value

def required(key):
    if key not in values:
        raise ValueError(f"Missing GNU time field: {key}")
    return values[key]

wall = required("Elapsed (wall clock) time (h:mm:ss or m:ss)")
parts = [float(part) for part in wall.split(":")]
if len(parts) == 2:
    wall_seconds = parts[0] * 60 + parts[1]
elif len(parts) == 3:
    wall_seconds = parts[0] * 3600 + parts[1] * 60 + parts[2]
else:
    raise ValueError(f"Unexpected GNU time elapsed value: {wall}")

cpu_model = platform.processor()
cpuinfo = Path("/proc/cpuinfo")
if cpuinfo.exists():
    for line in cpuinfo.read_text(encoding="utf-8").splitlines():
        if line.startswith("model name"):
            cpu_model = line.partition(":")[2].strip()
            break
user_seconds = float(required("User time (seconds)"))
system_seconds = float(required("System time (seconds)"))
record.update({
    "finished_utc": datetime.now(timezone.utc).isoformat(),
    "label": label,
    "cpus": int(cpus),
    "node": socket.gethostname(),
    "cpu_model": cpu_model,
    "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
    "slurm_array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID"),
    "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
    "slurm_step_id": os.environ.get("SLURM_STEP_ID"),
    "wall_seconds": wall_seconds,
    "user_cpu_seconds": user_seconds,
    "system_cpu_seconds": system_seconds,
    "cpu_seconds": user_seconds + system_seconds,
    "cpu_hours": (user_seconds + system_seconds) / 3600,
    "peak_rss_kib": int(required("Maximum resident set size (kbytes)")),
    "exit_status": int(status),
})
Path(output).write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
PY
exit "$status"
