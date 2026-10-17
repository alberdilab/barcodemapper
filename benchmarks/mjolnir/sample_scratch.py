"""Sample a scaling run's work-directory occupancy until its controller exits."""

import argparse
import csv
import json
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def size_bytes(path):
    if not path.exists():
        return 0
    result = subprocess.run(["du", "-sb", str(path)], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            universal_newlines=True, check=True)
    return int(result.stdout.split(None, 1)[0])


def running(job_id):
    result = subprocess.run(["squeue", "-h", "-j", job_id, "-o", "%T"],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            universal_newlines=True, check=True)
    return bool(result.stdout.strip())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work_dir", type=Path)
    parser.add_argument("output_prefix", type=Path)
    parser.add_argument("controller_job_id")
    parser.add_argument("--interval", type=int, default=30)
    args = parser.parse_args()
    if args.interval < 5:
        parser.error("Sampling interval must be at least five seconds")
    work = args.work_dir.expanduser().resolve()
    prefix = args.output_prefix.expanduser().resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    csv_path = Path(str(prefix) + ".scratch.csv")
    summary_path = Path(str(prefix) + ".scratch.json")
    if csv_path.exists() or summary_path.exists():
        parser.error("Scratch sample output already exists")
    observed_running = False
    peak = 0
    samples = 0
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("utc", "bytes"))
        while True:
            active = running(args.controller_job_id)
            observed_running = observed_running or active
            occupied = size_bytes(work)
            peak = max(peak, occupied)
            writer.writerow((utc_now(), occupied))
            handle.flush()
            samples += 1
            if not active:
                break
            time.sleep(args.interval)
    summary_path.write_text(json.dumps({
        "controller_job_id": args.controller_job_id,
        "work_dir": str(work),
        "sampling_interval_seconds": args.interval,
        "samples": samples,
        "observed_controller_in_queue": observed_running,
        "observed_peak_bytes": peak,
        "note": "Sampled work-directory occupancy; transient peaks between samples may be higher",
    }, indent=2) + "\n", encoding="utf-8")
    print(f"Observed peak scratch occupancy: {peak} bytes over {samples} samples")


if __name__ == "__main__":
    main()
