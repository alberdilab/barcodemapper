"""Verify a completed 1.2.3 paired-end pilot and write its benchmark run record.

Adapted from plans/finalize_pilot.py (1.1.2): Project 03 result/work paths
(pilots-1.2.3/), references-1.2.3/, 1.2.3 wheel and lock, mode standard|allhits.
"""

import argparse
import hashlib
import json
import re
import shlex
import subprocess
from pathlib import Path

import yaml


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def artifact(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Missing or empty pilot artifact: {path}")
    return {"path": str(path), "bytes": path.stat().st_size,
            "sha256": sha256(path)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bench_root", type=Path)
    parser.add_argument("run")
    parser.add_argument("reference", choices=("default", "all-fungi"))
    parser.add_argument("controller_job_id")
    parser.add_argument("controller_task")
    parser.add_argument("mode", choices=("standard", "allhits"))
    args = parser.parse_args()
    root = args.bench_root.expanduser().resolve()
    p01 = root / "projects" / "01-mock-communities"
    p03 = root / "projects" / "03-method-comparison"
    label = f"{args.run}-pilot-{args.mode}"
    result_base = p03 / "results" / "pilots-1.2.3" / label
    work = p03 / "work" / "pilots-1.2.3" / label
    task_dir = root / "tasks" / "1.2.3-pilots" / args.controller_task
    pilot_manifest = p01 / "inputs" / "pilot" / f"{args.run}.manifest.json"
    pilot = json.loads(pilot_manifest.read_text(encoding="utf-8"))
    if pilot["pairs"] != 100000 or pilot["run_accession"] != args.run:
        raise ValueError("Unexpected pilot input definition")
    for mate in pilot["mates"]:
        if sha256(Path(mate["path"])) != mate["sha256"]:
            raise ValueError("Pilot FASTQ checksum mismatch")
    reference_path = root / "references-1.2.3" / f"{args.reference}.fa.gz"
    snapshot_path = root / "references-1.2.3" / f"{args.reference}.snapshot.json"
    snapshot = json.loads(snapshot_path.read_text(encoding="utf-8"))
    reference_sha = next(item["sha256"] for item in snapshot["artifacts"]
                         if Path(item["path"]).name == reference_path.name)
    if sha256(reference_path) != reference_sha:
        raise ValueError("Reference FASTA checksum mismatch")
    output_paths = [
        result_base.with_suffix(".tsv"),
        result_base.with_suffix(".html"),
        result_base.with_suffix(".manifest.json"),
        work / "taxonomy" / f"{args.run}.tsv",
        work / "extract" / f"{args.run}.tsv",
        work / "fastp" / f"{args.run}_1.fq.gz",
        work / "fastp" / f"{args.run}_2.fq.gz",
        work / "fastp" / f"{args.run}.json",
        task_dir / f"controller-{args.controller_job_id}.out",
    ]
    output_records = [artifact(path) for path in output_paths]
    run_manifest = json.loads(result_base.with_suffix(".manifest.json").read_text(encoding="utf-8"))
    if run_manifest["database"]["sha256"] != reference_sha:
        raise ValueError("BarcodeMapper run used a different reference")
    if run_manifest["barcodemapper_version"] != "1.2.3":
        raise ValueError("Unexpected BarcodeMapper version")
    fastp = json.loads((work / "fastp" / f"{args.run}.json").read_text(encoding="utf-8"))
    retained_pairs = fastp["summary"]["after_filtering"]["total_reads"] // 2
    records_path = result_base.with_suffix(".results.json")
    checksums_path = result_base.with_suffix(".SHA256SUMS")
    record_path = result_base.with_suffix(".run.yaml")
    if any(path.exists() for path in (records_path, checksums_path, record_path)):
        parser.error("Pilot finalization artifacts already exist")
    records_path.write_text(json.dumps({"run": args.run, "mode": args.mode,
                                        "fastp_retained_pairs": retained_pairs,
                                        "artifacts": output_records},
                                       indent=2) + "\n", encoding="utf-8")
    checksums_path.write_text("".join(f"{item['sha256']}  {item['path']}\n"
                                     for item in output_records), encoding="utf-8")
    controller_log = output_paths[-1].read_text(encoding="utf-8", errors="replace")
    worker_ids = sorted(set(re.findall(
        r"(?:external jobid ['\"]?|SLURM jobid )(\d+)", controller_log, re.I
    )))
    sacct = subprocess.run([str(root / "software" / "sacct-wrapper.sh"),
                            "-j", args.controller_job_id,
                            "--format=NodeList", "-n", "-P"],
                           capture_output=True, text=True, check=False)
    node = next((line.strip() for line in sacct.stdout.splitlines()
                 if line.strip() and line.strip() != "None"), None)
    command = [
        "barcodemapper", "--slurm", "-1", str(pilot["mates"][0]["path"]),
        "-2", str(pilot["mates"][1]["path"]), "-d", str(reference_path),
        "-o", str(result_base.with_suffix(".tsv")), "-t", str(work),
        "--standard", "--adapter-r1", "auto", "--adapter-r2", "auto",
    ]
    if args.mode == "allhits":
        command.append("--all-hits")
    record = yaml.safe_load((root / "plans" / "run-record.template.yaml").read_text(encoding="utf-8"))
    wheel = p03 / "software" / "barcodemapper-1.2.3-py3-none-any.whl"
    wheel_sha = (p03 / "software" / "wheel.1.2.3.sha256").read_text().split()[0]
    if sha256(wheel) != wheel_sha:
        raise ValueError("Release wheel checksum mismatch")
    record["software"]["barcodemapper_version"] = "1.2.3"
    record["software"]["barcodemapper_wheel_sha256"] = wheel_sha
    record["software"]["release_commit"] = None
    record["project"] = "03-method-comparison"
    record["run_id"] = label
    record["run_utc"] = run_manifest["created_utc"]
    record["software"]["conda_explicit_lock_sha256"] = sha256(p03 / "software" / "conda-linux-64.1.2.3.explicit.txt")
    record["reference"] = {
        "name": args.reference,
        "fasta_sha256": reference_sha,
        "source_release_record": str(root / "references" / "sources" / "source-releases.json"),
        "source_release_record_sha256": sha256(root / "references" / "sources" / "source-releases.json"),
        "snapshot_manifest": str(snapshot_path),
        "snapshot_manifest_sha256": sha256(snapshot_path),
        "source_filter": snapshot.get("source_filter"),
    }
    record["inputs_manifest"] = str(pilot_manifest)
    record["input_manifest_sha256"] = sha256(pilot_manifest)
    record["random_seed"] = None
    record["command"] = shlex.join(command)
    record["slurm_job_ids"] = [args.controller_job_id] + worker_ids
    record["slurm_partition"] = "cpuqueue"
    record["compute_node"] = node
    record["threads_or_cpus"] = "controller 1; worker rules capped at 8"
    record["memory_limit_gb"] = "controller 4; worker rule-specific"
    record["wall_time_limit"] = "2-00:00:00 controller"
    record["results_manifest"] = str(records_path)
    record["results_sha256_manifest"] = str(checksums_path)
    record["notes"] = (f"1.2.3 pilot, mode={args.mode}; first 100,000 synchronized ENA read pairs; "
                       f"fastp retained pairs {retained_pairs}; preliminary pilot")
    record_path.write_text(yaml.safe_dump(record, sort_keys=False), encoding="utf-8")
    print(f"Verified pilot and wrote {record_path}")


if __name__ == "__main__":
    main()
