"""Write <run>.manifest.json and <run>.SHA256SUMS for one comparator run."""
import hashlib
import json
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

bench, run, tool = sys.argv[1:4]
bench = Path(bench)
P = bench / "projects/03-method-comparison"
out = P / "results/pilots-1.2.3" / "minimap2-sensitive"
work = P / "work/comparators-1.2.3"
reads = P / f"work/pilots-1.2.3/{run}-pilot-standard/fastp"
ref = {"SRR14351851": "all-fungi"}.get(run, "default")
ref_fa = bench / f"references-1.2.3/{ref}.fa.gz"


def sha(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 24), b""):
            digest.update(block)
    return digest.hexdigest()


ext = "tsv" if tool == "vsearch" else "paf"
envbin = P / "software/env-1.2.3/bin"
version = subprocess.run([str(envbin / tool), "--version"], capture_output=True, text=True)
version = (version.stdout + version.stderr).splitlines()[0]
calls = [line.rstrip("\n").split("\t") for line in open(out / f"{run}.calls.tsv")][1:]
statuses = Counter(row[2] for row in calls)
import gzip
with gzip.open(reads / f"{run}_1.fq.gz", "rt") as handle:
    retained = sum(1 for _ in handle) // 4
counts = {"pairs_in": retained, "converted_rows": len(calls), "status": dict(statuses),
          "pairs_with_hit": retained - statuses.get("no_hit", 0)}
for m in (1, 2):
    counts[f"R{m}_accept_log"] = (out / f"{run}.R{m}.accept.log").read_text().strip()
    counts[f"R{m}_raw_lines"] = sum(1 for _ in open(out / f"{run}.R{m}.raw.{ext}"))
files = [f"{run}.R{m}.{kind}" for m in (1, 2) for kind in (f"raw.{ext}", f"accepted.{ext}", "timing.json", "command.txt", "accept.log")]
files += [f"{run}.accepted.{ext}", f"{run}.hitmap.tsv", f"{run}.calls.tsv"]
manifest = {
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "run": run, "tool": tool, "tool_version": version,
    "reference": {"name": ref, "path": str(ref_fa), "sha256": sha(ref_fa)},
    "config": {"path": str(work / "comparator_config.minimap2-sensitive.yaml"), "sha256": sha(work / "comparator_config.minimap2-sensitive.yaml")},
    "scripts": {name: sha(bench / "tasks/1.2.3-comparators" / name) for name in
                ("accept_hits.py", "make_hit_map.py", "convert_calls.py", "time_classifier.sh",
                 "run_minimap2_sensitive.sh", "finish_minimap2_sensitive.sh", "write_manifest_sensitive.py", "run_comparator.sh", "common.sh")},
    "inputs": {f"fastq_R{m}": {"path": str(reads / f"{run}_{m}.fq.gz"), "sha256": sha(reads / f"{run}_{m}.fq.gz")} for m in (1, 2)},
    "query_fasta": {f"R{m}": {"path": str(work / f"queries/{run}_{m}.fa"), "sha256": sha(work / f"queries/{run}_{m}.fa")} for m in (1, 2)},
    "commands": {f"R{m}": (out / f"{run}.R{m}.command.txt").read_text().strip() for m in (1, 2)},
    "timing": {f"R{m}": json.loads((out / f"{run}.R{m}.timing.json").read_text()) for m in (1, 2)},
    "counts": counts,
    "outputs": {name: sha(out / name) for name in files},
}
(out / f"{run}.manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
with open(out / f"{run}.SHA256SUMS", "x") as handle:
    for name in files + [f"{run}.manifest.json"]:
        handle.write(f"{sha(out / name)}  {name}\n")
print(json.dumps(counts, indent=2))
