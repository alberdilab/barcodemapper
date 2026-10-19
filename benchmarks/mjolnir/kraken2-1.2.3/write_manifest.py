"""Write <run>.manifest.json and <run>.SHA256SUMS for one Kraken2 run (new files only)."""
import gzip, hashlib, json, subprocess, sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

bench, run, setting = sys.argv[1:4]
bench = Path(bench)
P = bench / "projects/03-method-comparison"
out = P / "results/pilots-1.2.3" / setting
work = P / "work/comparators-1.2.3"
reads = P / f"work/pilots-1.2.3/{run}-pilot-standard/fastp"
dbname = {"SRR14351851": "all-fungi"}.get(run, "default")
db = P / "work/kraken2" / dbname
tasks = bench / "tasks/1.2.3-kraken2"
conf = {"kraken2": "0.0", "kraken2-conf0.1": "0.1"}[setting]


def sha(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 24), b""):
            digest.update(block)
    return digest.hexdigest()


def fq_ids(path):
    with gzip.open(path, "rt") as handle:
        for i, line in enumerate(handle):
            if i % 4 == 0:
                tok = line[1:].split()[0]
                yield tok[:-2] if tok.endswith(("/1", "/2")) else tok


# input hashes must equal the pilot manifest SHA256SUMS
expected = {}
for line in (P / f"results/pilots-1.2.3/{run}-pilot-standard.SHA256SUMS").read_text().splitlines():
    digest, _, name = line.partition("  ")
    expected[Path(name).name] = digest
inputs = {}
for m in (1, 2):
    p = reads / f"{run}_{m}.fq.gz"
    h = sha(p)
    assert h == expected[p.name], f"{p} sha256 differs from pilot manifest"
    inputs[f"fastq_R{m}"] = {"path": str(p), "sha256": h, "matches_pilot_SHA256SUMS": True}
ids1, ids2 = list(fq_ids(reads / f"{run}_1.fq.gz")), list(fq_ids(reads / f"{run}_2.fq.gz"))
assert ids1 == ids2, "R1/R2 IDs differ or are out of order"
assert len(ids1) == len(set(ids1)), "duplicate IDs"

version = subprocess.run([str(P / "software/env-1.2.3/bin/kraken2"), "--version"], capture_output=True, text=True)
version = version.stdout.splitlines()[0]
rows = [l.rstrip("\n").split("\t") for l in open(out / f"{run}.calls.tsv")][1:]
statuses = Counter(r[2] for r in rows)
assert len(rows) == len(ids1), "converted rows != retained pairs"
assert [r[0] for r in rows] == ids1
with open(out / f"{run}.output.txt") as handle:
    codes = Counter(l[0] for l in handle)
counts = {"retained_pairs": len(ids1), "converted_rows": len(rows), "calls_status": dict(statuses),
          "kraken2_output_rows_by_code": dict(codes), "R1_R2_ids_identical_in_order": True}
files = [f"{run}.{k}" for k in ("output.txt", "report.txt", "command.txt", "timing.json", "kraken2.log", "calls.tsv")]
dbfiles = ("hash.k2d", "opts.k2d", "taxo.k2d", "taxid-to-path.tsv")
sc = ("common.sh", "run_kraken2.sbatch", "run_kraken2.sh", "finish_kraken2.sbatch", "finish_kraken2.sh",
      "write_manifest.py", "convert_calls.py", "time_classifier.sh", "submit_kraken2.sh")
manifest = {
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "run": run, "setting": setting, "confidence": conf, "tool_version": version,
    "database": {"name": dbname, "path": str(db), "files": {f: sha(db / f) for f in dbfiles}},
    "config": {"path": str(work / "comparator_config.kraken2.yaml"), "sha256": sha(work / "comparator_config.kraken2.yaml")},
    "scripts": {n: sha(tasks / n) for n in sc},
    "inputs": inputs,
    "command": (out / f"{run}.command.txt").read_text().strip(),
    "timing": json.loads((out / f"{run}.timing.json").read_text()),
    "counts": counts,
    "outputs": {n: sha(out / n) for n in files},
}
(out / f"{run}.manifest.json").open("x").write(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
with open(out / f"{run}.SHA256SUMS", "x") as handle:
    for name in files + [f"{run}.manifest.json"]:
        handle.write(f"{sha(out / name)}  {name}\n")
print(setting, run, json.dumps(counts))
