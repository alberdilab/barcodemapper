"""Verify one Kraken2 DB and write db-manifest.json + db-manifest.SHA256SUMS (new files only)."""
import hashlib, json, os, subprocess, sys
from datetime import datetime, timezone
from pathlib import Path
bench, name = Path(sys.argv[1]), sys.argv[2]
db = bench / "projects/03-method-comparison/work/kraken2" / name
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()
a = json.loads((db / "audit.json").read_text())
ok = (a["source_records"] == a["library_records"] == a["map_records"]
      and a["source_records_by_marker"] == a["library_records_by_marker"] == a["map_records_by_marker"])
fa = Path(a["source_fasta"].replace("/maps/", "/"))
src_ok = sha(fa) == a["source_fasta_sha256"] if fa.exists() else None
rec = {}
for f in ("hash.k2d", "opts.k2d", "taxo.k2d", "taxid-to-path.tsv", "seqid2taxid.map", "audit.json", "build-commands.sh"):
    rec[f] = {"bytes": (db / f).stat().st_size, "sha256": sha(db / f)}
# cross-check audit recorded hashes of derived files
cross = {k: rec[k]["sha256"] == v for k, v in a["outputs"].items() if k in rec}
insp = subprocess.run(["kraken2-inspect", "--db", str(db)], capture_output=True, text=True)
lines = insp.stdout.splitlines()
m = {"created_utc": datetime.now(timezone.utc).isoformat(), "db": name, "path": str(db),
     "audit_counts_equal": ok, "audit_source_fasta_sha256_recomputed_matches": src_ok,
     "audit_recorded_hashes_match": cross, "audit": a, "files": rec,
     "kraken2_inspect_returncode": insp.returncode, "kraken2_inspect_head": lines[:12],
     "kraken2_inspect_stderr": insp.stderr.strip().splitlines()[:5]}
(db / "db-manifest.json").open("x").write(json.dumps(m, indent=2, sort_keys=True) + "\n")
with open(db / "db-manifest.SHA256SUMS", "x") as h:
    for f in list(rec) + ["db-manifest.json"]:
        h.write(f"{sha(db / f)}  {f}\n")
for f in ("db-manifest.json", "db-manifest.SHA256SUMS"):
    os.chmod(db / f, 0o444)
print(name, "counts_equal", ok, "src_sha_ok", src_ok, "cross", cross, "inspect_rc", insp.returncode)
print("\n".join(lines[:8]))
print({k: v for k, v in rec.items()})
