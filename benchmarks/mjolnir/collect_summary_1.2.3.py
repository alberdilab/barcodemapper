"""Read-only summary of the 1.2.3 pilots and comparison with the 1.1.2 pilots."""
import json, subprocess, sys
from pathlib import Path
import yaml

root = Path("/projects/alberdilab/scratch/jpl786/projects/barcodemapper-benchmarks")
p01 = root / "projects/01-mock-communities"
p03 = root / "projects/03-method-comparison"
res = p03 / "results/pilots-1.2.3"

def tax(path):
    d = {}
    for line in path.read_text().splitlines()[1:]:
        r, t = line.split("\t")[:2]
        d[r] = t
    return d

def table(path):
    d = {}
    for line in path.read_text().splitlines()[1:]:
        t, n = line.split("\t")[:2]
        d[t] = int(n)
    return d

def sacct(jobs):
    out = subprocess.run(["sacct", "-n", "-P", "-j", ",".join(jobs),
        "--format=JobID,JobName%30,State,ExitCode,Elapsed,Timelimit,MaxRSS,ReqMem,AllocCPUS"],
        capture_output=True, text=True).stdout
    return out

for run, label in (("SRR14351851", "fungal"), ("ERR6001161", "animal"), ("ERR6001173", "blank")):
    print("=" * 30, run, label)
    old = tax(p01 / f"work/{run}-pilot-standard/taxonomy/{run}.tsv")
    oldt = table(p01 / f"results/{run}-pilot-standard.tsv")
    for mode in ("standard", "allhits"):
        lab = f"{run}-pilot-{mode}"
        rec = yaml.safe_load((res / f"{lab}.run.yaml").read_text())
        resj = json.loads((res / f"{lab}.results.json").read_text())
        print(f"--- {lab}: jobs {rec['slurm_job_ids']}, fastp retained pairs {resj['fastp_retained_pairs']}")
        for a in resj["artifacts"]:
            print(f"  sha256 {a['sha256']}  {a['bytes']:>12}  {a['path'].split('/projects/03-method-comparison/')[-1]}")
        print(sacct(rec["slurm_job_ids"]))
        new = tax(p03 / f"work/pilots-1.2.3/{lab}/taxonomy/{run}.tsv")
        newt = table(res / f"{lab}.tsv")
        only_new = set(new) - set(old); only_old = set(old) - set(new)
        chg = [r for r in set(new) & set(old) if new[r] != old[r]]
        print(f"  classified fragments: 1.1.2={len(old)} 1.2.3={len(new)}; only-in-1.2.3={len(only_new)} only-in-1.1.2={len(only_old)} changed-taxon={len(chg)}")
        print(f"  Unclassified in table: 1.1.2={oldt.get('Unclassified')} 1.2.3={newt.get('Unclassified')}; table rows 1.1.2={len(oldt)} 1.2.3={len(newt)}")
        diffs = sorted(((newt.get(t, 0) - oldt.get(t, 0), t) for t in set(oldt) | set(newt) if newt.get(t, 0) != oldt.get(t, 0)))
        for d, t in diffs[:12] + (diffs[-12:] if len(diffs) > 24 else diffs[12:]):
            print(f"    {d:+d}  old={oldt.get(t,0)} new={newt.get(t,0)}  {t}")
        for r in list(only_new)[:5]: print("    new-only", r, new[r])
        for r in list(only_old)[:5]: print("    old-only", r, old[r])
        for r in chg[:5]: print("    changed", r, old[r], "->", new[r])
    bh = json.loads((res / f"best-hit/{run}.manifest.json").read_text())
    print("--- best-hit:", {k: bh[k] for k in ("input_pairs_after_fastp","pairs_with_accepted_alignment","pairs_with_conflicting_taxonomy_score_tie","output_sha256")})
    old_bh = p03 / f"results/pilot-best-hit/{run}.manifest.json"
    if old_bh.exists():
        o = json.loads(old_bh.read_text()); print("    1.1.2 best-hit:", o["pairs_with_accepted_alignment"], "of", o["input_pairs_after_fastp"], "ties", o["pairs_with_conflicting_taxonomy_score_tie"])
        a = {l.split("\t")[0]: l.split("\t")[1] for l in (res / f"best-hit/{run}.tsv").read_text().splitlines()[1:]}
        b = {l.split("\t")[0]: l.split("\t")[1] for l in old_bh.with_suffix("").with_suffix(".tsv").read_text().splitlines()[1:]}
        print("    best-hit taxon differences vs 1.1.2:", sum(a.get(k) != b.get(k) for k in set(a) | set(b)))
