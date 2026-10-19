"""v3 summary: all methods incl. Kraken2; invalid vsearch only as appendix row."""
import csv, json, re, subprocess, sys
from collections import Counter
from pathlib import Path

bench = Path(sys.argv[1])
P = bench / "projects/03-method-comparison"; R = P / "results/pilots-1.2.3"
E = R / "evaluation/v3"; C0 = R / "evaluation/calls"
RUNS = [("SRR14351851", "fungal pooled DNA, all-fungi reference"), ("ERR6001161", "animal mock, default reference"),
        ("ERR6001173", "extraction blank (negative control), default reference")]
INVALID = "vsearch (invalid, N artifact)"
MAIN = ["barcodemapper_standard", "barcodemapper_allhits", "best_hit", "vsearch (n_mismatch)", "minimap2_preset", "minimap2_sensitive",
        "Kraken2 (conf 0.0)", "Kraken2 (conf 0.1)"]
METHODS = MAIN + [INVALID]
PATH = lambda run, m: {"Kraken2 (conf 0.0)": R / f"kraken2/{run}.calls.tsv", "Kraken2 (conf 0.1)": R / f"kraken2-conf0.1/{run}.calls.tsv",
                       "vsearch (n_mismatch)": C0 / f"{run}.vsearch_n_mismatch.tsv", INVALID: C0 / f"{run}.vsearch_invalid.tsv"}.get(m, C0 / f"{run}.{m}.tsv")
RANKS = ["species", "genus", "family"]


def tsv(path):
    with open(path, newline="") as h:
        return list(csv.DictReader(h, delimiter="\t"))


for kind, fname in (("metrics", "mock_metrics.tsv"), ("detection", "expected_detection.tsv"), ("unexpected", "unexpected_calls.tsv")):
    rows = [{"run": run, "in_main_table": "no (appendix)" if r["method"] == INVALID else "yes", **r}
            for run, _ in RUNS for mf in ((1, 2) if kind == "metrics" else (1,)) for r in tsv(E / f"{run}.minfrag{mf}" / fname)]
    with open(E / f"mock_tables_{kind}.tsv", "x", newline="") as h:
        w = csv.DictWriter(h, list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)

calls = {(run, m): {r["ReadID"]: r for r in tsv(PATH(run, m))} for run, _ in RUNS for m in METHODS}
summary = []
for run, _ in RUNS:
    m1 = tsv(E / f"{run}.minfrag1" / "mock_metrics.tsv"); m2 = tsv(E / f"{run}.minfrag2" / "mock_metrics.tsv")
    unexp = tsv(E / f"{run}.minfrag1" / "unexpected_calls.tsv")
    for m in METHODS:
        c = calls[run, m]; n = len(c); st = Counter(r["status"] for r in c.values())
        row = {"run": run, "method": m, "in_main_table": "no (appendix)" if m == INVALID else "yes", "pairs": n, "pairs_with_hit": st["classified"] if m.startswith("Kraken2") else n - st["no_hit"],
               "classified": st["classified"], "unclassified": st["unclassified"], "classified_fraction": f"{st['classified'] / n:.6f}"}
        for rank in RANKS:
            a = next(r for r in m1 if r["method"] == m and r["rank"] == rank); b = next(r for r in m2 if r["method"] == m and r["rank"] == rank)
            row.update({f"expected_{rank}": a["expected_taxa"], f"detected_{rank}_mf1": a["detected_taxa"], f"unexpected_{rank}_mf1": a["unexpected_taxa"],
                        f"detected_{rank}_mf2": b["detected_taxa"], f"unexpected_{rank}_mf2": b["unexpected_taxa"],
                        f"unexpected_{rank}_fragments": sum(int(r["fragments"]) for r in unexp if r["method"] == m and r["rank"] == rank)})
        summary.append(row)
with open(E / "summary.tsv", "x", newline="") as h:
    w = csv.DictWriter(h, list(summary[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(summary)

# classifier-only timings from manifests
timing = []
for d, label in (("vsearch-nmismatch", "vsearch (n_mismatch)"), ("minimap2", "minimap2 preset"), ("minimap2-sensitive", "minimap2 sensitive"),
                 ("kraken2", "Kraken2 (conf 0.0)"), ("kraken2-conf0.1", "Kraken2 (conf 0.1)")):
    for run, _ in RUNS:
        mf = json.load(open(R / d / f"{run}.manifest.json"))["timing"]
        recs = [mf[k] for k in ("R1", "R2")] if "R1" in mf else [mf]
        for k, t in zip(("R1", "R2") if "R1" in mf else ("pair",), recs):
            timing.append((label, run, k, t["wall_seconds"], t["cpu_hours"], t["peak_rss_kib"] / 2**20, t["node"].split(".")[0], t["cpus"]))
with open(E / "timings.tsv", "x") as h:
    h.write("method\trun\tunit\twall_seconds\tcpu_hours\tpeak_rss_GiB\tnode\tcpus\n")
    for t in timing: h.write("\t".join(f"{x:.3f}" if isinstance(x, float) else str(x) for x in t) + "\n")

# BarcodeMapper pilots: sacct rows from tasks/1.2.3-pilots/summary.txt (mapping worker = longest non-controller job)
def hms(s):
    d = 0
    if "-" in s: d, s = s.split("-"); d = int(d)
    p = [float(x) for x in s.split(":")]; p = [0] * (3 - len(p)) + p
    return d * 86400 + p[0] * 3600 + p[1] * 60 + p[2]
sections, cur = {}, None
for line in open(bench / "tasks/1.2.3-pilots/summary.txt"):
    m = re.match(r"--- (\S+-pilot-\w+): jobs", line)
    if m: cur = m[1]; sections[cur] = []; continue
    if line.startswith("---") or line.startswith("==="): cur = None; continue
    if cur and "|" in line and not line.startswith(" "): sections[cur].append(line.rstrip("\n").split("|"))
bm = []
for sec, rows in sections.items():
    main = [r for r in rows if "." not in r[0]]
    worker = max(main[1:], key=lambda r: hms(r[4]))
    step = next((r for r in rows if r[0] == worker[0] + ".0"), None)
    bm.append((sec, worker[0], worker[4], step[6] if step else worker[6], worker[8], main[0][4]))
with open(E / "timings_barcodemapper_sacct.tsv", "x") as h:
    h.write("pilot\tmapping_worker_job\telapsed\tMaxRSS_step\talloc_cpus\tcontroller_elapsed\n")
    for r in bm: h.write("\t".join(r) + "\n")
bh = []
for run, tag, jid in (("SRR14351851", "fungal", "45939557"), ("ERR6001161", "animal", "45939562"), ("ERR6001173", "blank", "45939567")):
    out = subprocess.run(["sacct", "-j", jid, "-n", "-P", "-o", "JobID,State,Elapsed,MaxRSS,AllocCPUS"], capture_output=True, text=True).stdout
    rows = [l.split("|") for l in out.splitlines()]
    top = next(r for r in rows if r[0] == jid); b = next((r for r in rows if r[0] == jid + ".batch"), top)
    bh.append((run, jid, top[1], top[2], b[3], top[4]))
with open(E / "timings_besthit_sacct.tsv", "x") as h:
    h.write("run\tjob\tstate\telapsed\tMaxRSS\talloc_cpus\n")
    for r in bh: h.write("\t".join(r) + "\n")

# consistency checks
checks = []
for run, _ in RUNS:
    hit = {m: {k for k, r in calls[run, m].items() if (r["status"] == "classified" if m.startswith("Kraken2") else r["status"] != "no_hit")} for m in METHODS}
    s = hit["barcodemapper_standard"]
    checks.append(f"{run}: best_hit hit set == BarcodeMapper standard: {hit['best_hit'] == s}; vsearch (n_mismatch) hit-only {len(hit['vsearch (n_mismatch)'] - s)}, BarcodeMapper-only {len(s - hit['vsearch (n_mismatch)'])}; "
                  f"minimap2 preset hit-only {len(hit['minimap2_preset'] - s)}, BarcodeMapper-only {len(s - hit['minimap2_preset'])}; "
                  f"Kraken2 0.0 hit-only {len(hit['Kraken2 (conf 0.0)'] - s)}, BarcodeMapper-only {len(s - hit['Kraken2 (conf 0.0)'])}; Kraken2 0.1 subset of 0.0: {hit['Kraken2 (conf 0.1)'] <= hit['Kraken2 (conf 0.0)']}")

# ---- markdown
def table(run, methods):
    rows = {r["method"]: r for r in summary if r["run"] == run}
    out = ["| method | pairs with hit | classified | classified fraction | species detected | genus detected | family detected | unexpected species (fragments) | unexpected genera | unexpected families |",
           "|---|---|---|---|---|---|---|---|---|---|"]
    for m in methods:
        r = rows[m]
        out.append(f"| {m} | {r['pairs_with_hit']} | {r['classified']} | {float(r['classified_fraction']):.5f} | "
                   + " | ".join(f"{r[f'detected_{k}_mf1']}/{r[f'expected_{k}']} ({r[f'detected_{k}_mf2']})" for k in RANKS)
                   + f" | {r['unexpected_species_mf1']} ({r['unexpected_species_fragments']}) | {r['unexpected_genus_mf1']} | {r['unexpected_family_mf1']} |")
    return out

md = ["# Pilots 1.2.3: comparator evaluation on real mocks, v3 (exploratory)", "",
      "Methods: BarcodeMapper 1.2.3 `--standard` (-k 200), `--all-hits`, single best-hit ablation, vsearch (n_mismatch rerun), minimap2 preset (primary minimap2 row), minimap2 sensitive (`-f 20000,100000`, labelled variant), Kraken2 (confidence 0.0, primary) and Kraken2 (confidence 0.1, predeclared variant). "
      "The first vsearch run is invalid (N-run artifact) and appears only in the appendix. BLASTn is still pending. Supersedes summary.md (v1) and v2.",
      "", "Truth: Project 01 `expected-taxa.tsv` (exact paper binomials). Fragments = read pairs. Table cells: detected/expected taxa at min 1 fragment (min 2 fragments in parentheses). "
      "Unexpected = called taxa not expected, with total fragments in parentheses. `classified fraction` = classified pairs / retained pairs. `pairs with hit`: pairs not `no_hit` (for Kraken2, which has no no_hit status, classified pairs). FASTQ and reference sha256 were verified for every method in each run (identical across methods).", ""]
for run, desc in RUNS:
    md += [f"## {run}: {desc}", ""] + table(run, MAIN) + [""]
    unexp = tsv(E / f"{run}.minfrag1" / "unexpected_calls.tsv")
    md.append("Top unexpected species (fragments):"); md.append("")
    for m in MAIN:
        top = sorted((r for r in unexp if r["method"] == m and r["rank"] == "species"), key=lambda r: -int(r["fragments"]))[:4]
        md.append(f"- {m}: " + ("; ".join(f"{r['taxon']} ({r['fragments']})" for r in top) if top else "none"))
    md.append("")
md += ["## Exploratory timings (not the final timing runs)", "",
       "Classifier-only wall time, CPU-hours and peak RSS from `time_classifier.sh` (`timing.json`), 16 CPUs, on mixed node classes (compute nodes comp01 to comp16 and GPU node gpu03 of the same partition, so wall times are not strictly comparable). "
       "Index or database load is included. vsearch and minimap2 ran R1 and R2 separately; Kraken2 ran the pair in one command.", "",
       "| method | run | unit | wall s | CPU-hours | peak RSS GiB | node |", "|---|---|---|---|---|---|---|"]
md += [f"| {t[0]} | {t[1]} | {t[2]} | {t[3]:.0f} | {t[4]:.3f} | {t[5]:.1f} | {t[6]} |" for t in timing]
md += ["", "BarcodeMapper and best-hit (sacct of the pilot jobs, **not comparable to classifier-only timings**: they include Snakemake, 8-CPU Bowtie2 mapping with a cap of 200 hits or `--all-hits`, and extraction; best-hit is a post-processing script on the BarcodeMapper BAM-derived extract, 1 CPU):", "",
       "| pilot | mapping worker job | elapsed | MaxRSS (step) | alloc CPUs | controller elapsed |", "|---|---|---|---|---|---|"]
md += [f"| {r[0]} | {r[1]} | {r[2]} | {r[3]} | {r[4]} | {r[5]} |" for r in bm]
md += ["", "| best-hit run | job | state | elapsed | MaxRSS | alloc CPUs |", "|---|---|---|---|---|---|"] + [f"| {r[0]} | {r[1]} | {r[2]} | {r[3]} | {r[4]} | {r[5]} |" for r in bh]
md += ["", "## Consistency checks", ""] + [f"- {c}" for c in checks]
md += ["", "## Exploratory notes (read-only, login node)", "", "(a) Why does vsearch (n_mismatch) hit 446 fungal pairs against 227 for BarcodeMapper? SRR14351851: 220 pairs are hit by vsearch only and 1 by BarcodeMapper only. For 10 evenly spaced vsearch-only pairs (full listing: exploratory_notes.txt), the vsearch alignments are partial: 81 to 214 aligned bases, 27 to 71% of the 300 bp read, at 97 to 100% identity. BarcodeMapper maps with Bowtie2 in end-to-end mode (no --local), so the whole read has to align. In 4 of the 10 pairs Bowtie2 reported no alignment at all; in the other 6 it aligned the full read (about 300 bases) to other references with NM 11 to 35 (3.7 to 11.7% of covered bases), which fails the 3% rule. So the extra vsearch hits come from local, partial alignments of 300 bp reads that overhang the marker; the shared absolute 80-base threshold admits them, and this is not a -k 200, soft-clipping or threshold effect in Bowtie2. The same logic probably applies to minimap2 (local alignments); this was not checked for minimap2 or the animal mock.",
       "",
       "(b) Why does Kraken2 classify few pairs? Pairs with at least one 35-mer present in the database (a non-zero, non-ambiguous entry in the k-mer hit string of output.txt), against pairs classified: SRR14351851 13,818 (15.0%) vs 7,151 (7.7%); ERR6001161 9,146 (9.3%) vs 757 (0.8%); ERR6001173 8,247 (8.3%) vs 1,162 (1.2%). Only 0.74%, 0.17% and 0.13% of read k-mers hit the database. So most pairs have no marker k-mer (shotgun reads are mostly off-target, as for BarcodeMapper), and many pairs with a k-mer hit stay unclassified; the command does not override Kraken2's default --minimum-hit-groups 2, which is the probable cause (not verified). At confidence 0.1 classified pairs fall to 1,096, 96 and 26. Kraken2 is not low relative to BarcodeMapper in the animal mock and blank (757 and 1,162 classified vs 21 and 4 pairs hit); those extra calls are counted in the tables above.", "",
       "## Caveats", "",
       "- One 100k-pair subsample per mock (92k to 99k pairs after fastp); counts of marker pairs are tiny (a few to a few hundred), so detection numbers are indicative only; no replicates or confidence intervals.",
       "- Expected-taxa matching uses exact paper binomials; species-name synonyms (for example Nakaseomyces glabratus for Candida glabrata) count as unexpected, so recall is a lower bound.",
       "- Matching ignores kingdom: off-target-kingdom calls count as unexpected. The `all-fungi` reference also contains Animalia records.",
       "- vsearch is not exhaustive (`--maxrejects 1000`, `--maxaccepts 0`, top-identity ties only) and its alignments are global with free terminal gaps; the acceptance rule (covered >= 80 bases, NM <= 3%) is applied post hoc.",
       "- minimap2 returns many more hit pairs than BarcodeMapper on the animal mock (78 vs 21) without a mechanism being established here.",
       "- Kraken2 classifies few pairs; it needs at least 2 hit groups by default (see note b), and the database hash is only about 0.4 GB.",
       "- BLASTn is pending. Synthetic per-fragment truth (Project 02) is not built; synthetic mode was not run.",
       "- Blank ERR6001173 uses a placeholder expected species so that every call is unexpected.",
       "", "## Appendix: invalid vsearch run (N-run artifact), not part of the comparison", "",
       "The first vsearch run scored N as match; 89 to 95% of its accepted alignments lie on reference stretches that are at least 50% N (see evaluation/v2/n_run_diagnostic.tsv). Retained for the record only.", ""]
for run, desc in RUNS:
    md += [f"{run}:", ""] + table(run, [INVALID]) + [""]
(E / "summary.md").write_text("\n".join(md) + "\n")
print("\n".join(checks))
