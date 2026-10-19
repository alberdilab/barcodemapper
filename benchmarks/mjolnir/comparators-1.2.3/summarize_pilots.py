"""Combine per-run evaluate_calls.py outputs into mock_tables_*.tsv, summary.tsv and summary.md."""
import csv, json, sys
from collections import Counter
from pathlib import Path

bench = Path(sys.argv[1])
P = bench / "projects/03-method-comparison"
R = P / "results/pilots-1.2.3"
E = R / "evaluation"
RUNS = [("SRR14351851", "fungal pooled DNA, all-fungi reference"), ("ERR6001161", "animal mock, default reference"),
        ("ERR6001173", "extraction blank (negative control), default reference")]
METHODS = ["barcodemapper_standard", "barcodemapper_allhits", "best_hit", "vsearch", "minimap2_preset", "minimap2_sensitive"]
RANKS = ["species", "genus", "family"]


def tsv(path):
    with open(path, newline="") as h:
        return list(csv.DictReader(h, delimiter="\t"))


for name in ("x",):
    pass
for kind, fname in (("metrics", "mock_metrics.tsv"), ("detection", "expected_detection.tsv"), ("unexpected", "unexpected_calls.tsv")):
    rows = []
    for run, _ in RUNS:
        for mf in (1, 2):
            if kind != "metrics" and mf == 2:
                continue
            for r in tsv(E / f"{run}.minfrag{mf}" / fname):
                rows.append({"run": run, **r})
    fields = list(rows[0])
    with open(E / f"mock_tables_{kind}.tsv", "x", newline="") as h:
        w = csv.DictWriter(h, fields, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)

calls = {}
for run, _ in RUNS:
    for m in METHODS:
        calls[run, m] = {r["ReadID"]: r for r in tsv(E / "calls" / f"{run}.{m}.tsv")}

summary = []
for run, _ in RUNS:
    for m in METHODS:
        c = calls[run, m]; n = len(c); st = Counter(r["status"] for r in c.values())
        sp = sum(1 for r in c.values() if ";s__" in r["Taxonomy"] and not r["Taxonomy"].endswith("s__"))
        met = {(r["rank"], r["min_fragments"]): r for r in tsv(E / f"{run}.minfrag1" / "mock_metrics.tsv") if r["method"] == m}
        met2 = {r["rank"]: r for r in tsv(E / f"{run}.minfrag2" / "mock_metrics.tsv") if r["method"] == m}
        unexp = [r for r in tsv(E / f"{run}.minfrag1" / "unexpected_calls.tsv") if r["method"] == m]
        row = {"run": run, "method": m, "pairs": n, "pairs_with_hit": n - st["no_hit"], "classified": st["classified"],
               "unclassified": st["unclassified"], "classified_fraction": f"{st['classified'] / n:.6f}", "species_resolved_pairs": sp}
        for rank in RANKS:
            a = met[(rank, "1")]
            row[f"expected_{rank}"] = a["expected_taxa"]; row[f"detected_{rank}_mf1"] = a["detected_taxa"]
            row[f"unexpected_{rank}_mf1"] = a["unexpected_taxa"]
            row[f"detected_{rank}_mf2"] = met2[rank]["detected_taxa"]; row[f"unexpected_{rank}_mf2"] = met2[rank]["unexpected_taxa"]
            row[f"unexpected_{rank}_fragments"] = sum(int(r["fragments"]) for r in unexp if r["rank"] == rank)
        summary.append(row)
with open(E / "summary.tsv", "x", newline="") as h:
    w = csv.DictWriter(h, list(summary[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(summary)

# timings
timing = []
for d, label in (("vsearch", "vsearch"), ("minimap2", "minimap2 preset"), ("minimap2-sensitive", "minimap2 sensitive")):
    for run, _ in RUNS:
        mf = json.load(open(R / d / f"{run}.manifest.json"))["timing"]
        for mate in ("R1", "R2"):
            t = mf[mate]
            timing.append((label, run, mate, t["wall_seconds"], t["cpu_hours"], t["peak_rss_kib"] / 2**20, t["node"].split(".")[0], t["exit_status"]))
with open(E / "timings.tsv", "x") as h:
    h.write("method\trun\tmate\twall_seconds\tcpu_hours\tpeak_rss_GiB\tnode\texit_status\n")
    for t in timing:
        h.write("\t".join(f"{x:.3f}" if isinstance(x, float) else str(x) for x in t) + "\n")

# consistency checks
checks = []
for run, _ in RUNS:
    hit = {m: {k for k, r in calls[run, m].items() if r["status"] != "no_hit"} for m in METHODS}
    checks.append(f"{run}: best_hit hit pairs = {len(hit['best_hit'])}, barcodemapper_standard = {len(hit['barcodemapper_standard'])}, "
                  f"identical hit sets: {hit['best_hit'] == hit['barcodemapper_standard']}; allhits-only pairs {len(hit['barcodemapper_allhits'] - hit['barcodemapper_standard'])}, "
                  f"standard-only pairs {len(hit['barcodemapper_standard'] - hit['barcodemapper_allhits'])}")
    checks.append(f"{run}: pairs hit by vsearch and not by BarcodeMapper standard: {len(hit['vsearch'] - hit['barcodemapper_standard'])}; "
                  f"BarcodeMapper standard hit and not vsearch: {len(hit['barcodemapper_standard'] - hit['vsearch'])}; "
                  f"minimap2 preset hit and not standard: {len(hit['minimap2_preset'] - hit['barcodemapper_standard'])}; "
                  f"minimap2 sensitive hit and not standard: {len(hit['minimap2_sensitive'] - hit['barcodemapper_standard'])}")
    pre, sen = hit["minimap2_preset"], hit["minimap2_sensitive"]
    checks.append(f"{run}: minimap2 preset hit pairs {len(pre)}, sensitive {len(sen)}, preset not in sensitive {len(pre - sen)}")

md = ["# Pilots 1.2.3: comparator evaluation on real mocks (exploratory)", "",
      "BarcodeMapper 1.2.3 `--standard` (-k 200), `--all-hits`, single best-hit ablation, vsearch, minimap2 preset (primary minimap2 row) and minimap2 sensitive variant (`-f 20000,100000`, separately labelled). "
      "Truth: Project 01 `expected-taxa.tsv` (exact paper binomials, synonyms unreviewed, so recall is a lower bound and synonym names count as unexpected). "
      "Matching is by name at species, genus and family, regardless of kingdom (so off-target kingdoms count as unexpected). Fragments = read pairs. "
      "Counts are one pilot per marker, 92k to 99k pairs, so they are indicative only. Per-fragment synthetic truth (Project 02) is not built; synthetic mode was not run. "
      "**Pending: Kraken2 (database builds queued) and BLASTn (subset awaiting decisions).**", ""]
for run, desc in RUNS:
    md += [f"## {run}: {desc}", ""]
    rows = [r for r in summary if r["run"] == run]
    md += ["Detected / expected (min 1 fragment; in parentheses min 2 fragments); unexpected = taxa called that are not expected, with total fragments.", "",
           "| method | pairs with hit | classified | classified fraction | species detected | genus detected | family detected | unexpected species (fragments) | unexpected genera | unexpected families |",
           "|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        md.append(f"| {r['method']} | {r['pairs_with_hit']} | {r['classified']} | {float(r['classified_fraction']):.5f} | "
                  + " | ".join(f"{r[f'detected_{k}_mf1']}/{r[f'expected_{k}']} ({r[f'detected_{k}_mf2']})" for k in RANKS)
                  + f" | {r['unexpected_species_mf1']} ({r['unexpected_species_fragments']}) | {r['unexpected_genus_mf1']} | {r['unexpected_family_mf1']} |")
    md.append("")
    unexp = tsv(E / f"{run}.minfrag1" / "unexpected_calls.tsv")
    md.append("Top unexpected species per method (fragments):"); md.append("")
    for m in METHODS:
        top = sorted((r for r in unexp if r["method"] == m and r["rank"] == "species"), key=lambda r: -int(r["fragments"]))[:5]
        md.append(f"- {m}: " + ("; ".join(f"{r['taxon']} ({r['fragments']})" for r in top) if top else "none"))
    md.append("")
md += ["## Exploratory timings (not the final timing runs)", "",
       "16 CPUs, classifier-only (index/database load included), `time_classifier.sh`; nodes vary (including a GPU node), so wall times are not strictly comparable. BarcodeMapper, best-hit and Kraken2 timings are not part of this table.", "",
       "| method | run | mate | wall s | CPU-hours | peak RSS GiB | node |", "|---|---|---|---|---|---|---|"]
for t in timing:
    md.append(f"| {t[0]} | {t[1]} | {t[2]} | {t[3]:.0f} | {t[4]:.2f} | {t[5]:.1f} | {t[6]} |")
md += ["", "## Consistency checks (for converter or evaluator bugs)", ""] + [f"- {c}" for c in checks]
md += ["", "## Notes", "",
       "- Blank ERR6001173 has no expected taxa; evaluate_calls.py needs at least one row, so a placeholder species `NO_EXPECTED_TAXA_BLANK_CONTROL` was used. All calls are therefore unexpected.",
       "- vsearch uses `--maxrejects 1000` (not exhaustive) and the exact BarcodeMapper acceptance rule applied post hoc; minimap2 rows use the same rule.",
       "- All conversions here use the committed (HEAD) converter; `script_and_truth_hashes.txt` lists hashes."]
(E / "summary.md").write_text("\n".join(md) + "\n")
print("\n".join(checks))
