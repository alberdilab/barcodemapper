#!/usr/bin/env python3
"""Evaluate the BarcodeMapper Project 01 pilot (mock communities), locally.

Reads only files under inputs/ and results/ next to this script and writes
tables, plots and README.md to analysis/. No network, no HPC.

Rerun:  python3 benchmarks/pilot-01/evaluate_pilot.py
"""
import argparse
import glob
import json
import os
import re
import sys
from collections import defaultdict

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
THRESHOLDS = (1, 2, 5)
RANKS = ("species", "genus", "family")
PREFIX = {"k": "kingdom", "p": "phylum", "c": "class", "o": "order",
          "f": "family", "g": "genus", "s": "species"}
RANK_ORDER = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
# Target kingdom of the mock, taken from sample-roles.tsv `role`.
ROLE_KINGDOM = {"animal_mock": "Animalia", "fungal_pooled_DNA": "Fungi"}
BLANK_ROLES = {"extraction_blank"}
MARKER_LABEL = {"animal_mock": "COI", "fungal_pooled_DNA": "ITS"}


def norm(name):
    """Case/underscore/space-insensitive key for an exact name match."""
    return re.sub(r"\s+", " ", str(name).replace("_", " ")).strip().lower()


# ---------------------------------------------------------------- parsing
def parse_taxonomy(path, run):
    """Return DataFrame(path, ranks dict parts, rank, name, kingdom, cum, excl)
    plus the Unclassified count. Counts are cumulative (verified in checks)."""
    df = pd.read_csv(path, sep="\t", dtype={"Taxonomy": str})
    col = run if run in df.columns else df.columns[1]
    df = df.rename(columns={col: "count"})
    unclassified = int(df.loc[df.Taxonomy == "Unclassified", "count"].sum())
    df = df[df.Taxonomy != "Unclassified"].copy()
    df["parts"] = df.Taxonomy.str.split(";")
    df["rank"] = df.parts.map(lambda p: PREFIX[p[-1][0]])
    df["name"] = df.parts.map(lambda p: p[-1].split("__", 1)[1])
    df["kingdom"] = df.parts.map(lambda p: p[0].split("__", 1)[1])
    df["depth"] = df.parts.map(len)
    cum = dict(zip(df.Taxonomy, df["count"]))
    child_sum = defaultdict(int)
    n_children = defaultdict(int)
    for tax, p in zip(df.Taxonomy, df.parts):
        if len(p) > 1:
            parent = ";".join(p[:-1])
            child_sum[parent] += cum[tax]
            n_children[parent] += 1
    df["child_sum"] = df.Taxonomy.map(lambda t: child_sum.get(t, 0))
    df["n_children"] = df.Taxonomy.map(lambda t: n_children.get(t, 0))
    df["excl"] = df["count"] - df["child_sum"]
    return df, unclassified


def build_expected(exp_run):
    """Expected sets per rank from the expected-taxa rows of one run.

    Returns dict rank -> DataFrame(key, label, in_exact_species_db,
    in_rank_db) and a list of (species, reason) that could not be derived."""
    rows = {r: {} for r in RANKS}
    underivable = []
    for _, r in exp_run.iterrows():
        sp = r.expected_species
        exact = r.coverage_status == "exact_species"
        items = [("species", sp, exact, exact)]
        genus = str(sp).split()[0] if str(sp).strip() else ""
        if genus:
            items.append(("genus", genus, exact, r.genus_reference_records > 0))
        else:
            underivable.append((sp, "genus"))
        fam = r.source_families
        if isinstance(fam, str) and fam.strip():
            items.append(("family", fam.strip(), exact, r.family_reference_records > 0))
        else:
            underivable.append((sp, "family"))
        for rank, label, ex, rk in items:
            d = rows[rank].setdefault(norm(label), dict(
                key=norm(label), label=label, in_exact_species_db=False,
                in_rank_db=False))
            d["in_exact_species_db"] |= bool(ex)
            d["in_rank_db"] |= bool(rk)
    return {r: pd.DataFrame(list(v.values())) for r, v in rows.items()}, underivable


def safe_div(a, b):
    return a / b if b else float("nan")


# ---------------------------------------------------------------- metrics
def rank_metrics(run, meta, tax, expected, kingdom_scopes):
    out = []
    for rank in RANKS:
        e = expected[rank]
        e_all = set(e.key)
        e_exact = set(e.key[e.in_exact_species_db])
        e_rank = set(e.key[e.in_rank_db])
        for scope in kingdom_scopes:
            sub = tax[tax["rank"] == rank]
            if scope != "all":
                sub = sub[sub.kingdom == scope]
            for t in THRESHOLDS:
                det = set(sub.loc[sub["count"] >= t, "name"].map(norm))
                tp = det & e_all
                out.append(dict(
                    run=run, marker=meta["marker_label"], reference=meta["reference"],
                    rank=rank, kingdom_scope=scope, min_fragments=t,
                    detected=len(det), expected_all=len(e_all),
                    expected_exact_species_db=len(e_exact),
                    expected_rank_db=len(e_rank),
                    TP=len(tp), FP=len(det - e_all),
                    FN_all=len(e_all - det),
                    FN_exact_species_db=len(e_exact - det),
                    FN_rank_db=len(e_rank - det),
                    precision=safe_div(len(tp), len(det)),
                    recall_all_expected=safe_div(len(tp), len(e_all)),
                    recall_exact_species_db=safe_div(len(det & e_exact), len(e_exact)),
                    recall_rank_db=safe_div(len(det & e_rank), len(e_rank)),
                ))
    return out


def expected_detail(run, tax, expected):
    rows = []
    for rank in RANKS:
        counts = tax[tax["rank"] == rank].assign(k=lambda d: d.name.map(norm)) \
            .groupby("k")["count"].sum()
        for _, r in expected[rank].iterrows():
            rows.append(dict(run=run, rank=rank, taxon=r.label,
                             in_exact_species_db=r.in_exact_species_db,
                             in_rank_db=r.in_rank_db,
                             fragments=int(counts.get(r.key, 0))))
    return rows


def unexpected_calls(run, meta, tax, expected):
    exp_gen = set(expected["genus"].key)
    exp_fam = set(expected["family"].key)
    target = meta["target_kingdom"]
    rows = []
    for rank in ("species", "genus", "family", "kingdom"):
        if rank == "kingdom":
            exp = {norm(target)}
        else:
            exp = set(expected[rank].key)
        for _, r in tax[tax["rank"] == rank].iterrows():
            if norm(r["name"]) in exp:
                continue
            parts = {PREFIX[p[0]]: p.split("__", 1)[1] for p in r.parts}
            genus = parts.get("genus", "")
            if rank == "species" and not genus:
                genus = r["name"].split("_")[0]
            rows.append(dict(
                run=run, marker=meta["marker_label"], taxon=r["name"], rank=rank,
                fragments=int(r["count"]), kingdom=r.kingdom,
                cross_kingdom=r.kingdom != target,
                genus_expected=norm(genus) in exp_gen if genus else "",
                family_expected=norm(parts.get("family", "")) in exp_fam
                if parts.get("family") else "",
                placeholder_name=bool(re.search(r"_sp\.|_cf\.|_aff\.|\d", r["name"])),
                lineage=r.Taxonomy))
    return rows


def read_fate(run, meta, tax, unclassified, man):
    before = man.get("fastp_reads_before")
    after = man.get("fastp_reads_after")
    assigned = man.get("assigned_read_ids")
    excl = tax.groupby("rank")["excl"].sum()
    g = lambda *rs: int(sum(excl.get(r, 0) for r in rs))
    roots = int(tax.loc[tax.depth == 1, "count"].sum())
    return dict(
        run=run, marker=meta["marker_label"], role=meta["role"],
        reference_db_file=man.get("_db"),
        reads_before_fastp=before, reads_after_fastp=after,
        pairs_input=before // 2 if before else None,
        pairs_retained_fastp=after // 2 if after else None,
        frac_pairs_retained=safe_div(after, before),
        marker_matched_assigned_ids=assigned,
        frac_input_pairs_matched=safe_div(assigned, before // 2) if before else None,
        frac_retained_pairs_matched=safe_div(assigned, after // 2) if after else None,
        resolved_species=g("species"), resolved_genus=g("genus"),
        resolved_family=g("family"),
        resolved_higher_order_to_kingdom=g("order", "class", "phylum", "kingdom"),
        unclassified=unclassified,
        taxonomy_table_total=roots + unclassified,
        matches_assigned=(roots + unclassified) == assigned,
        frac_matched_resolved_to_species=safe_div(g("species"), assigned),
        frac_matched_resolved_to_genus_or_better=safe_div(g("species", "genus"), assigned),
    )


def data_checks(run, tax, unclassified, man):
    """Cumulative-vs-exclusive test: if cumulative, parent >= sum(children)
    for every node with children. Exclusive counts would violate that."""
    par = tax[tax.n_children > 0]
    return dict(
        run=run, nodes=len(tax), parents_with_children=len(par),
        parents_ge_children_sum=int((par["count"] >= par.child_sum).sum()),
        parents_lt_children_sum=int((par["count"] < par.child_sum).sum()),
        parents_equal_children_sum=int((par["count"] == par.child_sum).sum()),
        parents_with_own_remainder=int((par["count"] > par.child_sum).sum()),
        root_sum_plus_unclassified=int(tax.loc[tax.depth == 1, "count"].sum()) + unclassified,
        manifest_assigned_read_ids=man.get("assigned_read_ids"),
        interpretation="cumulative" if (par["count"] >= par.child_sum).all() else "CHECK",
    )


# ---------------------------------------------------------------- plots
def plot_all(out_dir, fate, metrics):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        plt = None
    paths = []
    runs = fate.run.tolist()
    labels = [f"{r}\n({'blank' if m == 'none' else m})" for r, m in zip(fate.run, fate.marker)]
    cats = [("resolved_species", "species"), ("resolved_genus", "genus"),
            ("resolved_family", "family"),
            ("resolved_higher_order_to_kingdom", "higher (order or above)"),
            ("unclassified", "Unclassified")]
    sel = metrics[(metrics.kingdom_scope == "all") & metrics.marker.isin(["COI", "ITS"])]
    if plt is None:
        paths += svg_fate(out_dir, labels, fate, cats)
        paths += svg_recall(out_dir, sel)
        return paths, False
    fig, ax = plt.subplots(figsize=(6.5, 4))
    bottom = [0] * len(runs)
    for col, lab in cats:
        ax.bar(labels, fate[col].tolist(), bottom=bottom, label=lab)
        bottom = [b + v for b, v in zip(bottom, fate[col])]
    ax.set_ylabel("marker-matched fragments (pairs)")
    ax.set_title("Read fate: matched fragments by resolution (100k input pairs)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    p = os.path.join(out_dir, "plot_read_fate.png")
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(p)
    markers = sorted(sel.marker.unique())
    fig, axes = plt.subplots(1, len(markers), figsize=(6 * len(markers), 4), squeeze=False)
    for ax, mk in zip(axes[0], markers):
        d = sel[sel.marker == mk]
        w = 0.25
        for i, rank in enumerate(RANKS):
            for j, (col, hatch) in enumerate([("recall_exact_species_db", ""),
                                              ("recall_all_expected", "//")]):
                dd = d[d["rank"] == rank].sort_values("min_fragments")
                xs = [k + i * w + j * w / 2.4 * 0 for k in range(len(THRESHOLDS))]
                ax.bar([x + (j - 0.5) * w / 2 for x in xs], dd[col], w / 2,
                       color=f"C{i}", hatch=hatch, edgecolor="white",
                       label=f"{rank} {'in-DB' if j == 0 else 'all'}")
        ax.set_xticks([k + w for k in range(len(THRESHOLDS))])
        ax.set_xticklabels([f">={t}" for t in THRESHOLDS])
        ax.set_ylim(0, 1)
        ax.set_ylabel("recall")
        ax.set_xlabel("min fragments")
        ax.set_title(f"{mk} (single run)")
        ax.legend(fontsize=7, ncol=3)
    fig.tight_layout()
    p = os.path.join(out_dir, "plot_recall.png")
    fig.savefig(p, dpi=150)
    plt.close(fig)
    paths.append(p)
    return paths, True


def _svg_write(path, w, h, body):
    with open(path, "w") as fh:
        fh.write(f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" '
                 f'font-family="sans-serif" font-size="11">'
                 f'<rect width="100%" height="100%" fill="white"/>{body}</svg>')


COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#7f7f7f", "#9467bd"]


def svg_fate(out_dir, labels, fate, cats):
    w, h, x0, y0, ph = 760, 330, 60, 40, 220
    tot = max(1, max(fate[[c for c, _ in cats]].sum(axis=1)))
    body = [f'<text x="{x0}" y="20">Matched fragments by resolution (input 100k pairs/run)</text>']
    bw = 90
    for i, (lab, (_, r)) in enumerate(zip(labels, fate.iterrows())):
        x = x0 + 20 + i * (bw + 40)
        y = y0 + ph
        for k, (col, _) in enumerate(cats):
            hh = ph * r[col] / tot
            y -= hh
            body.append(f'<rect x="{x}" y="{y:.1f}" width="{bw}" height="{hh:.1f}" fill="{COLORS[k]}"/>')
        body.append(f'<text x="{x}" y="{y0 + ph + 15}">{lab.splitlines()[0]} {lab.splitlines()[1]}</text>')
        body.append(f'<text x="{x}" y="{y - 4:.1f}">n={int(sum(r[c] for c, _ in cats))}</text>')
    for k, (_, lab) in enumerate(cats):
        body.append(f'<rect x="{w - 150}" y="{40 + 16 * k}" width="10" height="10" fill="{COLORS[k]}"/>'
                    f'<text x="{w - 135}" y="{49 + 16 * k}">{lab}</text>')
    p = os.path.join(out_dir, "plot_read_fate.svg")
    _svg_write(p, w, h, "".join(body))
    return [p]


def svg_recall(out_dir, sel):
    markers = sorted(sel.marker.unique())
    w, h, x0, y0, pw, ph = 380 * len(markers) + 40, 320, 50, 40, 300, 200
    body = []
    for m, mk in enumerate(markers):
        ox = x0 + m * 380
        body.append(f'<text x="{ox}" y="22">{mk} recall (solid: in-DB denominator, '
                    f'light: all expected)</text>')
        body.append(f'<line x1="{ox}" y1="{y0 + ph}" x2="{ox + pw}" y2="{y0 + ph}" stroke="black"/>')
        body.append(f'<line x1="{ox}" y1="{y0}" x2="{ox}" y2="{y0 + ph}" stroke="black"/>')
        for tick in (0, 0.5, 1):
            body.append(f'<text x="{ox - 28}" y="{y0 + ph - ph * tick + 4:.1f}">{tick:.1f}</text>')
        d = sel[sel.marker == mk]
        bw = 11
        for ti, t in enumerate(THRESHOLDS):
            gx = ox + 10 + ti * 100
            body.append(f'<text x="{gx + 20}" y="{y0 + ph + 15}">&gt;={t}</text>')
            for ri, rank in enumerate(RANKS):
                row = d[(d["rank"] == rank) & (d.min_fragments == t)]
                if row.empty:
                    continue
                row = row.iloc[0]
                for j, col in enumerate(["recall_exact_species_db", "recall_all_expected"]):
                    v = 0 if pd.isna(row[col]) else row[col]
                    x = gx + ri * 2 * bw + j * bw
                    body.append(f'<rect x="{x}" y="{y0 + ph - ph * v:.1f}" width="{bw}" '
                                f'height="{ph * v:.1f}" fill="{COLORS[ri]}" '
                                f'opacity="{1 if j == 0 else 0.45}"/>')
        for k, rank in enumerate(RANKS):
            body.append(f'<rect x="{ox}" y="{y0 + ph + 30 + 14 * k}" width="10" height="10" fill="{COLORS[k]}"/>'
                        f'<text x="{ox + 15}" y="{y0 + ph + 39 + 14 * k}">{rank}</text>')
    p = os.path.join(out_dir, "plot_recall.svg")
    _svg_write(p, w, h + 40, "".join(body))
    return [p]


# ---------------------------------------------------------------- README
def fmt(x, nd=2):
    return "NA" if x is None or (isinstance(x, float) and pd.isna(x)) else f"{x:.{nd}f}"


def write_readme(out_dir, metrics, fate, checks, unexp, blank, underiv, plots, png_ok, runs_meta):
    meta_k = {r: m["target_kingdom"] for r, m in runs_meta.items()}
    L = []
    L.append("# BarcodeMapper Project 01 pilot: local evaluation\n")
    L.append("**Scope and caveats.** Preliminary 3-run pilot: one animal mock (COI, run ERR6001161), "
             "one fungal pooled-DNA sample (ITS, run SRR14351851) and one extraction blank "
             "(ERR6001173), each subsampled to 100k read pairs (200k reads) and run with the "
             "'pilot-standard' BarcodeMapper settings. There is exactly one run per marker, so "
             "precision and recall are indicative only; no confidence intervals or "
             "between-replicate variation are available. Species recall is a lower bound: "
             "`synonym_status` is `unreviewed` for every expected taxon (exact paper binomial "
             "only, no synonyms invented), so a taxon reported under a current synonym is "
             "counted as a miss and its reported name as an unexpected call. The blank role "
             "(ERR6001173) is an inference in `sample-roles.tsv`.\n")
    L.append("## Method\n")
    L.append("- Counts in the taxonomy TSVs are **cumulative** (a taxon includes all its "
             "descendants); see `data_checks.tsv`. Exclusive counts are derived as "
             "parent minus sum of children.\n"
             "- A taxon is *detected* at threshold t if its cumulative fragment count is >= t "
             "(t = 1, 2, 5). `Unclassified` is never a taxon.\n"
             "- Matching: exact name, case-insensitive, `_` equals space. Expected genus = first "
             "word of the expected binomial; expected family = `source_families`.\n"
             "- Precision = TP / detected at that rank (all detected taxa, including other "
             "kingdoms; a `kingdom_scope` row restricted to the target kingdom is also given).\n"
             "- Recall has three denominators: `recall_all_expected` (all expected taxa at the "
             "rank), `recall_exact_species_db` (expected taxa derived from species with "
             "`coverage_status == exact_species`, i.e. the species is in the reference) and "
             "`recall_rank_db` (expected taxa with >0 reference records at that rank, "
             "from the genus/family record counts).\n")
    L.append("## Headline numbers (kingdom_scope = all)\n")
    L.append("| marker | rank | min frag | detected | expected (all / in-DB sp.) | TP | FP | precision | recall all | recall in-DB |")
    L.append("|---|---|---|---|---|---|---|---|---|---|")
    m = metrics[metrics.kingdom_scope == "all"]
    for _, r in m.iterrows():
        L.append(f"| {r.marker} | {r['rank']} | {r.min_fragments} | {r.detected} | "
                 f"{r.expected_all} / {r.expected_exact_species_db} | {r.TP} | {r.FP} | "
                 f"{fmt(r.precision)} | {fmt(r.recall_all_expected)} | "
                 f"{fmt(r.recall_exact_species_db)} |")
    L.append("")
    L.append("## Read fate\n")
    L.append("Pairs retained by fastp use `fastp_reads_after / 2` (assumes paired filtering). "
             "'Matched' is `assigned_read_ids` from the manifest, which equals the sum of "
             "kingdom-level counts plus `Unclassified` (checked), i.e. matched fragments include "
             "those left `Unclassified`; that reading of `Unclassified` is inferred, not "
             "documented in the files.\n")
    cols = ["run", "marker", "pairs_input", "pairs_retained_fastp", "marker_matched_assigned_ids",
            "frac_input_pairs_matched", "resolved_species", "resolved_genus", "resolved_family",
            "resolved_higher_order_to_kingdom", "unclassified"]
    L.append("| " + " | ".join(cols) + " |")
    L.append("|" + "---|" * len(cols))
    for _, r in fate.iterrows():
        L.append("| " + " | ".join(
            fmt(r[c], 5) if c.startswith("frac") else str(r[c]) for c in cols) + " |")
    L.append("")
    L.append("Fields available in the inputs: manifest `fastp_reads_before`, `fastp_reads_after`, "
             "`assigned_read_ids`; `*.results.json` holds only artifact paths/hashes (no counts). "
             "Per-fragment assignments and the pre-match-filter counts are not available locally.\n")
    L.append("## Surprises worth a look\n")
    L.append("- Very few pairs reach an accepted marker match (21, 226 and 4 of 100,000 pairs for the "
             "animal, fungal and blank runs); 100k pairs of shotgun data hold only a handful of "
             "marker fragments, so the species-level evidence per expected taxon is 0-3 fragments "
             "and precision/recall are dominated by sampling noise.\n"
             "- The fungal run was mapped with the `all-fungi` reference yet received "
             "Animalia calls (26 fragments, mostly Cecidomyiidae, 21) and a Dipodidae call; "
             "check what non-fungal ITS records the all-fungi database contains.\n"
             "- The blank (default reference) has Viridiplantae (Trebouxiophyceae, 1 fragment), "
             "Hymenoptera (2) and a Sminthurididae (1) fragment although it is a blank; "
             "Sminthurididae is also an expected family of the animal mock "
             "(*Sphaeridia pumilis*), so sample-to-blank carry-over cannot be excluded.\n"
             "- The animal mock run received one Fungi/Polyporaceae fragment.\n")
    L.append("## Unexpected calls\n")
    if len(unexp):
        for run, d in unexp.groupby("run"):
            s = d[d["rank"] == "species"]
            L.append(f"- {run}: {len(s)} unexpected species, {int(s.fragments.sum())} fragments; "
                     f"{int(d[(d['rank']=='kingdom')].fragments.sum())} fragments assigned to a "
                     f"kingdom other than the target ({meta_k[run]}). "
                     f"Full list: `unexpected_calls.tsv`.")
    L.append("")
    L.append("Names that look like synonyms of expected taxa (flagged for manual review, **not** "
             "applied to any metric): `Scedosporium_prolificans` vs expected "
             "`Lomentospora prolificans`; `Nakaseomyces glabratus` vs expected `Candida glabrata`. "
             "`Eupelops_sp._SFH2024_01` is a placeholder species name in a genus whose expected "
             "species is `Eupelops occultus`. This is why species recall is a lower bound.\n")
    L.append("## Blank control (ERR6001173)\n")
    if len(blank):
        top = blank[blank.depth == 1]
        L.append("Kingdom-level calls: " + ", ".join(f"{r['name']}={r['count']}" for _, r in top.iterrows())
                 + ". Full lineage list in `blank_calls.tsv`. No precision/recall computed; "
                 "calls are not interpreted as contamination without further evidence. "
                 "Family-level calls also present in the animal-mock expected set are flagged in "
                 "`in_animal_mock_expected_family`.\n")
    L.append("## Expected taxa without a derivable rank\n")
    if underiv:
        for run, items in underiv.items():
            fam = sorted({s for s, r in items if r == "family"})
            gen = sorted({s for s, r in items if r == "genus"})
            L.append(f"- {run}: genus not derivable for {len(gen)} taxa; family not derivable "
                     f"(empty `source_families`) for {len(fam)} taxa: {', '.join(fam) or 'none'}. "
                     f"These are absent from the family-level expected set.")
    L.append("")
    L.append("## Files\n")
    for f in ["rank_metrics.tsv", "expected_detection.tsv", "unexpected_calls.tsv",
              "read_fate.tsv", "blank_calls.tsv", "data_checks.tsv"] + [os.path.basename(p) for p in plots]:
        L.append(f"- `{f}`")
    if not png_ok:
        L.append("\nmatplotlib is not installed in this environment, so plots were written as "
                 "SVG fallbacks; install matplotlib and rerun to get PNGs.")
    L.append("")
    with open(os.path.join(out_dir, "README.md"), "w") as fh:
        fh.write("\n".join(L))


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--inputs", default=os.path.join(HERE, "inputs"))
    ap.add_argument("--results", default=os.path.join(HERE, "results"))
    ap.add_argument("--out", default=os.path.join(HERE, "analysis"))
    ap.add_argument("--suffix", default="-pilot-standard")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    exp = pd.read_csv(os.path.join(a.inputs, "expected-taxa.tsv"), sep="\t")
    roles = pd.read_csv(os.path.join(a.inputs, "sample-roles.tsv"), sep="\t").set_index("run_accession")

    metrics, detail, unexp, fate, checks, blank, underiv = [], [], [], [], [], [], {}
    runs_meta = {}
    all_expected_keys = set()
    for tsv in sorted(glob.glob(os.path.join(a.results, f"*{a.suffix}.tsv"))):
        run = os.path.basename(tsv)[: -len(a.suffix + ".tsv")]
        role = roles.loc[run, "role"] if run in roles.index else "unknown"
        with open(os.path.join(a.results, f"{run}{a.suffix}.manifest.json")) as fh:
            mj = json.load(fh)
        man = dict(mj["samples"][run])
        man["_db"] = os.path.basename(mj["database"]["path"])
        tax, uncl = parse_taxonomy(tsv, run)
        e_run = exp[exp.run_accession == run]
        marker_label = MARKER_LABEL.get(role, "none")
        meta = dict(role=role, marker_label=marker_label,
                    reference=e_run.reference.iloc[0] if len(e_run) else
                    os.path.basename(mj["database"]["path"]).split(".")[0],
                    target_kingdom=ROLE_KINGDOM.get(role, ""))
        runs_meta[run] = meta
        fate.append(read_fate(run, meta, tax, uncl, man))
        checks.append(data_checks(run, tax, uncl, man))
        if role in BLANK_ROLES or e_run.empty:
            b = tax.copy()
            b.insert(0, "run", run)
            b["unclassified_run_total"] = uncl
            blank.append(b[["run", "Taxonomy", "rank", "name", "kingdom", "depth", "count",
                            "excl", "unclassified_run_total"]])
            continue
        expected, und = build_expected(e_run)
        underiv[run] = und
        all_expected_keys.update((rk, k) for rk in RANKS for k in expected[rk].key)
        scopes = ["all", meta["target_kingdom"]]
        metrics += rank_metrics(run, meta, tax, expected, scopes)
        detail += expected_detail(run, tax, expected)
        unexp += unexpected_calls(run, meta, tax, expected)

    if blank and underiv:
        pass
    metrics = pd.DataFrame(metrics)
    fate = pd.DataFrame(fate)
    checks = pd.DataFrame(checks)
    unexp = pd.DataFrame(unexp).sort_values(["run", "rank", "fragments"], ascending=[True, True, False])
    blank = pd.concat(blank) if blank else pd.DataFrame()
    if len(blank):
        blank = blank.sort_values(["run", "Taxonomy"])
        blank["in_animal_mock_expected_family"] = [
            (rk == "family" and ("family", norm(n)) in all_expected_keys)
            if rk == "family" else "" for rk, n in zip(blank["rank"], blank["name"])]
    metrics.to_csv(os.path.join(a.out, "rank_metrics.tsv"), sep="\t", index=False, float_format="%.4f")
    pd.DataFrame(detail).to_csv(os.path.join(a.out, "expected_detection.tsv"), sep="\t", index=False)
    unexp.to_csv(os.path.join(a.out, "unexpected_calls.tsv"), sep="\t", index=False)
    fate.to_csv(os.path.join(a.out, "read_fate.tsv"), sep="\t", index=False, float_format="%.6f")
    blank.to_csv(os.path.join(a.out, "blank_calls.tsv"), sep="\t", index=False)
    checks.to_csv(os.path.join(a.out, "data_checks.tsv"), sep="\t", index=False)

    # plots: only runs with a marker (blank included in read-fate plot)
    for old in glob.glob(os.path.join(a.out, "plot_*")):
        os.remove(old)
    plots, png_ok = plot_all(a.out, fate, metrics)
    write_readme(a.out, metrics, fate, checks, unexp, blank, underiv, plots, png_ok, runs_meta)
    print(f"wrote outputs to {a.out}")
    with pd.option_context("display.width", 250, "display.max_columns", 40):
        print(metrics[metrics.kingdom_scope == "all"][
            ["marker", "rank", "min_fragments", "detected", "expected_all",
             "expected_exact_species_db", "TP", "FP", "precision",
             "recall_all_expected", "recall_exact_species_db"]].to_string(index=False))
        print(fate.T.to_string())
        print(checks.T.to_string())


if __name__ == "__main__":
    sys.exit(main())
