import argparse
import base64
import colorsys
import json
from datetime import UTC, datetime
from functools import cache
from html import escape
from pathlib import Path

import pandas as pd
from plotly.offline import get_plotlyjs

try:
    from barcodemapper.workflow.scripts.level_stats import (
        load_merged,
        load_stats,
        overview_stats,
    )
except ModuleNotFoundError:
    from level_stats import load_merged, load_stats, overview_stats


RANK_NAMES = ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
LOGO_PATH = Path(__file__).resolve().parents[1] / "assets" / "logo.png"


def _logo_uri():
    try:
        return "data:image/png;base64," + base64.b64encode(LOGO_PATH.read_bytes()).decode("ascii")
    except OSError:
        return ""


def _version_label(version):
    version = str(version)
    return version if version in ("", "unknown") else "v" + version


LEVEL_PALETTES = {
    "light": ("#1d568f", "#407cbd", "#70a5d6", "#a9cbea"),
    "dark": ("#a9cefa", "#80b4ed", "#5396dc", "#2b72b7"),
}


def parse_args():
    parser = argparse.ArgumentParser(description="Generate a static HTML sunburst plot with sidebar sample selection.")
    parser.add_argument("--input", required=True, action="append", help="NAME=PATH for each level, or PATH for one legacy level")
    parser.add_argument("--stats", action="append", default=[], help="Per-sample extraction stats JSON")
    parser.add_argument("--limits", action="append", default=[], metavar="NAME=MAXPCT:MINCOV",
                        help="Stringency thresholds for an overview level")
    parser.add_argument("--output", required=True, help="Output HTML file")
    parser.add_argument("--version", default="unknown", help="Version of BarcodeMapper")
    return parser.parse_args()


@cache
def assign_color(kingdom, depth):
    if kingdom.startswith("k__Animalia"):
        base_hue = 211
    elif kingdom.startswith("k__Fungi"):
        base_hue = 17
    elif kingdom.startswith("k__Viridiplantae"):
        base_hue = 160
    else:
        base_hue = 0

    lightness = min(0.42 + depth * 0.055, 0.76)
    saturation = 0.65 if kingdom.startswith("k__") else 0
    rgb = colorsys.hls_to_rgb(base_hue / 360, lightness, saturation)
    red, green, blue = (int(channel * 255) for channel in rgb)
    return f"#{red:02x}{green:02x}{blue:02x}"


def prepare_data(df):
    samples = df.columns[1:]
    all_data = {}
    all_phylum_rows = []
    paths = [tuple(level for level in row if pd.notnull(level))
             for row in df["Taxonomy"].str.split(";", expand=True).itertuples(index=False, name=None)]

    for sample in samples:
        nodes = {}
        for abundance, levels in zip(df[sample].to_numpy(), paths):
            if abundance <= 0:
                continue
            full_path = []
            for depth, label in enumerate(levels):
                full_path.append(label)
                path_id = ";".join(full_path)
                parent_id = ";".join(full_path[:-1]) if depth > 0 else ""
                if path_id not in nodes:
                    nodes[path_id] = {
                        "id": path_id,
                        "label": label,
                        "parent": parent_id,
                        "depth": depth,
                        "value": 0,
                        "kingdom": levels[0]
                    }
                if depth == len(levels) - 1:
                    nodes[path_id]["value"] += abundance

        frame = pd.DataFrame(nodes.values(), columns=["id", "label", "parent", "depth", "value", "kingdom"])
        frame["color"] = [assign_color(kingdom, depth)
                          for kingdom, depth in zip(frame["kingdom"], frame["depth"])]
        all_data[sample] = frame.to_dict(orient="list")

        if frame.empty:
            phylum_rows = pd.DataFrame({
                "Phylum": pd.Series(dtype="string"),
                "Kingdom": pd.Series(dtype="string"),
                sample: pd.Series(dtype="int64"),
                "color": pd.Series(dtype="string"),
            })
        else:
            phylum_rows = frame[(frame["depth"] == 1) | (frame["label"] == "Unclassified")][["label", "kingdom", "value", "color"]].copy()
            phylum_rows.columns = ["Phylum", "Kingdom", sample, "color"]
            phylum_rows = phylum_rows.reset_index(drop=True)
            child_totals = frame.groupby("parent")["value"].sum()
            kingdoms = frame[frame["depth"] == 0]
            for label, path_id, value, color in zip(kingdoms["label"], kingdoms["id"],
                                                    kingdoms["value"], kingdoms["color"]):
                if label == "Unclassified":
                    continue
                unresolved = value - child_totals.get(path_id, 0)
                if unresolved > 0:
                    phylum_rows.loc[len(phylum_rows)] = [
                        "Unresolved within " + label,
                        label, unresolved, color,
                    ]
        all_phylum_rows.append(phylum_rows)

    merged = all_phylum_rows[0] if all_phylum_rows else pd.DataFrame(columns=["Phylum", "Kingdom", "color"])
    for phylum_df in all_phylum_rows[1:]:
        merged = pd.merge(merged, phylum_df, on=["Phylum", "Kingdom", "color"], how="outer")
    merged = merged.fillna(0)
    merged = merged.sort_values(["Kingdom", "Phylum"])

    rel = merged.copy()
    for sample in samples:
        total = rel[sample].sum()
        rel[sample] = rel[sample] / total if total > 0 else 0

    return list(samples), all_data, merged, rel




def safe_json(value):
    return json.dumps(value).replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")


def _percent(value):
    return f"{value * 100:.1f}%" if value is not None else "n/a"


def _compact_level(df):
    """Store a DFS taxonomy with sparse cumulative counts per taxon and sample."""
    samples = list(df.columns[1:])
    paths = df["Taxonomy"].fillna("").astype(str).tolist()
    node_index = {}
    labels, parents = [], []
    terminal_indices = []
    for path in paths:
        parent = -1
        for label in path.split(";"):
            if not label:
                continue
            key = (parent, label)
            if key not in node_index:
                node_index[key] = len(labels)
                labels.append(label)
                parents.append(parent)
            parent = node_index[key]
        terminal_indices.append(parent)

    children = [[] for _ in labels]
    roots = []
    for index, parent in enumerate(parents):
        (roots if parent < 0 else children[parent]).append(index)
    kingdom_order = {"k__Animalia": 0, "k__Fungi": 1, "k__Viridiplantae": 2}
    roots.sort(key=lambda index: (kingdom_order.get(labels[index], 3), labels[index]))
    for siblings in children:
        siblings.sort(key=lambda index: labels[index])
    ordered, dfs_parents, depths = [], [], []

    def visit(index, parent, depth):
        ordered.append(index)
        dfs_parents.append(parent)
        depths.append(depth)
        position = len(ordered) - 1
        for child in children[index]:
            visit(child, position, depth + 1)

    for root in roots:
        visit(root, -1, 0)
    old_to_new = {old: new for new, old in enumerate(ordered)}
    labels = [labels[index] for index in ordered]
    parents = dfs_parents
    children = [[] for _ in labels]
    for index, parent in enumerate(parents):
        if parent >= 0:
            children[parent].append(index)

    values = df.iloc[:, 1:].fillna(0).to_numpy(dtype="int64", copy=False)
    node_values = [{} for _ in labels]
    explicit = set()
    for row, old_index in enumerate(terminal_indices):
        if old_index < 0:
            continue
        index = old_to_new[old_index]
        explicit.add(index)
        for col, value in enumerate(values[row]):
            if value > 0:
                node_values[index][col] = node_values[index].get(col, 0) + int(value)
    for index in range(len(labels) - 1, -1, -1):
        if index not in explicit:
            for child in children[index]:
                for sample_index, value in node_values[child].items():
                    node_values[index][sample_index] = node_values[index].get(sample_index, 0) + value

    totals, classified = [], []
    raw_composition = []
    for col in range(len(samples)):
        root_counts = {index: node_values[index].get(col, 0)
                       for index, parent in enumerate(parents) if parent < 0}
        total = sum(root_counts.values())
        totals.append(total)
        classified.append(sum(count for index, count in root_counts.items() if labels[index] != "Unclassified"))
        phyla = {}
        for index, depth in enumerate(depths):
            if depth == 1 and node_values[index].get(col, 0) > 0:
                root_index = parents[index]
                count = node_values[index].get(col, 0)
                key = (labels[root_index], labels[index])
                phyla[key] = phyla.get(key, 0) + count
        for index, count in root_counts.items():
            if labels[index] == "Unclassified" and count > 0:
                phyla[("Unclassified", "Unclassified")] = count
            else:
                unresolved = count - sum(node_values[child].get(col, 0) for child in children[index])
                if unresolved > 0:
                    key = (labels[index], "Unresolved within " + labels[index])
                    phyla[key] = unresolved
        raw_composition.append(phyla)

    by_kingdom = {}
    for composition in raw_composition:
        for (kingdom, phylum), value in composition.items():
            if phylum == "Unclassified" or phylum.startswith("Unresolved within "):
                continue
            by_kingdom.setdefault(kingdom, {})[phylum] = by_kingdom.setdefault(kingdom, {}).get(phylum, 0) + value
    categories = []
    for kingdom in sorted(by_kingdom, key=lambda item: (kingdom_order.get(item, 3), item)):
        ranked = sorted(by_kingdom[kingdom], key=lambda item: (-by_kingdom[kingdom][item], item))
        categories.extend((kingdom, phylum) for phylum in ranked[:8])
        if any((kingdom, "Unresolved within " + kingdom) in item for item in raw_composition):
            categories.append((kingdom, "Unresolved within " + kingdom))
        if len(ranked) > 8:
            categories.append((kingdom, "Other " + kingdom.removeprefix("k__")))
    if any(("Unclassified", "Unclassified") in item for item in raw_composition):
        categories.append(("Unclassified", "Unclassified"))
    composition_data = []
    for kingdom, phylum in categories:
        if phylum.startswith("Other "):
            selected = {name for k, name in categories if k == kingdom and not name.startswith("Other ")
                        and not name.startswith("Unresolved within ")}
            series = [sum(value for (k, name), value in row.items()
                          if k == kingdom and name not in selected and not name.startswith("Unresolved within "))
                      for row in raw_composition]
        else:
            series = [row.get((kingdom, phylum), 0) for row in raw_composition]
        composition_data.append({"kingdom": kingdom, "name": phylum, "values": series})
    sparse = [[[sample_index, count] for sample_index, count in sorted(counts.items())]
              for counts in node_values]
    return {"samples": samples, "taxa": {"labels": labels, "parents": parents,
            "depths": depths, "values": sparse}, "totals": totals,
            "classified": classified, "composition": composition_data}


def _number(value):
    return f"{value:,}" if value is not None else "n/a"


def _limit_text(limit):
    if limit is None:
        return ""
    return f"≤{limit[0]:g}% edits · ≥{limit[1]} bp"


def _overview_table(overview, limits):
    headings = ["Level", "Thresholds", "Fragments retained", "% of aligned fragments",
                "Classified fragments", "% classified"] + [f"Distinct {rank.lower()} taxa" for rank in RANK_NAMES]
    rows = []
    max_retained = max((item["totals"]["fragments_retained"] or 0 for item in overview["levels"].values()), default=0)
    for index, (name, item) in enumerate(overview["levels"].items()):
        totals = item["totals"]
        classified = totals["classified_fragments"]
        denominator = classified + totals["unclassified_fragments"]
        fraction = classified / denominator if denominator else None
        retained = totals["fragments_retained"]
        cells = [f'<th scope="row">{escape(name)}</th>',
                 f'<td>{escape(_limit_text(limits.get(name)) or "n/a")}</td>']
        for value, proportion, formatted in (
            (retained, retained / max_retained if retained is not None and max_retained else 0, _number(retained)),
            (totals["retained_fraction"], totals["retained_fraction"] or 0, _percent(totals["retained_fraction"])),
            (classified, classified / denominator if denominator else 0, _number(classified)),
            (fraction, fraction or 0, _percent(fraction)),
        ):
            cells.append(f'<td class="bar-cell"><span class="cell-bar" style="width:{max(0, min(proportion, 1)) * 100:.1f}%"></span><span>{escape(formatted)}</span></td>')
        cells.extend(f'<td>{count:,}</td>' for count in totals["distinct_taxa"])
        rows.append(f'<tr data-level-index="{index}">' + "".join(cells) + "</tr>")
    return ('<div class="table-wrap"><table id="comparison-table"><thead><tr>'
            + "".join(f'<th scope="col">{escape(heading)}</th>' for heading in headings)
            + '</tr></thead><tbody>' + "".join(rows) + '</tbody></table></div>')


STYLE = r'''
:root { color-scheme: light; --page:#f7f6f3; --surface:#fff; --text:#202624; --muted:#66716d; --border:#e5e8e4; --accent:#2a78d6; --soft:#eaf2fd; --grid:#e8ece8; --bar-strip:rgba(42,120,214,.11); --shadow:0 3px 16px rgba(32,38,36,.035); --neutral:#9ba5a1; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) { color-scheme:dark; --page:#141a1b; --surface:#1d2527; --text:#edf1ee; --muted:#a1b0aa; --border:#344140; --accent:#3987e5; --soft:#1c344b; --grid:#35413f; --bar-strip:rgba(57,135,229,.19); --shadow:0 4px 20px rgba(0,0,0,.13); --neutral:#889792; } }
:root[data-theme="dark"] { color-scheme:dark; --page:#141a1b; --surface:#1d2527; --text:#edf1ee; --muted:#a1b0aa; --border:#344140; --accent:#3987e5; --soft:#1c344b; --grid:#35413f; --bar-strip:rgba(57,135,229,.19); --shadow:0 4px 20px rgba(0,0,0,.13); --neutral:#889792; }
* { box-sizing:border-box; }
body { margin:0; background:var(--page); color:var(--text); font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Inter,Helvetica,Arial,sans-serif; }
button,input,select { font:inherit; }
button { cursor:pointer; }
button:focus-visible,input:focus-visible,select:focus-visible,a:focus-visible,summary:focus-visible { outline:2px solid var(--accent); outline-offset:2px; }
.topbar { display:flex; align-items:center; justify-content:space-between; flex-wrap:wrap; gap:12px; padding:15px calc(max(0px,(100vw - 1600px)/2) + 24px); background:var(--surface); border-bottom:1px solid var(--border); }
.brand { display:flex; align-items:flex-end; gap:11px; }
.brand a { color:var(--text); text-decoration:none; font-size:18px; font-weight:700; letter-spacing:-.04em; display:inline-flex; align-items:center; }
.brand img { height:44px; width:auto; display:block; }
:root[data-theme="dark"] .brand img { background:#fff; border-radius:5px; padding:3px 6px; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) .brand img { background:#fff; border-radius:5px; padding:3px 6px; } }
.brand a:hover { color:var(--accent); }
.brand span,.run-meta { color:var(--muted); font-size:12px; }
.run-meta { display:flex; align-items:center; flex-wrap:wrap; gap:14px; font-variant-numeric:tabular-nums; }
.icon-button,.small-button,.segment button { border:1px solid var(--border); background:var(--surface); color:var(--text); border-radius:8px; padding:7px 11px; }
.icon-button:hover,.small-button:hover,.segment button:hover { border-color:var(--accent); }
.shell { max-width:1600px; padding:0 24px 48px; margin:auto; }
#level-tabs { display:flex; gap:2px; overflow-x:auto; border-bottom:1px solid var(--border); margin-bottom:25px; }
.level-tab { border:0; border-bottom:2px solid transparent; background:transparent; color:var(--muted); padding:14px 18px 12px; white-space:nowrap; text-align:left; }
.level-tab strong { display:block; color:var(--text); font-size:13px; font-weight:600; }
.level-tab small { display:block; font-size:11px; }
.level-tab[aria-selected="true"] { border-color:var(--accent); color:var(--accent); }
.level-tab[aria-selected="true"] strong { color:var(--accent); }
.level-tab[data-level]:not([data-level="overview"])[aria-selected="true"] { border-color:var(--level-accent); }
.level-tab[data-level]:not([data-level="overview"])[aria-selected="true"] strong { color:var(--level-accent); }
.page-head { margin:28px 0 20px; }
h1,h2,h3,p { margin-top:0; }
h1 { font-size:28px; letter-spacing:-.045em; line-height:1.2; margin-bottom:6px; }
h2 { font-size:16px; letter-spacing:-.015em; margin-bottom:2px; }
.eyebrow { color:var(--muted); font-size:10px; font-weight:700; text-transform:uppercase; letter-spacing:.13em; margin-bottom:7px; }
.subtle { color:var(--muted); font-size:12px; }
.card { min-width:0; background:var(--surface); border:1px solid var(--border); border-radius:11px; box-shadow:var(--shadow); padding:20px; }
.card-head { display:flex; align-items:center; justify-content:space-between; flex-wrap:wrap; gap:10px; margin-bottom:16px; }
.kpis { display:grid; grid-template-columns:repeat(auto-fit,minmax(160px,1fr)); gap:12px; margin-bottom:18px; }
.kpi { padding:17px 19px; }
.kpi[data-level-index] { border-top:2px solid var(--level-accent); }
.kpi .eyebrow { white-space:nowrap; overflow:hidden; text-overflow:ellipsis; }
.kpi-value { display:block; font-size:26px; font-weight:650; letter-spacing:-.05em; line-height:1.2; font-variant-numeric:tabular-nums; }
.kpi-note { color:var(--muted); font-size:11px; }
.overview-grid { display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:16px; }
.span-all { grid-column:1/-1; }
.chart { height:340px; width:100%; }
#classified-chart { min-height:340px; height:auto; }
.table-wrap { width:100%; overflow:auto; }
table { border-collapse:collapse; min-width:100%; font-size:12px; font-variant-numeric:tabular-nums; }
th,td { padding:9px 11px; border-bottom:1px solid var(--border); white-space:nowrap; text-align:right; }
th { color:var(--muted); font-size:10px; letter-spacing:.045em; font-weight:700; }
th:first-child,td:first-child { text-align:left; position:sticky; left:0; background:var(--surface); }
tbody tr:hover td { background:var(--page); }
#comparison-table tbody tr[data-level-index] th:first-child { border-left:3px solid var(--level-accent); }
.bar-cell { position:relative; min-width:100px; }
.bar-cell span:last-child { position:relative; z-index:1; }
.cell-bar { position:absolute; right:6px; top:50%; transform:translateY(-50%); height:62%; max-width:calc(100% - 12px); background:var(--bar-strip); border-radius:3px; }
.toolbar { display:flex; align-items:center; flex-wrap:wrap; gap:9px; }
input[type="search"],select { background:var(--surface); border:1px solid var(--border); border-radius:7px; color:var(--text); padding:7px 10px; min-width:0; }
input[type="search"] { width:220px; }
.table-sort { border:0; background:none; color:inherit; padding:0; font-size:inherit; font-weight:inherit; letter-spacing:inherit; }
.table-sort:hover { color:var(--accent); }
#sample-table-count,.navigator-count { color:var(--muted); font-size:11px; }
.note { margin:18px 0; border-left:2px solid var(--border); padding:6px 0 6px 13px; color:var(--muted); font-size:12px; max-width:950px; }
#detail-panel { display:block; }
.sample-section { display:grid; grid-template-columns:230px minmax(0,1fr); gap:18px; align-items:start; }
.navigator { position:sticky; top:12px; padding:15px; }
.navigator-head { display:flex; align-items:center; justify-content:space-between; }
.navigator-head h2 { margin:0; }
.nav-buttons { display:flex; gap:4px; }
.nav-buttons button { padding:2px 8px; }
#sample-search { width:100%; margin:12px 0 8px; }
#sample-buttons { max-height:min(64vh,700px); overflow-y:auto; margin:10px -6px 0; padding:0 4px; }
.sample-button { display:block; border:0; border-radius:6px; background:transparent; color:var(--text); padding:7px 8px; margin:1px 0; width:100%; text-align:left; }
.sample-button:hover { background:var(--page); }
.sample-button.active { background:var(--soft); color:var(--accent); font-weight:600; }
.sample-name { display:block; overflow:hidden; white-space:nowrap; text-overflow:ellipsis; font-size:12px; }
.sample-meter { display:block; height:2px; border-radius:2px; background:var(--border); margin-top:5px; }
.sample-meter span { display:block; height:100%; border-radius:2px; background:var(--accent); }
.detail-main { min-width:0; display:grid; gap:16px; }
.segment { display:inline-flex; gap:2px; background:var(--page); padding:3px; border-radius:8px; }
.segment button { border:0; background:transparent; font-size:11px; padding:5px 9px; }
.segment button[aria-pressed="true"] { background:var(--surface); box-shadow:0 1px 3px rgba(0,0,0,.11); color:var(--accent); }
.plot-scroll { overflow-x:auto; }
#composition-legend { display:flex; align-items:flex-start; gap:14px; max-height:112px; overflow:auto; padding:3px 1px 10px; margin-bottom:4px; border-bottom:1px solid var(--border); }
.legend-group { flex:0 0 auto; min-width:160px; max-width:430px; }
.legend-heading { display:flex; align-items:center; gap:5px; color:var(--muted); font-size:10px; font-weight:700; text-transform:uppercase; letter-spacing:.08em; margin-bottom:4px; }
.legend-items { display:flex; flex-wrap:wrap; gap:3px; }
.legend-chip { display:inline-flex; align-items:center; gap:5px; max-width:145px; border:1px solid var(--border); background:var(--surface); color:var(--text); border-radius:5px; padding:2px 5px; font-size:10px; line-height:1.3; }
.legend-chip:hover { border-color:var(--accent); }
.legend-chip[aria-pressed="false"] { opacity:.48; text-decoration:line-through; }
.legend-chip-label { overflow:hidden; text-overflow:ellipsis; white-space:nowrap; }
.legend-swatch { flex:0 0 8px; width:8px; height:8px; border-radius:2px; background:var(--swatch); }
#barplot { min-height:560px; }
#sunburst { width:100%; height:480px; }
.taxon-toolbar { display:flex; align-items:center; flex-wrap:wrap; gap:8px; margin-bottom:12px; }
.taxon-toolbar input[type="search"] { width:180px; }
.taxon-grid { --tree-width:360px; max-height:70vh; min-height:180px; overflow:auto; border:1px solid var(--border); border-radius:8px; position:relative; }
.taxon-grid-inner { width:max-content; min-width:100%; }
.taxon-grid-head,.taxon-row { display:grid; align-items:center; }
.taxon-grid-head { position:sticky; top:0; z-index:3; height:96px; background:var(--surface); border-bottom:1px solid var(--border); }
.taxon-row { height:30px; border-bottom:1px solid var(--border); }
.taxon-row:hover { background:var(--page); }
.taxon-tree-cell { position:sticky; left:0; z-index:2; align-self:stretch; display:flex; align-items:center; width:var(--tree-width); min-width:var(--tree-width); max-width:var(--tree-width); overflow:hidden; background:var(--surface); border-right:1px solid var(--border); }
.taxon-row:hover .taxon-tree-cell { background:var(--page); }
.taxon-grid-head .taxon-tree-cell { z-index:4; font-size:11px; font-weight:700; color:var(--muted); padding-left:10px; }
.taxon-node { position:relative; display:flex; align-items:center; gap:4px; width:100%; min-width:0; padding:3px 5px 3px calc(5px + var(--depth)*10px); border:0; color:var(--text); background:transparent; text-align:left; font-size:11px; }
.taxon-node::before { content:""; width:9px; flex:none; border-top:1px solid var(--border); }
.taxon-node::after { content:""; position:absolute; left:calc(5px + var(--depth)*10px); top:0; height:50%; border-left:1px solid var(--border); }
.taxon-resize { position:absolute; right:0; top:0; bottom:0; width:10px; cursor:col-resize; touch-action:none; background:linear-gradient(90deg,transparent 40%,var(--border) 40%,var(--border) 60%,transparent 60%); }
.taxon-resize:hover,.taxon-resize:focus-visible { background:linear-gradient(90deg,transparent 40%,var(--accent) 40%,var(--accent) 60%,transparent 60%); outline:2px solid var(--accent); outline-offset:-2px; }
.taxon-node-label { overflow:hidden; white-space:nowrap; text-overflow:ellipsis; }
.taxon-rank { color:var(--muted); font-size:9px; min-width:12px; }
.taxon-toggle { flex:none; width:18px; height:18px; padding:0; border:1px solid var(--border); border-radius:4px; background:var(--surface); color:var(--muted); font-size:12px; line-height:1; }
.taxon-sample-head { position:relative; height:96px; width:22px; border:0; padding:0; background:transparent; color:var(--muted); }
.taxon-sample-head span { position:absolute; bottom:5px; left:12px; display:block; width:85px; overflow:hidden; white-space:nowrap; text-overflow:ellipsis; text-align:left; transform:rotate(-58deg); transform-origin:0 100%; font-size:10px; }
.taxon-sample-head:hover { color:var(--accent); }
.taxon-sample-head.selected,.taxon-cell.selected { background:var(--soft); }
.taxon-sample-head.selected { color:var(--accent); font-weight:700; }
.taxon-cell { display:flex; align-items:center; justify-content:center; width:22px; height:29px; }
.taxon-bubble { display:block; border-radius:50%; background:var(--accent); }
.taxon-heat { display:block; width:18px; height:23px; border-radius:2px; }
.taxon-legend { color:var(--muted); font-size:11px; margin:8px 0 0; }
.taxon-legend-swatch { display:inline-block; width:55px; height:8px; border-radius:3px; vertical-align:middle; background:linear-gradient(90deg,transparent,var(--accent)); margin:0 6px; }
.taxon-selection-note { display:flex; align-items:center; flex-wrap:wrap; gap:8px; margin-top:8px; color:var(--muted); font-size:11px; }
.taxon-selection-note button { border:0; padding:0; background:none; color:var(--accent); text-decoration:underline; }
.taxon-spacer { height:0; }
.sample-metrics { display:flex; gap:18px; color:var(--muted); font-size:11px; font-variant-numeric:tabular-nums; }
.sample-metrics strong { color:var(--text); font-size:13px; }
details { margin-top:12px; }
summary { color:var(--accent); cursor:pointer; font-size:12px; }
.empty { color:var(--muted); padding:55px 16px; text-align:center; }
[hidden] { display:none !important; }
@media (max-width:850px) { .overview-grid { grid-template-columns:1fr; } .sample-section { grid-template-columns:1fr; } .navigator { position:static; } #sample-buttons { display:flex; max-height:140px; overflow:auto; } .sample-button { min-width:130px; } .shell { padding:0 14px 30px; } }
@media (prefers-reduced-motion:reduce) { *,*::before,*::after { scroll-behavior:auto !important; transition:none !important; animation:none !important; } }
@media print { body { background:white; } .topbar,#level-tabs,.navigator,.toolbar,.taxon-toolbar,details,.segment,.note { display:none !important; } .shell { max-width:none; padding:0; } #detail-panel { display:block; } .card { box-shadow:none; break-inside:avoid; } .plot-scroll { overflow:visible; } }
'''


def _write_html(level_data, overview, output_file, version, limits=None):
    limits = limits or {}
    multiple = len(level_data) > 1
    names = list(level_data)
    sample_count = len(set().union(*(set(item["samples"]) for item in level_data.values())))
    tabs = ''
    if multiple:
        tabs = '<nav id="level-tabs" role="tablist" aria-label="Stringency levels"><button type="button" class="level-tab" role="tab" id="tab-overview" aria-controls="overview-panel" aria-selected="true" data-level="overview"><strong>Overview</strong><small>Across all levels</small></button>'
        for index, name in enumerate(names):
            subtitle = _limit_text(limits.get(name))
            tabs += (f'<button type="button" class="level-tab" role="tab" id="tab-level-{index}" aria-controls="detail-panel" aria-selected="false" data-level="{index}"><strong>{escape(name)}</strong>'
                     + (f'<small>{escape(subtitle)}</small>' if subtitle else '') + '</button>')
        tabs += '</nav>'
    overview_panel = ''
    if multiple:
        kpis = f'<div class="card kpi"><div class="eyebrow">Samples</div><span class="kpi-value">{sample_count:,}</span><span class="kpi-note">Across this run</span></div>'
        for index, (name, item) in enumerate(overview["levels"].items()):
            totals = item["totals"]
            denom = totals["classified_fragments"] + totals["unclassified_fragments"]
            pct = _percent(totals["classified_fragments"] / denom if denom else None)
            kpis += (f'<div class="card kpi" data-level-index="{index}"><div class="eyebrow">{escape(name)}</div>'
                     f'<span class="kpi-value">{_number(totals["fragments_retained"])}</span>'
                     f'<span class="kpi-note">Fragments retained · {pct} classified</span></div>')
        overview_panel = f'''<section id="overview-panel" role="tabpanel" aria-labelledby="tab-overview">
          <div class="page-head"><div class="eyebrow">Run comparison</div><h1>Overview</h1><p class="subtle">How stringency changes retention, resolution and composition.</p></div>
          <div class="kpis">{kpis}</div>
          <div class="overview-grid">
            <article class="card span-all"><div class="card-head"><div><div class="eyebrow">01 / Summary</div><h2>Level comparison</h2></div><button type="button" id="download-comparison-csv" class="small-button">Download CSV</button></div>{_overview_table(overview, limits)}</article>
            <article class="card"><div class="card-head"><div><div class="eyebrow">02 / Taxonomy</div><h2>Resolution by rank</h2></div></div><div id="depth-chart" class="chart"></div></article>
            <article class="card"><div class="card-head"><div><div class="eyebrow">03 / Composition</div><h2>Kingdom composition</h2></div></div><div id="kingdom-chart" class="chart"></div></article>
            <article class="card span-all"><div class="card-head"><div><div class="eyebrow">04 / Samples</div><h2>Classified fragments by sample</h2></div><div class="toolbar"><label class="subtle" for="overview-search">Search</label><input id="overview-search" type="search" placeholder="Find a sample"><button type="button" id="download-sample-csv" class="small-button">Download CSV</button></div></div><div id="sample-table-count"></div><div class="table-wrap"><table id="sample-table"></table></div></article>
            <article class="card span-all"><div class="card-head"><div><div class="eyebrow">05 / Distribution</div><h2>Per-sample comparison</h2></div><div class="toolbar"><label class="subtle" for="dot-sort">Sort by</label><select id="dot-sort"></select><div class="segment" aria-label="Dot plot measure"><button type="button" id="dot-count" aria-pressed="true">Fragments</button><button type="button" id="dot-percent" aria-pressed="false">% classified</button></div></div></div><div id="classified-chart"></div></article>
          </div><p class="note">Ultrastrict requires at least 120 aligned bases while the default minimum read length is 100, so some reads can never pass that level. Stricter levels can lower the assigned rank; LCA over fewer alignments can also make assignments deeper. Counts are not simple subsets at every rank.</p>
        </section>'''
    now = datetime.now(UTC).astimezone().strftime("%Y-%m-%d %H:%M")
    first = names[0]
    first_sample = level_data[first]["samples"][0] if level_data[first]["samples"] else ""
    with open(output_file, "w", encoding="utf-8") as handle:
        logo = f'<img src="{_logo_uri()}" alt="BarcodeMapper">' if LOGO_PATH.exists() else 'BarcodeMapper'
        handle.write(f'''<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>BarcodeMapper · Report</title><script>{get_plotlyjs()}</script><style>{STYLE}</style></head><body>
<header class="topbar"><div class="brand"><a href="https://github.com/alberdilab/barcodemapper" target="_blank" rel="noopener noreferrer" aria-label="BarcodeMapper">{logo}</a><span>{escape(_version_label(version))}</span></div><div class="run-meta"><span>{sample_count:,} samples · {len(names)} {'level' if len(names)==1 else 'levels'}</span><time>{now}</time><button type="button" id="theme-toggle" class="icon-button" aria-label="Switch color theme">◐</button></div></header>
<div class="shell">{tabs}{overview_panel}
<section id="detail-panel" role="tabpanel"{' aria-labelledby="tab-level-0"' if multiple else ''}{' hidden' if multiple else ''}>
<div class="detail-main"><div class="page-head"><div class="eyebrow">Stringency level</div><h1 id="level-title">{escape(first)}</h1><p id="level-subtitle" class="subtle"></p></div>
<article class="card"><div class="card-head"><div><div class="eyebrow">01 / Composition</div><h2>Composition across samples</h2></div><div class="toolbar"><label class="subtle" for="composition-sort">Sort</label><select id="composition-sort"><option value="name">Name</option><option value="total">Total fragments</option><option value="kingdom">Largest kingdom share</option></select><div class="segment" aria-label="Composition measure"><button type="button" id="tab-absolute" aria-pressed="true">Absolute</button><button type="button" id="tab-relative" aria-pressed="false" title="Fraction of marker-matching fragments">Share</button></div></div></div><div id="composition-legend" class="composition-legend" role="group" aria-label="Toggle composition categories"></div><div class="plot-scroll"><div id="barplot"></div></div><details id="composition-table-view"><summary>Table view of composition data</summary><button type="button" id="download-composition-csv" class="small-button" title="Download the folded composition categories shown in this table">Download CSV</button><div class="table-wrap"><table id="composition-table"></table></div></details></article>
<article class="card"><div class="card-head"><div><div class="eyebrow">02 / Detection</div><h2>Taxonomy cladogram and detection grid</h2></div><span class="subtle">Cumulative fragments at each rank</span></div><div class="taxon-toolbar"><label class="subtle" for="taxon-search">Find taxon</label><input type="search" id="taxon-search" placeholder="Search taxonomy"><button type="button" id="taxon-expand" class="small-button">Expand all</button><button type="button" id="taxon-collapse" class="small-button">Collapse all</button><label class="subtle" for="taxon-depth">Depth</label><select id="taxon-depth"><option value="0">Kingdom</option><option value="1" selected>Phylum</option><option value="2">Class</option><option value="3">Order</option><option value="4">Family</option><option value="5">Genus</option><option value="6">Species</option></select><label class="subtle" for="taxon-sort">Samples</label><select id="taxon-sort"><option value="name">Name</option><option value="total">Total fragments</option><option value="taxon">Focused taxon</option></select><div class="segment" aria-label="Detection display"><button type="button" id="taxon-bubbles" aria-pressed="true">Bubbles</button><button type="button" id="taxon-heatmap" aria-pressed="false">Heatmap</button></div><div class="segment" aria-label="Detection value"><button type="button" id="taxon-fragments" aria-pressed="true">Fragments</button><button type="button" id="taxon-share" aria-pressed="false">Share of sample</button><button type="button" id="taxon-presence" aria-pressed="false">Presence</button></div><button type="button" id="download-taxon-visible-csv" class="small-button" title="Download currently visible taxonomy rows">Download visible CSV</button><button type="button" id="download-taxon-all-csv" class="small-button" title="Download every taxon at this level">Download all taxa CSV</button></div><div id="taxon-grid" class="taxon-grid" aria-label="Taxonomy detection grid"><div id="taxon-grid-inner" class="taxon-grid-inner"><div id="taxon-grid-head" class="taxon-grid-head"></div><div id="taxon-grid-body" role="tree" aria-label="Taxonomy tree"></div></div></div><div id="taxon-legend" class="taxon-legend"></div><div id="taxon-selection-note" class="taxon-selection-note" aria-live="polite" hidden><span id="taxon-selection-text"></span><button type="button" id="taxon-jump-sunburst">View sunburst</button></div></article>
<div class="sample-section"><aside class="navigator card" aria-label="Sample navigator"><div class="navigator-head"><h2>Samples</h2><div class="nav-buttons"><button type="button" id="sample-prev" class="small-button" aria-label="Previous sample">↑</button><button type="button" id="sample-next" class="small-button" aria-label="Next sample">↓</button></div></div><label class="subtle" for="sample-search">Search samples</label><input type="search" id="sample-search" placeholder="Filter samples"><span id="navigator-count" class="navigator-count"></span><div id="sample-buttons"></div></aside>
<article class="card" id="sunburst-card"><div class="card-head"><div><div class="eyebrow">03 / Taxonomy</div><h2 id="sample-title">{escape('Taxonomy of ' + first_sample) if first_sample else 'Taxonomy'}</h2><p class="subtle">Click a segment to zoom into its taxonomy.</p></div><div class="sample-metrics"><span><strong id="sample-total"></strong><br>Fragments</span><span><strong id="sample-classified"></strong><br>Classified</span></div></div><div id="sunburst"></div></article></div></div></section></div>''')
        for index, data in enumerate(level_data.values()):
            handle.write(f'<script type="application/json" id="data-level-{index}">{safe_json(data)}</script>\n')
        if multiple:
            ordered_overview = {**overview, "levels": [[name, overview["levels"][name]] for name in names]}
            handle.write(f'<script type="application/json" id="data-overview">{safe_json(ordered_overview)}</script>\n')
        handle.write(f'<script type="application/json" id="data-level-order">{safe_json(names)}</script>\n')
        handle.write(f'<script type="application/json" id="data-limits">{safe_json([_limit_text(limits.get(name)) for name in names])}</script>\n')
        handle.write(f'<script type="application/json" id="data-level-colors">{safe_json(LEVEL_PALETTES)}</script>\n')
        handle.write(REPORT_SCRIPT)
        handle.write('\n</body></html>\n')


def generate_html(samples, data_by_sample, output_file, abs_barplot, rel_barplot, version):
    """Keep the existing single-level Python entry point for callers."""
    taxonomy = list(dict.fromkeys(item for nodes in data_by_sample.values()
                                  for item in nodes.get("id", []) if item))
    frame = pd.DataFrame({"Taxonomy": taxonomy})
    for sample in samples:
        nodes = data_by_sample.get(sample, {})
        mapping = dict(zip(nodes.get("id", []), nodes.get("value", [])))
        frame[sample] = [mapping.get(path, 0) for path in taxonomy]
    _write_html({"custom": _compact_level(frame)}, None, output_file, version)


def generate_report(merged_by_level, extract_stats, output_file, version="unknown", limits=None):
    if not merged_by_level:
        raise ValueError("At least one input level is required")
    level_data = {name: _compact_level(frame) for name, frame in merged_by_level.items()}
    overview = overview_stats(merged_by_level, extract_stats) if len(level_data) > 1 else None
    _write_html(level_data, overview, output_file, version, limits)


def main():
    args = parse_args()
    limits = {}
    for item in args.limits:
        name, separator, thresholds = item.partition("=")
        parts = thresholds.split(":")
        if not separator or not name or len(parts) != 2 or name in limits:
            raise ValueError(f"Invalid or duplicate --limits: {item}")
        try:
            limits[name] = (float(parts[0]), int(parts[1]))
        except ValueError:
            raise ValueError(f"Invalid --limits thresholds: {item}") from None
    inputs = {}
    for item in args.input:
        if "=" in item:
            name, path = item.split("=", 1)
            if not name or not path:
                raise ValueError("--input requires NAME=PATH")
        else:
            if len(args.input) != 1:
                raise ValueError("A bare --input PATH is only valid for one level")
            name, path = "custom", item
        if name == "overview":
            raise ValueError("'overview' is reserved for the report tab")
        if name in inputs:
            raise ValueError(f"Duplicate level name: {name}")
        inputs[name] = load_merged(path)
    generate_report(inputs, load_stats(args.stats), args.output, args.version, limits)

REPORT_SCRIPT = r'''<script>
(() => {
    'use strict';
    const cache = {};
    const levelNames = JSON.parse(document.getElementById('data-level-order').textContent);
    const limits = JSON.parse(document.getElementById('data-limits').textContent);
    const levelPalettes = JSON.parse(document.getElementById('data-level-colors').textContent);
    const overviewElement = document.getElementById('data-overview');
    const overview = overviewElement ? JSON.parse(overviewElement.textContent) : null;
    const overviewLevels = overview ? overview.levels.map(([,data]) => data) : null;
    const ranks = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'];
    let currentLevel = null, currentSample = null, currentBarplot = 'absolute';
    let currentPanel = overview ? 'overview' : 'detail';
    let overviewDrawn = false, tableSort = { key: 'sample', descending: false };
    let dotSort = 0, dotMeasure = 'count';
    const hiddenComposition = {};
    const themePreference = window.matchMedia('(prefers-color-scheme: dark)');
    const taxonScroller=document.getElementById('taxon-grid');
    let treeWidth=0;
    try {
        const saved=Number(window.sessionStorage.getItem('barcodemapper-tree-width'));
        if (Number.isFinite(saved) && saved>=320) treeWidth=saved;
    } catch { /* File URLs may not grant session storage. */ }
    function clampTreeWidth(value) {
        return Math.max(320,Math.min(Math.max(320,Math.floor(taxonScroller.clientWidth*.45)),Math.round(value)));
    }
    function applyTreeWidth(value,persist=false) {
        treeWidth=clampTreeWidth(value || taxonScroller.clientWidth*.38);
        taxonScroller.style.setProperty('--tree-width',treeWidth+'px');
        const samples=currentLevel==null ? 0 : ensureTree(levelData(currentLevel)).order.length;
        document.getElementById('taxon-grid-inner').style.minWidth=(treeWidth+22*samples)+'px';
        const handle=document.getElementById('taxon-resize');
        if (handle) {
            handle.setAttribute('aria-valuenow',String(treeWidth));
            handle.setAttribute('aria-valuemax',String(Math.max(320,Math.floor(taxonScroller.clientWidth*.45))));
        }
        if (persist) try { window.sessionStorage.setItem('barcodemapper-tree-width',String(treeWidth)); } catch { /* Keep the width in memory. */ }
    }
    const count = value => Number(value || 0).toLocaleString();
    const pct = value => value == null ? 'n/a' : (value * 100).toFixed(1) + '%';
    const esc = value => String(value).replace(/[&<>"']/g, char => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[char]));
    function csvCell(value) {
        if (typeof value === 'number' && Number.isFinite(value)) return String(value);
        let cell = String(value ?? '');
        if (/^[=+\-@"\t\r\n]/.test(cell)) cell = "'" + cell;
        return '"' + cell.replace(/"/g,'""') + '"';
    }
    function csvRows(rows) {
        return '\ufeff' + rows.map(row => row.map(csvCell).join(',')).join('\r\n') + '\r\n';
    }
    function downloadCsv(level, table, rows) {
        const slug = value => String(value).replace(/[^a-zA-Z0-9._-]+/g,'_').replace(/^_+|_+$/g,'') || 'level';
        const url=URL.createObjectURL(new Blob([csvRows(rows)],{type:'text/csv;charset=utf-8'}));
        const link=document.createElement('a'); link.href=url;
        link.download='barcodemapper_'+slug(level)+'_'+slug(table)+'.csv';
        link.click(); setTimeout(()=>URL.revokeObjectURL(url),1000);
    }
    function effectiveTheme() {
        const explicit = document.documentElement.getAttribute('data-theme');
        return explicit === 'light' || explicit === 'dark' ? explicit : (themePreference.matches ? 'dark' : 'light');
    }
    const colors = () => effectiveTheme() === 'dark' ?
        {surface:'#1d2527',text:'#edf1ee',muted:'#a1b0aa',grid:'#35413f',accent:'#3987e5',neutral:'#889792'} :
        {surface:'#ffffff',text:'#202624',muted:'#66716d',grid:'#e8ece8',accent:'#2a78d6',neutral:'#9ba5a1'};
    function levelColor(index) {
        const palette = levelPalettes[effectiveTheme()];
        if (levelNames.length <= palette.length) return palette[index];
        const position = index / (levelNames.length - 1);
        const lightness = effectiveTheme() === 'dark' ? 82 - 52*position : 28 + 52*position;
        return `hsl(211 65% ${lightness.toFixed(1)}%)`;
    }
    function applyLevelAccents() {
        for (let index=0; index<levelNames.length; index++) {
            document.querySelectorAll(`[data-level-index="${index}"]`).forEach(element =>
                element.style.setProperty('--level-accent',levelColor(index)));
            document.getElementById('tab-level-'+index)?.style.setProperty('--level-accent',levelColor(index));
        }
    }
    function kingdomColor(kingdom) {
        const palette = effectiveTheme() === 'dark' ? {'k__Animalia':'#3987e5','k__Fungi':'#d95926','k__Viridiplantae':'#199e70'}
                                                    : {'k__Animalia':'#2a78d6','k__Fungi':'#eb6834','k__Viridiplantae':'#1baf7a'};
        return palette[kingdom] || colors().neutral;
    }
    function taxonShade(kingdom, index, total) {
        const hues = {'k__Animalia':211,'k__Fungi':17,'k__Viridiplantae':160};
        if (!(kingdom in hues)) return colors().neutral;
        const start = effectiveTheme() === 'dark' ? 40 : 30;
        const end = effectiveTheme() === 'dark' ? 85 : 75;
        const lightness = total <= 1 ? (start + end)/2 : start + (end-start)*index/(total-1);
        return `hsl(${hues[kingdom]} 67% ${lightness.toFixed(1)}%)`;
    }
    function compositionColors(data) {
        const sizes = {}, seen = {}, map = new Map();
        data.composition.forEach(item => { sizes[item.kingdom] = (sizes[item.kingdom] || 0) + 1; });
        data.composition.forEach(item => {
            const index = seen[item.kingdom] || 0; seen[item.kingdom] = index+1;
            map.set(item.kingdom + '\u0000' + item.name, item.name === 'Unclassified' ? colors().neutral : taxonShade(item.kingdom,index,sizes[item.kingdom]));
        });
        return map;
    }
    function taxonGridColor(data, palette, index) {
        const taxa = data.taxa;
        let node = index;
        if (taxa.labels[index] === 'Unclassified') return colors().neutral;
        while (taxa.parents[node] >= 0 && taxa.depths[node] > 1) node = taxa.parents[node];
        const root = taxa.parents[node] >= 0 ? taxa.parents[node] : node;
        if (taxa.parents[node] < 0) return kingdomColor(taxa.labels[node]);
        return palette.get(taxa.labels[root] + '\u0000' + taxa.labels[node]) || kingdomColor(taxa.labels[root]);
    }
    function layout(extra={}) {
        const c = colors();
        return Object.assign({paper_bgcolor:'rgba(0,0,0,0)', plot_bgcolor:'rgba(0,0,0,0)', font:{family:'-apple-system, Segoe UI, Arial, sans-serif',color:c.text,size:11},
            margin:{l:54,r:22,t:18,b:50}, hoverlabel:{bgcolor:c.surface,font:{color:c.text},bordercolor:c.grid},
            xaxis:{gridcolor:c.grid,zerolinecolor:c.grid}, yaxis:{gridcolor:c.grid,zerolinecolor:c.grid},
            showlegend:true, legend:{orientation:'h',y:1.15,x:0,font:{size:10}}, autosize:true}, extra);
    }
    const plotConfig = {displaylogo:false,responsive:true,modeBarButtonsToRemove:['lasso2d','select2d']};
    function drawPlot(element, traces, figureLayout) {
        Plotly.purge(element);
        element.replaceChildren();
        return Plotly.newPlot(element,traces,figureLayout,plotConfig);
    }
    function showEmpty(element, message) {
        Plotly.purge(element);
        element.replaceChildren();
        const paragraph = document.createElement('p');
        paragraph.className='empty'; paragraph.textContent=message;
        element.append(paragraph);
    }
    function levelData(index) {
        if (!cache[index]) cache[index] = JSON.parse(document.getElementById('data-level-' + index).textContent);
        return cache[index];
    }
    function setTabs(panel, index) {
        document.querySelectorAll('.level-tab').forEach(tab => {
            const selected = panel === 'overview' ? tab.dataset.level === 'overview' : tab.dataset.level === String(index);
            tab.setAttribute('aria-selected', String(selected));
            tab.tabIndex = selected ? 0 : -1;
        });
        const overviewPanel = document.getElementById('overview-panel');
        if (overviewPanel) overviewPanel.hidden = panel !== 'overview';
        const detail = document.getElementById('detail-panel');
        detail.hidden = panel !== 'detail';
        if (panel === 'detail' && overview) detail.setAttribute('aria-labelledby', 'tab-level-' + index);
        currentPanel = panel;
    }
    function putText(id, value) { document.getElementById(id).textContent = value; }
    function sampleRows() {
        if (!overview) return [];
        const sampleNames = [...new Set(overview.classified_by_sample.map(row => row.sample))];
        return sampleNames.map(sample => ({name:sample, values:levelNames.map((_,index) => overviewLevels[index].samples[sample]?.classified_fragments || 0),
            fractions:levelNames.map((_,index) => { const item = overviewLevels[index].samples[sample];
                const denominator = (item?.classified_fragments || 0) + (item?.unclassified_fragments || 0);
                return denominator ? item.classified_fragments / denominator : 0; })}));
    }
    function filteredRows() {
        const query = document.getElementById('overview-search').value.toLocaleLowerCase();
        return sampleRows().filter(row => row.name.toLocaleLowerCase().includes(query));
    }
    function renderSampleTable() {
        const rows = filteredRows();
        const key = tableSort.key;
        rows.sort((a,b) => {
            const result = key === 'sample' ? a.name.localeCompare(b.name) :
                a.values[Number(key.slice('level:'.length))] - b.values[Number(key.slice('level:'.length))];
            return tableSort.descending ? -result : result;
        });
        const table = document.getElementById('sample-table');
        table.replaceChildren();
        const header = table.createTHead().insertRow();
        const columns = [{key:'sample',label:'Sample'},...levelNames.map((label,index) => ({key:'level:'+index,label}))];
        for (const {key,label} of columns) {
            const th = document.createElement('th'); th.scope = 'col';
            if (tableSort.key === key) th.setAttribute('aria-sort',tableSort.descending ? 'descending' : 'ascending');
            const button = document.createElement('button'); button.type = 'button'; button.className = 'table-sort';
            button.textContent = label;
            button.setAttribute('aria-label','Sort by ' + (key === 'sample' ? 'sample name' : label));
            button.addEventListener('click', () => { tableSort.descending = tableSort.key === key ? !tableSort.descending : key !== 'sample'; tableSort.key = key; renderSampleTable(); });
            th.append(button); header.append(th);
        }
        const body = table.createTBody();
        const max = levelNames.map((_,index) => Math.max(1,...rows.map(row => row.values[index])));
        for (const row of rows) {
            const tr = body.insertRow();
            const name = tr.insertCell(); name.textContent = row.name; name.title = row.name;
            row.values.forEach((value,index) => {
                const cell = tr.insertCell(); cell.className = 'bar-cell';
                const bar = document.createElement('span'); bar.className = 'cell-bar'; bar.style.width = (100*value/max[index]) + '%';
                const label = document.createElement('span'); label.textContent = count(value);
                cell.append(bar,label);
            });
        }
        const allCount = sampleRows().length;
        putText('sample-table-count', rows.length + ' of ' + allCount + ' samples');
        const download = document.getElementById('download-sample-csv');
        const filtered = document.getElementById('overview-search').value.length > 0;
        download.textContent = filtered ? 'Download filtered CSV' : 'Download all CSV';
        download.title = `Export ${rows.length} ${filtered ? 'filtered' : 'total'} samples`;
        return rows;
    }
    function drawDots() {
        const rows = filteredRows();
        const sortIndex = dotSort;
        rows.sort((a,b) => (dotMeasure === 'count' ? b.values[sortIndex]-a.values[sortIndex] : b.fractions[sortIndex]-a.fractions[sortIndex]) || a.name.localeCompare(b.name));
        const height = Math.max(340, rows.length*22 + 95);
        const el = document.getElementById('classified-chart'); el.style.height = height + 'px';
        if (!rows.length) { showEmpty(el,'No matching samples.'); return; }
        const traces = levelNames.map((level,index) => ({type:'scatter',mode:'markers',name:level,
            x:rows.map(row => dotMeasure === 'count' ? row.values[index] : row.fractions[index]),
            y:rows.map(row => row.name), marker:{size:9,color:levelColor(index),symbol:['circle','square','diamond','triangle-up'][index % 4]},
            customdata:rows.map(row => esc(row.name)), hovertemplate:'%{customdata}<br>'+esc(level)+': %{x'+(dotMeasure === 'count' ? ':,.0f' : ':.1%')+'}<extra></extra>'}));
        drawPlot(el,traces,layout({height,margin:{l:Math.min(230,Math.max(115,...rows.map(r => r.name.length*7))),r:25,t:20,b:50},
            xaxis:{title:dotMeasure === 'count' ? 'Classified fragments' : 'Classified share',tickformat:dotMeasure === 'count' ? ',d' : '.0%',gridcolor:colors().grid},
            yaxis:{type:'category',categoryorder:'array',categoryarray:rows.map(row => row.name).reverse(),gridcolor:colors().grid,automargin:false},
            legend:{orientation:'h',y:1.04,x:0,font:{size:10}}}));
    }
    function drawOverview() {
        if (!overview) return;
        const c = colors();
        const depth = levelNames.map((level,index) => {
            const totals = overviewLevels[index].totals;
            const denominator = totals.classified_fragments + totals.unclassified_fragments;
            return {type:'scatter',mode:'lines+markers',name:level,x:ranks,
                y:totals.depth_distribution.map(value => denominator ? value/denominator : 0),
                line:{color:levelColor(index),width:2.5},marker:{size:8,symbol:['circle','square','diamond','triangle-up'][index % 4]},
                hovertemplate:'%{x}: %{y:.1%}<extra>'+esc(level)+'</extra>'};
        });
        drawPlot(document.getElementById('depth-chart'),depth,layout({margin:{l:50,r:20,t:48,b:45},
            xaxis:{type:'category',categoryorder:'array',categoryarray:ranks,gridcolor:c.grid},
            yaxis:{tickformat:'.0%',range:[0,1],gridcolor:c.grid},
            legend:{orientation:'h',x:0,y:1.12,font:{size:10}}}));
        const kingdoms = [...new Set(overviewLevels.flatMap(item => Object.keys(item.totals.kingdom_composition)))];
        const kingdomOrder = ['k__Animalia','k__Fungi','k__Viridiplantae'];
        kingdoms.sort((a,b) => (kingdomOrder.indexOf(a) < 0 ? 3 : kingdomOrder.indexOf(a)) -
            (kingdomOrder.indexOf(b) < 0 ? 3 : kingdomOrder.indexOf(b)) || a.localeCompare(b));
        const composition = kingdoms.map(kingdom => ({type:'bar',name:kingdom.replace(/^k__/,''),x:levelNames,
            y:levelNames.map((_,index) => { const totals = overviewLevels[index].totals;
                return totals.classified_fragments ? (totals.kingdom_composition[kingdom] || 0)/totals.classified_fragments : 0; }),
            marker:{color:kingdomColor(kingdom),line:{color:c.surface,width:2}},hovertemplate:'%{x}<br>%{y:.1%}<extra>'+esc(kingdom)+'</extra>'}));
        drawPlot(document.getElementById('kingdom-chart'),composition,layout({barmode:'stack',margin:{l:50,r:20,t:30,b:45},
            xaxis:{type:'category',categoryorder:'array',categoryarray:levelNames,gridcolor:c.grid},
            yaxis:{tickformat:'.0%',range:[0,1],gridcolor:c.grid}}));
        renderSampleTable(); drawDots(); overviewDrawn = true;
    }
    function resizeVisibleCharts() {
        const ids = currentPanel === 'overview' ? ['depth-chart','kingdom-chart','classified-chart'] : ['barplot','sunburst'];
        for (const id of ids) {
            const chart = document.getElementById(id);
            if (chart && !chart.closest('[hidden]') && chart._fullLayout) Plotly.Plots.resize(chart);
        }
    }
    function showOverview() {
        setTabs('overview');
        if (!overviewDrawn) drawOverview();
        else requestAnimationFrame(resizeVisibleCharts);
    }
    function compositionOrder(data) {
        const indices = data.samples.map((_,index) => index);
        const sort = document.getElementById('composition-sort').value;
        const shares = data.samples.map((_,index) => {
            const byKingdom = {};
            for (const item of data.composition) {
                if (item.kingdom !== 'Unclassified') byKingdom[item.kingdom] = (byKingdom[item.kingdom] || 0) + item.values[index];
            }
            return data.totals[index] ? Math.max(0,...Object.values(byKingdom))/data.totals[index] : 0;
        });
        indices.sort((a,b) => sort === 'total' ? data.totals[b]-data.totals[a] || data.samples[a].localeCompare(data.samples[b]) :
            sort === 'kingdom' ? shares[b]-shares[a] || data.samples[a].localeCompare(data.samples[b]) : data.samples[a].localeCompare(data.samples[b]));
        return indices;
    }
    function drawBarplot() {
        if (currentLevel == null) return;
        const data = levelData(currentLevel), order = compositionOrder(data), c = colors();
        const samples = order.map(index => data.samples[index]);
        const width = Math.max(document.querySelector('.plot-scroll').clientWidth - 4, samples.length*19 + 135);
        const el = document.getElementById('barplot'); el.style.width = width+'px';
        const legend = document.getElementById('composition-legend');
        legend.replaceChildren();
        if (!data.composition.length || !samples.length) {
            document.getElementById('composition-table').replaceChildren();
            showEmpty(el,'No composition data for this level.');
            return;
        }
        const hidden = hiddenComposition[currentLevel] ||= new Set();
        const palette = compositionColors(data);
        const traceColors = [];
        const traces = data.composition.map(item => {
            const color = palette.get(item.kingdom + '\u0000' + item.name);
            traceColors.push(color);
            const key = item.kingdom + '\u0000' + item.name;
            return {type:'bar',name:item.name.replace(/^[pk]__/,''),x:samples.map((_,position) => position),
                y:order.map(i => currentBarplot === 'absolute' ? item.values[i] : (data.totals[i] ? item.values[i]/data.totals[i] : 0)),
                customdata:order.map(i => [esc(data.samples[i]),esc(item.name),esc(item.kingdom.replace(/^k__/,'')),data.totals[i] ? item.values[i]/data.totals[i] : 0,item.values[i]]),
                visible:!hidden.has(key), marker:{color,line:{color:c.surface,width:2}},
                hovertemplate:'<b>%{customdata[0]}</b><br>%{customdata[1]} · %{customdata[2]}<br>%{customdata[4]:,} fragments · %{customdata[3]:.1%} of sample<extra></extra>'};
        });
        renderCompositionLegend(data,traceColors,hidden);
        drawPlot(el,traces,layout({width,height:560,barmode:'stack',showlegend:false,margin:{l:58,r:12,t:12,b:108},
            xaxis:{tickmode:'array',tickvals:samples.map((_,i) => i),ticktext:samples.map(name => name.length > 17 ? name.slice(0,15)+'…' : name),tickangle:-45,gridcolor:c.grid},
            yaxis:{title:currentBarplot === 'absolute' ? 'Fragments' : 'Share of sample',tickformat:currentBarplot === 'absolute' ? ',d' : '.0%',gridcolor:c.grid}}));
        drawCompositionTable(data,order);
    }
    function renderCompositionLegend(data,traceColors,hidden) {
        const legend = document.getElementById('composition-legend');
        const groups = new Map();
        data.composition.forEach((item,index) => {
            if (!groups.has(item.kingdom)) groups.set(item.kingdom,[]);
            groups.get(item.kingdom).push({item,index});
        });
        for (const [kingdom,entries] of groups) {
            const group = document.createElement('div'); group.className='legend-group';
            group.setAttribute('role','group'); group.setAttribute('aria-label',kingdom.replace(/^k__/,''));
            const heading = document.createElement('div'); heading.className='legend-heading';
            const headingSwatch = document.createElement('span'); headingSwatch.className='legend-swatch';
            headingSwatch.style.setProperty('--swatch',kingdomColor(kingdom));
            const headingLabel = document.createElement('span'); headingLabel.textContent=kingdom.replace(/^k__/, '');
            heading.append(headingSwatch,headingLabel); group.append(heading);
            const items = document.createElement('div'); items.className='legend-items';
            for (const {item,index} of entries) {
                const key = item.kingdom + '\u0000' + item.name;
                const chip = document.createElement('button'); chip.type='button'; chip.className='legend-chip';
                chip.setAttribute('aria-pressed',String(!hidden.has(key)));
                chip.title=item.name + ' · ' + kingdom.replace(/^k__/,'');
                const swatch = document.createElement('span'); swatch.className='legend-swatch';
                swatch.style.setProperty('--swatch',traceColors[index]);
                const label = document.createElement('span'); label.className='legend-chip-label';
                label.textContent=item.name.replace(/^p__/, '').replace(/k__/g,'');
                chip.append(swatch,label);
                chip.addEventListener('click',()=>{
                    const visible = hidden.has(key);
                    if (visible) hidden.delete(key); else hidden.add(key);
                    chip.setAttribute('aria-pressed',String(visible));
                    Plotly.restyle('barplot',{visible},[index]);
                });
                items.append(chip);
            }
            group.append(items); legend.append(group);
        }
    }
    function drawCompositionTable(data,order) {
        const table = document.getElementById('composition-table'); table.replaceChildren();
        const header = table.createTHead().insertRow();
        for (const text of ['Sample',...data.composition.map(item => item.name)]) {
            const th = document.createElement('th'); th.scope='col'; th.textContent=text; header.append(th);
        }
        const body = table.createTBody();
        for (const index of order) {
            const row = body.insertRow(); row.insertCell().textContent = data.samples[index];
            data.composition.forEach(item => { row.insertCell().textContent = count(item.values[index]); });
        }
    }
    function taxonCount(pairs, sample) {
        let low=0, high=pairs.length-1;
        while (low<=high) {
            const middle=(low+high)>>1;
            if (pairs[middle][0]===sample) return pairs[middle][1];
            if (pairs[middle][0]<sample) low=middle+1; else high=middle-1;
        }
        return 0;
    }
    function ensureTree(data) {
        if (data._tree) return data._tree;
        const taxa=data.taxa, children=taxa.labels.map(()=>[]), roots=[], lineages=[];
        taxa.parents.forEach((parent,index)=>{
            if (parent<0) roots.push(index); else children[parent].push(index);
            lineages[index]=(parent<0 ? '' : lineages[parent]+'; ') + taxa.labels[index].replace(/^[kpcofgs]__/, '');
        });
        data._tree={children,roots,lineages,expanded:new Set(roots),effectiveExpanded:new Set(roots),
            visible:[],order:[],focused:roots[0] ?? -1,window:null};
        return data._tree;
    }
    function taxonSampleOrder(data,tree,sort) {
        const order=data.samples.map((_,index)=>index);
        order.sort((a,b)=>{
            const difference=sort==='total' ? data.totals[b]-data.totals[a] :
                sort==='taxon' && tree.focused>=0 ?
                    taxonCount(data.taxa.values[tree.focused],b)-taxonCount(data.taxa.values[tree.focused],a) : 0;
            return difference || data.samples[a].localeCompare(data.samples[b],undefined,{numeric:false});
        });
        return order;
    }
    function visibleTaxa(data,tree,search) {
        const query=search.toLocaleLowerCase().trim();
        const included=new Set();
        if (query) data.taxa.labels.forEach((label,index)=>{
            if (!tree.lineages[index].toLocaleLowerCase().includes(query)) return;
            for (let ancestor=index; ancestor>=0; ancestor=data.taxa.parents[ancestor]) included.add(ancestor);
        });
        const visible=[], expanded=query ? new Set() : tree.expanded;
        function walk(index) {
            if (query && !included.has(index)) return;
            visible.push(index);
            if (query) {
                if (tree.children[index].some(child=>included.has(child))) expanded.add(index);
                tree.children[index].forEach(walk);
            } else if (expanded.has(index)) tree.children[index].forEach(walk);
        }
        tree.roots.forEach(walk);
        return {visible,expanded,searching:Boolean(query)};
    }
    function expandedToDepth(taxa,children,depth) {
        return new Set(taxa.depths.flatMap((value,index)=>value<depth && children[index].length ? [index] : []));
    }
    function changeTaxonExpansion(data,tree,index,expand) {
        for (let ancestor=data.taxa.parents[index];ancestor>=0;ancestor=data.taxa.parents[ancestor]) tree.expanded.add(ancestor);
        if (expand) tree.expanded.add(index); else tree.expanded.delete(index);
    }
    function taxonWindow(length,samples,scrollTop,height) {
        const budget=Math.max(1,Math.floor(6000/Math.max(1,samples)));
        if (length<=budget) return {first:0,last:length,virtual:false};
        const viewport=Math.max(0,Math.floor(scrollTop/30));
        const rows=Math.ceil(height/30)+2;
        const span=Math.min(budget,rows+12);
        const first=Math.max(0,Math.min(length-span,viewport-6));
        return {first,last:Math.min(length,first+span),virtual:true};
    }
    function windowCovers(window,scrollTop,height,length) {
        if (!window || !window.virtual) return Boolean(window);
        const first=Math.max(0,Math.floor(scrollTop/30));
        const last=Math.min(length,first+Math.ceil(height/30)+2);
        return window.first<=first && window.last>=last;
    }
    function rovingTaxon(visible,first,last,focused) {
        if (first>=last) return -1;
        const position=visible.indexOf(focused);
        return visible[Math.min(last-1,Math.max(first,position<0 ? first : position))];
    }
    function taxonShare(amount,total) {
        return total>0 ? amount/total : 0;
    }
    function taxonCsvData(data,tree,all,ranks) {
        const indices=all ? data.taxa.labels.map((_,index)=>index) : tree.visible;
        return [['Taxon lineage','Rank',...tree.order.map(index=>data.samples[index])],
            ...indices.map(index=>[tree.lineages[index],ranks[data.taxa.depths[index]] || 'Taxon',
                ...tree.order.map(sample=>taxonCount(data.taxa.values[index],sample))])];
    }
    function clearTaxonSearch() {
        const input=document.getElementById('taxon-search');
        if (input.value) input.value='';
    }
    function setTaxonDepth(depth) {
        const data=levelData(currentLevel), tree=ensureTree(data);
        clearTaxonSearch();
        tree.expanded=expandedToDepth(data.taxa,tree.children,depth);
        document.getElementById('taxon-depth').value=String(Math.min(6,Math.max(0,depth)));
        document.getElementById('taxon-grid').scrollTop=0;
        drawTaxonGrid();
    }
    function focusTaxon(index) {
        const data=levelData(currentLevel), tree=ensureTree(data), position=tree.visible.indexOf(index);
        if (position<0) return;
        tree.focused=index;
        const scroller=document.getElementById('taxon-grid');
        const top=96+position*30;
        if (top<scroller.scrollTop+96) scroller.scrollTop=Math.max(0,top-96);
        else if (top+30>scroller.scrollTop+scroller.clientHeight) scroller.scrollTop=top+30-scroller.clientHeight;
        if (document.getElementById('taxon-sort').value==='taxon') drawTaxonGrid();
        else if (!document.querySelector(`#taxon-grid-body [data-taxon="${index}"]`)) renderTaxonRows();
        document.querySelectorAll('#taxon-grid-body .taxon-node').forEach(button=>
            { button.tabIndex=Number(button.dataset.taxon)===tree.focused ? 0 : -1; });
        document.querySelector(`#taxon-grid-body [data-taxon="${index}"]`)?.focus({preventScroll:true});
    }
    function renderTaxonRows() {
        if (currentLevel==null) return;
        const data=levelData(currentLevel), tree=ensureTree(data), body=document.getElementById('taxon-grid-body');
        const visible=tree.visible, order=tree.order, scroller=document.getElementById('taxon-grid');
        const scrollTop=scroller.scrollTop, clientHeight=scroller.clientHeight;
        const hadTreeFocus=body.contains(document.activeElement);
        body.replaceChildren();
        if (!visible.length || !order.length || !data.totals.some(value=>value>0)) {
            const empty=document.createElement('p'); empty.className='empty';
            empty.textContent=!data.samples.length ? 'No samples in this level.' :
                !data.taxa.labels.length || !data.totals.some(value=>value>0) ? 'No taxonomy data for this level.' : 'No matching taxa.';
            body.append(empty); return;
        }
        const window=taxonWindow(visible.length,order.length,scrollTop,clientHeight);
        const {first,last}=window;
        tree.window=window;
        tree.focused=rovingTaxon(visible,first,last,tree.focused);
        const spacer=height=>{ const div=document.createElement('div'); div.className='taxon-spacer'; div.style.height=height+'px'; return div; };
        if (first) body.append(spacer(first*30));
        const fragment=document.createDocumentFragment(), mode=document.getElementById('taxon-heatmap').getAttribute('aria-pressed')==='true' ? 'heatmap':'bubbles';
        const valueMode=['fragments','share','presence'].find(key=>document.getElementById('taxon-'+key).getAttribute('aria-pressed')==='true');
        const maxTotal=Math.max(1,...data.totals), dark=effectiveTheme()==='dark', palette=compositionColors(data);
        for (let position=first;position<last;position++) {
            const index=visible[position], label=data.taxa.labels[index], depth=data.taxa.depths[index];
            const row=document.createElement('div'); row.className='taxon-row'; row.style.gridTemplateColumns=`var(--tree-width) repeat(${order.length},22px)`;
            const taxonColor=taxonGridColor(data,palette,index);
            const left=document.createElement('div'); left.className='taxon-tree-cell';
            const node=document.createElement('button'); node.type='button'; node.className='taxon-node'; node.role='treeitem';
            node.dataset.taxon=String(index); node.style.setProperty('--depth',String(depth));
            node.setAttribute('aria-level',String(depth+1)); node.title=tree.lineages[index]+' · '+(ranks[depth] || 'Taxon');
            node.tabIndex=index===tree.focused ? 0 : -1;
            const parent=data.taxa.parents[index], siblings=parent<0 ? tree.roots : tree.children[parent];
            node.setAttribute('aria-setsize',String(siblings.length));
            node.setAttribute('aria-posinset',String(siblings.indexOf(index)+1));
            if (tree.children[index].length) node.setAttribute('aria-expanded',String(tree.effectiveExpanded.has(index)));
            const toggle=document.createElement('span'); toggle.className='taxon-toggle'; toggle.setAttribute('aria-hidden','true');
            toggle.textContent=tree.children[index].length ? tree.effectiveExpanded.has(index) ? '−' : '+' : '·';
            const rank=document.createElement('span'); rank.className='taxon-rank'; rank.textContent=(ranks[depth] || 'T')[0]; rank.style.color=taxonColor; rank.style.fontWeight='700';
            const name=document.createElement('span'); name.className='taxon-node-label';
            const plain=label.replace(/^[kpcofgs]__/, '');
            name.textContent=depth>=4 && plain.length>34 ? plain.slice(0,13)+'…'+plain.slice(-18) : plain;
            node.append(toggle,rank,name); left.append(node); row.append(left);
            node.addEventListener('click',()=>{
                tree.focused=index;
                if (tree.children[index].length) {
                    const wasExpanded=tree.effectiveExpanded.has(index);
                    clearTaxonSearch();
                    changeTaxonExpansion(data,tree,index,!wasExpanded);
                    drawTaxonGrid();
                }
                focusTaxon(index);
            });
            for (const sample of order) {
                const amount=taxonCount(data.taxa.values[index],sample), share=taxonShare(amount,data.totals[sample]);
                const cell=document.createElement('div'); cell.className='taxon-cell'; cell.role='img';
                cell.dataset.sampleIndex=String(sample); cell.classList.toggle('selected',data.samples[sample]===currentSample);
                const description=`${tree.lineages[index]}; ${data.samples[sample]}; ${amount} fragments; ${(share*100).toFixed(1)}% of sample`;
                cell.title=description; cell.setAttribute('aria-label',description);
                if (amount>0) {
                    const strength=valueMode==='presence' ? 1 : valueMode==='share' ? Math.sqrt(share) : Math.log1p(amount)/Math.log1p(maxTotal);
                    const mark=document.createElement('span'); mark.className=mode==='heatmap' ? 'taxon-heat':'taxon-bubble';
                    mark.style.background=taxonColor;
                    if (mode==='heatmap') mark.style.opacity=String(.15+.85*Math.max(.12,strength));
                    else { const radius=2.5+6.5*Math.sqrt(strength); mark.style.width=2*radius+'px'; mark.style.height=2*radius+'px'; }
                    cell.append(mark);
                }
                row.append(cell);
            }
            fragment.append(row);
        }
        body.append(fragment);
        if (last<visible.length) body.append(spacer((visible.length-last)*30));
        if (hadTreeFocus) body.querySelector(`[data-taxon="${tree.focused}"]`)?.focus({preventScroll:true});
    }
    function drawTaxonGrid() {
        if (currentLevel==null) return;
        const data=levelData(currentLevel), tree=ensureTree(data), head=document.getElementById('taxon-grid-head');
        const search=document.getElementById('taxon-search').value;
        const visibility=visibleTaxa(data,tree,search);
        tree.visible=visibility.visible; tree.effectiveExpanded=visibility.expanded; tree.searching=visibility.searching;
        tree.order=taxonSampleOrder(data,tree,document.getElementById('taxon-sort').value);
        head.replaceChildren(); head.style.gridTemplateColumns=`var(--tree-width) repeat(${tree.order.length},22px)`;
        const corner=document.createElement('div'); corner.className='taxon-tree-cell'; corner.textContent='Taxonomy';
        const handle=document.createElement('span'); handle.id='taxon-resize'; handle.className='taxon-resize';
        handle.role='separator'; handle.tabIndex=0; handle.setAttribute('aria-orientation','vertical');
        handle.setAttribute('aria-label','Resize taxonomy column'); handle.setAttribute('aria-valuemin','320');
        let startX=0, startWidth=0;
        handle.addEventListener('pointerdown',event=>{startX=event.clientX;startWidth=treeWidth;handle.setPointerCapture(event.pointerId);event.preventDefault();});
        handle.addEventListener('pointermove',event=>{if (handle.hasPointerCapture(event.pointerId)) applyTreeWidth(startWidth+event.clientX-startX);});
        handle.addEventListener('pointerup',event=>{if (handle.hasPointerCapture(event.pointerId)) {handle.releasePointerCapture(event.pointerId);applyTreeWidth(treeWidth,true);}});
        handle.addEventListener('keydown',event=>{
            if (!['ArrowLeft','ArrowRight'].includes(event.key)) return;
            event.preventDefault(); applyTreeWidth(treeWidth+(event.key==='ArrowRight' ? 20:-20),true);
        });
        corner.append(handle); head.append(corner);
        for (const index of tree.order) {
            const button=document.createElement('button'); button.type='button'; button.className='taxon-sample-head';
            button.dataset.sampleIndex=String(index);
            button.title=data.samples[index]; button.setAttribute('aria-label','Show sunburst for '+data.samples[index]);
            const name=document.createElement('span'); name.textContent=data.samples[index]; button.append(name);
            button.addEventListener('click',()=>loadSample(data.samples[index])); head.append(button);
        }
        applyTreeWidth(treeWidth);
        renderTaxonRows();
        const mode=['fragments','share','presence'].find(key=>document.getElementById('taxon-'+key).getAttribute('aria-pressed')==='true');
        const display=document.getElementById('taxon-heatmap').getAttribute('aria-pressed')==='true' ? 'Heatmap colour' : 'Bubble area';
        document.getElementById('taxon-legend').innerHTML='';
        const swatch=document.createElement('span'); swatch.className='taxon-legend-swatch';
        document.getElementById('taxon-legend').append(swatch,document.createTextNode(display+' · '+({fragments:'log fragments',share:'share of sample',presence:'detected or absent'}[mode])+' · '+tree.visible.length+' of '+data.taxa.labels.length+(tree.searching ? ' taxa matching search' : ' taxa visible')));
        syncTaxonSelection();
    }
    function taxonCsvRows(all) {
        return taxonCsvData(levelData(currentLevel),ensureTree(levelData(currentLevel)),all,ranks);
    }
    function syncTaxonSelection() {
        if (currentLevel==null) return;
        const data=levelData(currentLevel), selected=data.samples.indexOf(currentSample);
        document.querySelectorAll('#taxon-grid [data-sample-index]').forEach(element=>{
            const active=Number(element.dataset.sampleIndex)===selected && selected>=0;
            element.classList.toggle('selected',active);
            if (element.classList.contains('taxon-sample-head')) element.setAttribute('aria-pressed',String(active));
        });
        const note=document.getElementById('taxon-selection-note'); note.hidden=selected<0;
        if (selected>=0) putText('taxon-selection-text','Showing '+currentSample+' in the sunburst below');
    }
    function renderNavigator() {
        const data = levelData(currentLevel), list = document.getElementById('sample-buttons');
        const query = document.getElementById('sample-search').value.toLocaleLowerCase();
        list.replaceChildren();
        const maximum = Math.max(1,...data.totals);
        data.samples.forEach((sample,index) => {
            if (!sample.toLocaleLowerCase().includes(query)) return;
            const button = document.createElement('button'); button.type='button'; button.className='sample-button';
            button.dataset.sample = sample; button.title=sample; button.setAttribute('aria-current',String(sample === currentSample));
            button.classList.toggle('active',sample === currentSample);
            const label = document.createElement('span'); label.className='sample-name'; label.textContent=sample;
            const meter = document.createElement('span'); meter.className='sample-meter';
            const fill = document.createElement('span'); fill.style.width=(100*data.totals[index]/maximum)+'%'; meter.append(fill);
            button.append(label,meter); button.addEventListener('click',()=>loadSample(sample)); list.append(button);
        });
        const at = data.samples.indexOf(currentSample);
        putText('navigator-count', data.samples.length ? (at+1)+' of '+data.samples.length+' samples · '+list.children.length+' shown' : '0 samples');
        document.getElementById('sample-prev').disabled = at <= 0;
        document.getElementById('sample-next').disabled = at < 0 || at >= data.samples.length-1;
    }
    function loadSample(sample) {
        const data = levelData(currentLevel), index = data.samples.indexOf(sample);
        currentSample = index >= 0 ? sample : data.samples[0] || null;
        renderNavigator();
        syncTaxonSelection();
        putText('sample-title',currentSample ? 'Taxonomy of '+currentSample : 'Taxonomy');
        const active = data.samples.indexOf(currentSample);
        putText('sample-total',active >= 0 ? count(data.totals[active]) : '0');
        putText('sample-classified',active >= 0 && data.totals[active] ? pct(data.classified[active]/data.totals[active]) : 'n/a');
        const el = document.getElementById('sunburst');
        const taxa=data.taxa;
        if (active < 0 || !data.totals[active]) { showEmpty(el,'No data for this sample.'); return; }
        const present=index=>{ for (let node=index; node>=0; node=taxa.parents[node]) if (!(taxonCount(taxa.values[node],active)>0)) return false; return true; };
        const indices=taxa.labels.map((_,index)=>index).filter(present);
        const palette=compositionColors(data);
        const trace = {type:'sunburst',ids:indices.map(String),labels:indices.map(i=>taxa.labels[i].replace(/^[kpcofgs]__/,'')),
            parents:indices.map(i=>taxa.parents[i] < 0 ? '' : String(taxa.parents[i])),
            values:indices.map(i=>taxonCount(taxa.values[i],active)),branchvalues:'total',maxdepth:5,
            marker:{colors:indices.map(i=>taxonGridColor(data,palette,i)),line:{color:colors().surface,width:1.5}},
            insidetextorientation:'radial',hovertemplate:'%{label}<br>%{value:,} fragments<extra></extra>'};
        drawPlot(el,[trace],layout({height:480,margin:{l:0,r:0,t:0,b:0},showlegend:false,uniformtext:{minsize:10,mode:'hide'}}));
    }
    function showLevel(index) {
        currentLevel=index; setTabs('detail',index);
        const data=levelData(index), name=levelNames[index];
        putText('level-title',name); putText('level-subtitle',limits[index] || '');
        if (!data.samples.includes(currentSample)) currentSample=data.samples[0] || null;
        document.getElementById('taxon-grid').scrollTop=0;
        drawBarplot(); drawTaxonGrid(); loadSample(currentSample);
    }
    function moveSample(step, focusSelection=false) {
        if (currentLevel == null) return;
        const samples=focusSelection ? [...document.querySelectorAll('#sample-buttons .sample-button')].map(button => button.dataset.sample) :
            levelData(currentLevel).samples;
        const index=samples.indexOf(currentSample);
        const next=Math.max(0,Math.min(samples.length-1,index+step));
        if (next>=0 && next!==index) {
            loadSample(samples[next]);
            if (focusSelection) {
                const button=document.querySelector('#sample-buttons .sample-button.active');
                button?.focus(); button?.scrollIntoView({block:'nearest',inline:'nearest'});
            }
        }
    }
    document.querySelectorAll('.level-tab').forEach(tab => tab.addEventListener('click',()=>
        tab.dataset.level === 'overview' ? showOverview() : showLevel(Number(tab.dataset.level))));
    document.getElementById('level-tabs')?.addEventListener('keydown',event=>{
        const tabs=[...document.querySelectorAll('.level-tab')];
        const at=tabs.indexOf(document.activeElement);
        if (at < 0 || !['ArrowLeft','ArrowRight','Home','End'].includes(event.key)) return;
        event.preventDefault();
        const next=event.key === 'Home' ? 0 : event.key === 'End' ? tabs.length-1 :
            (at+(event.key === 'ArrowRight' ? 1 : -1)+tabs.length)%tabs.length;
        tabs[next].focus(); tabs[next].click();
    });
    document.getElementById('tab-absolute').addEventListener('click',()=>{currentBarplot='absolute';document.getElementById('tab-absolute').setAttribute('aria-pressed','true');document.getElementById('tab-relative').setAttribute('aria-pressed','false');drawBarplot();});
    document.getElementById('tab-relative').addEventListener('click',()=>{currentBarplot='relative';document.getElementById('tab-absolute').setAttribute('aria-pressed','false');document.getElementById('tab-relative').setAttribute('aria-pressed','true');drawBarplot();});
    document.getElementById('composition-sort').addEventListener('change',drawBarplot);
    document.getElementById('download-composition-csv').addEventListener('click',()=>{
        const data=levelData(currentLevel), order=compositionOrder(data);
        downloadCsv(levelNames[currentLevel],'composition',[
            ['Sample',...data.composition.map(item=>item.name)],
            ...order.map(index=>[data.samples[index],...data.composition.map(item=>item.values[index])])]);
    });
    document.getElementById('taxon-grid').addEventListener('scroll',()=>{
        if (currentLevel==null) return;
        const tree=ensureTree(levelData(currentLevel)), scroller=document.getElementById('taxon-grid');
        if (!windowCovers(tree.window,scroller.scrollTop,scroller.clientHeight,tree.visible.length)) renderTaxonRows();
    });
    document.getElementById('taxon-search').addEventListener('input',()=>{
        document.getElementById('taxon-grid').scrollTop=0; drawTaxonGrid();
    });
    document.getElementById('taxon-sort').addEventListener('change',drawTaxonGrid);
    document.getElementById('taxon-depth').addEventListener('change',event=>setTaxonDepth(Number(event.target.value)));
    document.getElementById('taxon-expand').addEventListener('click',()=>setTaxonDepth(7));
    document.getElementById('taxon-collapse').addEventListener('click',()=>setTaxonDepth(0));
    document.getElementById('taxon-jump-sunburst').addEventListener('click',()=>
        document.getElementById('sunburst-card').scrollIntoView({block:'start'}));
    for (const mode of ['bubbles','heatmap','fragments','share','presence']) {
        document.getElementById('taxon-'+mode).addEventListener('click',()=>{
            const group=['bubbles','heatmap'].includes(mode) ? ['bubbles','heatmap'] : ['fragments','share','presence'];
            group.forEach(item=>document.getElementById('taxon-'+item).setAttribute('aria-pressed',String(item===mode)));
            drawTaxonGrid();
        });
    }
    document.getElementById('taxon-grid-body').addEventListener('keydown',event=>{
        if (!['ArrowUp','ArrowDown','ArrowLeft','ArrowRight'].includes(event.key) ||
            event.ctrlKey || event.metaKey || event.altKey || event.shiftKey || !event.target.matches('.taxon-node')) return;
        event.preventDefault();
        const tree=ensureTree(levelData(currentLevel)), index=Number(event.target.dataset.taxon);
        if (event.key==='ArrowUp' || event.key==='ArrowDown') {
            const position=tree.visible.indexOf(index)+(event.key==='ArrowUp' ? -1 : 1);
            if (position>=0 && position<tree.visible.length) focusTaxon(tree.visible[position]);
        } else if (tree.children[index].length) {
            clearTaxonSearch();
            changeTaxonExpansion(levelData(currentLevel),tree,index,event.key==='ArrowRight');
            drawTaxonGrid(); focusTaxon(index);
        }
    });
    document.getElementById('download-taxon-visible-csv').addEventListener('click',()=>
        downloadCsv(levelNames[currentLevel],'visible_taxa',taxonCsvRows(false)));
    document.getElementById('download-taxon-all-csv').addEventListener('click',()=>
        downloadCsv(levelNames[currentLevel],'all_taxa',taxonCsvRows(true)));
    document.getElementById('sample-search').addEventListener('input',renderNavigator);
    document.getElementById('sample-prev').addEventListener('click',()=>moveSample(-1));
    document.getElementById('sample-next').addEventListener('click',()=>moveSample(1));
    document.addEventListener('keydown',event=>{
        if (currentPanel !== 'detail' || !['ArrowUp','ArrowDown'].includes(event.key)) return;
        if (event.ctrlKey || event.metaKey || event.altKey || event.shiftKey) return;
        if (!document.getElementById('sample-buttons').contains(event.target)) return;
        event.preventDefault(); moveSample(event.key === 'ArrowUp' ? -1 : 1, true);
    });
    window.addEventListener('resize',()=>requestAnimationFrame(()=>{
        if (currentPanel === 'detail') drawBarplot();
        resizeVisibleCharts();
    }));
    const themeButton=document.getElementById('theme-toggle');
    function redrawForTheme() {
        applyLevelAccents();
        themeButton.setAttribute('aria-label','Switch to '+(effectiveTheme() === 'dark' ? 'light' : 'dark')+' theme');
        if (currentPanel === 'overview') drawOverview();
        else { overviewDrawn=false; drawBarplot(); drawTaxonGrid(); loadSample(currentSample); }
    }
    themeButton.addEventListener('click',()=>{
        document.documentElement.setAttribute('data-theme',effectiveTheme() === 'dark' ? 'light' : 'dark');
    });
    new MutationObserver(redrawForTheme).observe(document.documentElement,{attributes:true,attributeFilter:['data-theme']});
    themePreference.addEventListener('change',()=>{
        if (!document.documentElement.hasAttribute('data-theme')) redrawForTheme();
    });
    if (overview) {
        document.getElementById('download-comparison-csv').addEventListener('click',()=>{
            const headers=['Level','Thresholds','Fragments retained','% of aligned fragments',
                'Classified fragments','% classified',...ranks.map(rank=>'Distinct '+rank.toLowerCase()+' taxa')];
            const rows=overviewLevels.map((item,index)=>{
                const total=item.totals, denominator=total.classified_fragments+total.unclassified_fragments;
                return [levelNames[index],limits[index] || 'n/a',total.fragments_retained,
                    total.retained_fraction==null ? 'n/a' : Number((100*total.retained_fraction).toFixed(1)),
                    total.classified_fragments,denominator ? Number((100*total.classified_fragments/denominator).toFixed(1)) : 'n/a',
                    ...total.distinct_taxa];
            });
            downloadCsv('overview','level_comparison',[headers,...rows]);
        });
        const sort=document.getElementById('dot-sort');
        levelNames.forEach((level,index)=>{const option=document.createElement('option');option.value=String(index);option.textContent=level;sort.append(option);});
        sort.addEventListener('change',()=>{dotSort=Number(sort.value);drawDots();});
        for (const measure of ['count','percent']) document.getElementById('dot-'+measure).addEventListener('click',()=>{
            dotMeasure=measure;
            document.getElementById('dot-count').setAttribute('aria-pressed',String(measure==='count'));
            document.getElementById('dot-percent').setAttribute('aria-pressed',String(measure==='percent'));
            drawDots();
        });
        document.getElementById('overview-search').addEventListener('input',()=>{renderSampleTable();drawDots();});
        document.getElementById('download-sample-csv').addEventListener('click',()=>{
            const rows=renderSampleTable();
            downloadCsv('overview','per_sample',[['Sample',...levelNames],...rows.map(row=>[row.name,...row.values])]);
        });
    }
    applyLevelAccents();
    themeButton.setAttribute('aria-label','Switch to '+(effectiveTheme() === 'dark' ? 'light' : 'dark')+' theme');
    if (overview) showOverview(); else showLevel(0);
})();
</script>'''


if __name__ == "__main__":
    main()
