import pandas as pd
import argparse
import json
import colorsys
from functools import lru_cache
from html import escape
from datetime import datetime
from plotly.offline import get_plotlyjs

try:
    from barcodemapper.workflow.scripts.level_stats import load_merged, load_stats, overview_stats
except ModuleNotFoundError:
    from level_stats import load_merged, load_stats, overview_stats


RANK_NAMES = ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")


def parse_args():
    parser = argparse.ArgumentParser(description="Generate a static HTML sunburst plot with sidebar sample selection.")
    parser.add_argument("--input", required=True, action="append", help="NAME=PATH for each level, or PATH for one legacy level")
    parser.add_argument("--stats", action="append", default=[], help="Per-sample extraction stats JSON")
    parser.add_argument("--limits", action="append", default=[], metavar="NAME=MAXPCT:MINCOV",
                        help="Stringency thresholds for an overview level")
    parser.add_argument("--output", required=True, help="Output HTML file")
    parser.add_argument("--version", default="unknown", help="Version of BarcodeMapper")
    return parser.parse_args()


@lru_cache(maxsize=None)
def assign_color(kingdom, depth):
    if kingdom.startswith("k__Animalia"):
        base_hue = 0
    elif kingdom.startswith("k__Fungi"):
        base_hue = 220
    elif kingdom.startswith("k__Viridiplantae"):
        base_hue = 120
    else:
        base_hue = 0

    lightness = 0.3 + (depth * 0.1)
    lightness = min(lightness, 0.85)
    rgb = colorsys.hls_to_rgb(base_hue / 360, lightness, 0.6)
    return '#%02x%02x%02x' % tuple(int(x * 255) for x in rgb)


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


def _overview_table(overview, limits):
    levels = overview["levels"]
    headings = ["Level", "Max mismatch %", "Min aligned bases"]
    headings += ["Fragments retained", "% of aligned fragments", "Classified fragments",
                 "% classified"]
    headings += [f"Distinct {rank.lower()} taxa" for rank in RANK_NAMES]
    rows = []
    for name, item in levels.items():
        totals = item["totals"]
        classified = totals["classified_fragments"]
        denominator = classified + totals["unclassified_fragments"]
        values = [name]
        level_limits = limits.get(name)
        values += [f"{level_limits[0]:g}" if level_limits else "n/a",
                   str(level_limits[1]) if level_limits else "n/a"]
        values += [str(totals["fragments_retained"]) if totals["fragments_retained"] is not None else "n/a",
                   _percent(totals["retained_fraction"]), str(classified),
                   _percent(classified / denominator if denominator else None)]
        values += [str(count) for count in totals["distinct_taxa"]]
        rows.append("<tr>" + "".join(f"<td>{escape(value)}</td>" for value in values) + "</tr>")
    return ("<table><thead><tr>" + "".join(f"<th>{escape(heading)}</th>" for heading in headings)
            + "</tr></thead><tbody>" + "".join(rows) + "</tbody></table>")


def _write_html(level_data, overview, output_file, version, limits=None):
    multiple = len(level_data) > 1
    names = list(level_data)
    first_name = names[0]
    first_samples = level_data[first_name]["samples"]
    menu_html = "".join(
        f'<button type="button" class="sample-button" data-sample="{escape(sample, quote=True)}">{escape(sample)}</button>'
        for sample in first_samples
    )
    tabs = ""
    if multiple:
        tabs = '<nav id="level-tabs" aria-label="Stringency levels">'
        tabs += '<button type="button" class="level-tab active" id="tab-overview" data-level="overview">Overview</button>'
        tabs += "".join(
            f'<button type="button" class="level-tab" id="tab-level-{index}" data-level="{index}">{escape(name)}</button>'
            for index, name in enumerate(names)
        )
        tabs += "</nav>"
    overview_panel = ""
    if multiple:
        overview_panel = f'''<section id="overview-panel">
            <h2>Overview</h2>
            <div class="table-wrap">{_overview_table(overview, limits or {})}</div>
            <h3>Classified fragments per sample</h3><div id="classified-chart" class="overview-chart"></div>
            <h3>Share resolved to each rank</h3><div id="depth-chart" class="overview-chart"></div>
            <h3>Kingdom composition by level</h3><div id="kingdom-chart" class="overview-chart"></div>
            <p class="footnote">Ultrastrict requires at least 120 aligned bases while the default minimum read length is 100, so some reads can never pass that level. Stricter levels can lower the assigned rank; LCA over fewer alignments can also make assignments deeper. Counts are not simple subsets at every rank.</p>
        </section>'''
    generation_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(output_file, "w", encoding="utf-8") as f:
        f.write(f'''<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>BarcodeMapper – Report</title>
    <script>{get_plotlyjs()}</script>
    <style>
        body {{ margin: 0; font-family: Arial, sans-serif; display: flex; height: 100vh; flex-direction: column; }}
        #header {{ background-color: #2c3e50; color: white; padding: 12px 20px; font-size: 22px; font-weight: bold; display: flex; justify-content: space-between; align-items: center; }}
        #header-link {{ text-decoration: none; color: white; }}
        #header-link:hover {{ color: #ecf0f1; }}
        #version {{ font-size: 18px; opacity: 0.8; }}
        #date {{ font-size: 14px; opacity: 0.7; }}
        #level-tabs {{ display: flex; gap: 5px; padding: 8px 12px; background: #e9edf0; overflow-x: auto; }}
        .level-tab {{ border: 0; border-radius: 16px; padding: 8px 16px; cursor: pointer; white-space: nowrap; background: #ddd; }}
        .level-tab.active {{ background: #b6dbfc; }}
        #main {{ display: flex; flex: 1; overflow: hidden; }}
        #sidebar {{ width: 220px; background: #f5f5f5; padding: 10px; overflow-y: auto; border-right: 1px solid #ccc; }}
        .sample-button {{ padding: 10px; cursor: pointer; border-radius: 4px; margin-bottom: 5px; background-color: #fff; border: 0; display: block; text-align: left; width: 100%; }}
        .sample-button:hover {{ background-color: #e0e0e0; }}
        .sample-button.active {{ background-color: #b6dbfc; }}
        #plot {{ flex-grow: 1; padding: 10px; overflow-y: auto; background-color: #f0f0f0; }}
        #barplot-tabs {{ margin-bottom: 10px; }}
        .tab {{ display: inline-block; padding: 8px 16px; cursor: pointer; background: #ddd; margin-right: 5px; border-radius: 16px; border: 0; }}
        .tab.active {{ background: #b6dbfc; }}
        #sample-title {{ font-size: 24px; font-weight: bold; margin-top: 20px; margin-bottom: 10px; }}
        #barplot {{ width: 100%; height: 400px; }}
        #sunburst {{ width: 100%; height: 700px; }}
        #overview-panel h2 {{ margin-top: 8px; }}
        .table-wrap {{ overflow-x: auto; }}
        table {{ border-collapse: collapse; background: white; white-space: nowrap; }}
        th, td {{ border: 1px solid #ccc; padding: 7px 10px; text-align: right; }}
        th:first-child, td:first-child {{ text-align: left; }}
        .overview-chart {{ width: 100%; height: 380px; }}
        .footnote {{ max-width: 900px; font-size: 13px; color: #555; }}
    </style>
</head>
<body>
    <div id="header">
        <div><a href="https://github.com/alberdilab/barcodemapper" target="_blank" rel="noopener noreferrer" id="header-link">BarcodeMapper</a> <span id="version">v{escape(version)}</span></div>
        <div id="date">{generation_date}</div>
    </div>
    {tabs}
    <div id="main">
        <div id="sidebar"{' style="display:none"' if multiple else ''}><h3>Samples</h3><div id="sample-buttons">{menu_html}</div></div>
        <div id="plot">
            {overview_panel}
            <section id="detail-panel"{' style="display:none"' if multiple else ''}>
                <div id="barplot-tabs">
                    <button type="button" class="tab active" id="tab-absolute">Marker-matching fragments</button>
                    <button type="button" class="tab" id="tab-relative">Share of marker-matching fragments</button>
                </div>
                <div id="barplot"></div>
                <div id="sample-title">{escape(first_samples[0]) if first_samples else ''}</div>
                <div id="sunburst"></div>
            </section>
        </div>
    </div>
''')
        for index, data in enumerate(level_data.values()):
            f.write(f'<script type="application/json" id="data-level-{index}">{safe_json(data)}</script>\n')
        if multiple:
            f.write(f'<script type="application/json" id="data-overview">{safe_json(overview)}</script>\n')
        f.write(REPORT_SCRIPT)
        f.write("\n</body>\n</html>\n")


def generate_html(samples, data_by_sample, output_file, abs_barplot, rel_barplot, version):
    """Keep the existing single-level Python entry point for callers."""
    data = {"samples": samples, "dataBySample": data_by_sample,
            "absData": abs_barplot.to_dict(orient="list"),
            "relData": rel_barplot.to_dict(orient="list")}
    _write_html({"custom": data}, None, output_file, version)


def generate_report(merged_by_level, extract_stats, output_file, version="unknown", limits=None):
    level_data = {}
    for name, frame in merged_by_level.items():
        samples, sunbursts, absolute, relative = prepare_data(frame)
        level_data[name] = {"samples": samples, "dataBySample": sunbursts,
                            "absData": absolute.to_dict(orient="list"),
                            "relData": relative.to_dict(orient="list")}
    if not level_data:
        raise ValueError("At least one input level is required")
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


REPORT_SCRIPT = '''<script>
    const levelCache = {};
    let currentLevel = null;
    let currentSample = null;
    let currentBarplot = "absolute";
    let overviewDrawn = false;

    function levelData(index) {
        if (!levelCache[index]) {
            levelCache[index] = JSON.parse(document.getElementById('data-level-' + index).textContent);
        }
        return levelCache[index];
    }

    function drawBarplot(type) {
        const data = type === "absolute" ? levelData(currentLevel).absData : levelData(currentLevel).relData;
        const yLabel = type === "absolute" ? "Number of fragments" : "Fraction of marker-matching fragments";
        const phyla = data.Phylum;
        const kingdoms = data.Kingdom;
        const colors = data.color;
        const samples = Object.keys(data).filter(k => !["Phylum", "Kingdom", "color"].includes(k));
        const traces = phyla.map((phylum, idx) => {
            const label = phylum + " (" + kingdoms[idx] + ")";
            return {
                x: samples,
                y: samples.map(s => data[s][idx]),
                type: 'bar',
                name: label,
                customdata: samples.map(() => label),
                marker: { color: colors[idx], line: { color: '#ffffff', width: 1 } },
                hovertemplate: '<b>%{customdata}</b><br>Value: %{y}<extra></extra>'
            };
        });
        Plotly.newPlot('barplot', traces, {
            barmode: 'stack', xaxis: { title: '' }, yaxis: { title: yLabel }, margin: { t: 10 },
            plot_bgcolor: '#f0f0f0', paper_bgcolor: '#f0f0f0',
            legend: { traceorder: 'normal', y: 0.925, yanchor: 'top' }
        });
    }

    function switchBarplot(type) {
        currentBarplot = type;
        document.querySelectorAll('#barplot-tabs .tab').forEach(t => t.classList.remove('active'));
        document.getElementById('tab-' + type).classList.add('active');
        drawBarplot(type);
    }

    function loadSample(sample) {
        currentSample = sample;
        document.getElementById('sample-title').textContent = sample || '';
        document.querySelectorAll('.sample-button').forEach(button => {
            button.classList.toggle('active', button.dataset.sample === sample);
        });
        const d = levelData(currentLevel).dataBySample[sample];
        if (!d || d.label.length === 0) {
            Plotly.purge('sunburst');
            document.getElementById('sunburst').innerHTML = '<p>No data for this sample.</p>';
            return;
        }
        Plotly.newPlot('sunburst', [{
            type: 'sunburst', ids: d.id, labels: d.label, parents: d.parent, values: d.value,
            marker: { colors: d.color }, branchvalues: 'total', hoverinfo: 'label+value+percent entry'
        }], { margin: { t: 0, l: 0, r: 0, b: 0 } });
    }

    function showLevel(index) {
        const data = levelData(index);
        currentLevel = index;
        document.getElementById('overview-panel')?.style.setProperty('display', 'none');
        document.getElementById('detail-panel').style.display = '';
        document.getElementById('sidebar').style.display = '';
        document.querySelectorAll('.level-tab').forEach(tab => tab.classList.remove('active'));
        document.getElementById('tab-level-' + index)?.classList.add('active');
        const buttons = document.getElementById('sample-buttons');
        buttons.replaceChildren();
        for (const sample of data.samples) {
            const button = document.createElement('button');
            button.type = 'button';
            button.className = 'sample-button';
            button.dataset.sample = sample;
            button.textContent = sample;
            button.addEventListener('click', () => loadSample(sample));
            buttons.appendChild(button);
        }
        if (!data.samples.includes(currentSample)) currentSample = data.samples[0] || null;
        drawBarplot(currentBarplot);
        loadSample(currentSample);
    }

    function drawOverview() {
        const overview = JSON.parse(document.getElementById('data-overview').textContent);
        const levels = Object.keys(overview.levels);
        const samples = [...new Set(overview.classified_by_sample.map(row => row.sample))];
        const classified = levels.map(level => ({
            type: 'bar', name: level, x: samples,
            y: samples.map(sample => overview.levels[level].samples[sample]?.classified_fragments || 0)
        }));
        Plotly.newPlot('classified-chart', classified, {
            barmode: 'group', yaxis: { title: 'Classified fragments' }, margin: { t: 10 },
            plot_bgcolor: '#f0f0f0', paper_bgcolor: '#f0f0f0'
        });
        const ranks = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'];
        const depth = ranks.map((rank, index) => ({
            type: 'bar', name: rank, x: levels,
            y: levels.map(level => {
                const totals = overview.levels[level].totals;
                const fragments = totals.classified_fragments + totals.unclassified_fragments;
                return fragments ? totals.depth_distribution[index] / fragments : 0;
            })
        }));
        Plotly.newPlot('depth-chart', depth, {
            barmode: 'group', yaxis: { title: 'Share of fragments', tickformat: '.0%', range: [0, 1] },
            margin: { t: 10 }, plot_bgcolor: '#f0f0f0', paper_bgcolor: '#f0f0f0'
        });
        const kingdoms = [...new Set(levels.flatMap(level =>
            Object.keys(overview.levels[level].totals.kingdom_composition)))];
        const composition = kingdoms.map(kingdom => ({
            type: 'bar', name: kingdom, x: levels,
            y: levels.map(level => {
                const totals = overview.levels[level].totals;
                return totals.classified_fragments ?
                    (totals.kingdom_composition[kingdom] || 0) / totals.classified_fragments : 0;
            })
        }));
        Plotly.newPlot('kingdom-chart', composition, {
            barmode: 'stack', yaxis: { title: 'Share of classified fragments', tickformat: '.0%', range: [0, 1] },
            margin: { t: 10 }, plot_bgcolor: '#f0f0f0', paper_bgcolor: '#f0f0f0'
        });
    }

    function showOverview() {
        document.getElementById('overview-panel').style.display = '';
        document.getElementById('detail-panel').style.display = 'none';
        document.getElementById('sidebar').style.display = 'none';
        document.querySelectorAll('.level-tab').forEach(tab => tab.classList.remove('active'));
        document.getElementById('tab-overview').classList.add('active');
        if (!overviewDrawn) {
            drawOverview();
            overviewDrawn = true;
        }
    }

    document.getElementById('tab-absolute').addEventListener('click', () => switchBarplot('absolute'));
    document.getElementById('tab-relative').addEventListener('click', () => switchBarplot('relative'));
    document.querySelectorAll('.level-tab').forEach(tab => {
        tab.addEventListener('click', () => tab.dataset.level === 'overview' ?
            showOverview() : showLevel(Number(tab.dataset.level)));
    });
    if (document.getElementById('tab-overview')) showOverview();
    else showLevel(0);
</script>'''


if __name__ == "__main__":
    main()
