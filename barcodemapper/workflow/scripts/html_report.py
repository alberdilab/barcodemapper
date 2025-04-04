import pandas as pd
import argparse
import json
import colorsys
from datetime import datetime


def parse_args():
    parser = argparse.ArgumentParser(description="Generate a static HTML sunburst plot with sidebar sample selection.")
    parser.add_argument("--input", required=True, help="Input TSV file with taxonomy and abundance per sample")
    parser.add_argument("--output", required=True, help="Output HTML file")
    parser.add_argument("--version", default="1.0.0", help="Version of BarcodeMapper")
    return parser.parse_args()


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

    for sample in samples:
        subset = df.copy()
        subset = subset.rename(columns={sample: 'Abundance'})
        taxonomy_levels = subset['Taxonomy'].str.split(';', expand=True)
        max_depth = taxonomy_levels.shape[1]
        colnames = [f'Level_{i}' for i in range(max_depth)]
        taxonomy_levels.columns = colnames
        subset = pd.concat([subset, taxonomy_levels], axis=1)

        nodes = {}
        for _, row in subset.iterrows():
            if row["Abundance"] <= 0:
                continue
            levels = [lvl for lvl in row[colnames] if pd.notnull(lvl)]
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
                    nodes[path_id]["value"] += row["Abundance"]

        frame = pd.DataFrame(nodes.values())
        frame["color"] = [assign_color(row["kingdom"], row["depth"]) for _, row in frame.iterrows()]
        all_data[sample] = frame.to_dict(orient="list")

        phylum_rows = frame[frame["depth"] == 1][["label", "kingdom", "value", "color"]].copy()
        phylum_rows.columns = ["Phylum", "Kingdom", sample, "color"]
        all_phylum_rows.append(phylum_rows)

    merged = all_phylum_rows[0]
    for phylum_df in all_phylum_rows[1:]:
        merged = pd.merge(merged, phylum_df, on=["Phylum", "Kingdom", "color"], how="outer")
    merged = merged.fillna(0)
    merged = merged.sort_values(["Kingdom", "Phylum"])

    rel = merged.copy()
    for sample in samples:
        total = rel[sample].sum()
        rel[sample] = rel[sample] / total if total > 0 else 0

    return list(samples), all_data, merged, rel


def generate_html(samples, data_by_sample, output_file, abs_barplot, rel_barplot, version):
    menu_html = "".join([f'<div onclick="loadSample(\'{s}\')" class="sample-button">{s}</div>' for s in samples])
    abs_json = abs_barplot.to_dict(orient="list")
    rel_json = rel_barplot.to_dict(orient="list")
    generation_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    with open(output_file, "w") as f:
        f.write(f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>BarcodeMapper – Sunburst Plot</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        body {{
            margin: 0;
            font-family: Arial, sans-serif;
            display: flex;
            height: 100vh;
            flex-direction: column;
        }}
        #header {{
            background-color: #2c3e50;
            color: white;
            padding: 12px 20px;
            font-size: 22px;
            font-weight: bold;
            display: flex;
            justify-content: space-between;
            align-items: center;
        }}

	#header-link {{
    		text-decoration: none;
    		color: white; /* Or whatever color you want */
	}}

	#header-link:hover {{
  		  color: #ecf0f1; /* Optional: slightly lighter on hover */
	}}
        #version {{
            font-size: 18px;
            opacity: 0.8;
        }}
        #date {{
            font-size: 14px;
            opacity: 0.7;
        }}
        #main {{
            display: flex;
            flex: 1;
            overflow: hidden;
        }}
        #sidebar {{
            width: 220px;
            background: #f5f5f5;
            padding: 10px;
            overflow-y: auto;
            border-right: 1px solid #ccc;
        }}
        .sample-button {{
            padding: 10px;
            cursor: pointer;
            border-radius: 4px;
            margin-bottom: 5px;
            background-color: #fff;
        }}
        .sample-button:hover {{
            background-color: #e0e0e0;
        }}
        #plot {{
            flex-grow: 1;
            padding: 10px;
            overflow-y: auto;
            background-color: #f0f0f0;
        }}
        #barplot-tabs {{
            margin-bottom: 10px;
        }}
        .tab {{
            display: inline-block;
            padding: 8px 16px;
            cursor: pointer;
            background: #ddd;
            margin-right: 5px;
            border-radius: 16px;
        }}
        .tab.active {{
            background: #b6dbfc;
        }}
        #sample-title {{
            font-size: 24px;
            font-weight: bold;
            margin-top: 20px;
            margin-bottom: 10px;
        }}
        #barplot {{
            width: 100%;
            height: 400px;
        }}
        #sunburst {{
            width: 100%;
            height: 700px;
        }}
    </style>
</head>
<body>
    <div id="header">
        <div><a href="https://github.com/alberdilab/barcodemapper" target="_blank" rel="BarcodeMapper website" id="header-link">BarcodeMapper</a> <span id="version">v{version}</span></div>
        <div id="date">{generation_date}</div>
    </div>
    <div id="main">
        <div id="sidebar">
            <h3>Samples</h3>
            {menu_html}
        </div>
        <div id="plot">
            <div id="barplot-tabs">
		<div class="tab active" id="tab-absolute" onclick="switchBarplot('absolute')">Absolute number of reads</div>
		<div class="tab" id="tab-relative" onclick="switchBarplot('relative')">Relative number of reads</div>
            </div>
            <div id="barplot"></div>
            <div id="sample-title">{samples[0]}</div>
            <div id="sunburst"></div>
        </div>
    </div>
    <script>
        const dataBySample = {json.dumps(data_by_sample)};
        const absData = {json.dumps(abs_json)};
        const relData = {json.dumps(rel_json)};
        let currentBarplot = "absolute";

        function drawBarplot(type) {{
            const data = type === "absolute" ? absData : relData;
            const yLabel = type === "absolute" ? "Number of reads" : "Relative abundance";

            const phyla = data["Phylum"];
            const kingdoms = data["Kingdom"];
            const colors = data["color"];
            const samples = Object.keys(data).filter(k => !["Phylum", "Kingdom", "color"].includes(k));

            const traces = phyla.map((phylum, idx) => {{
                const label = phylum + " (" + kingdoms[idx] + ")";
                return {{
                    x: samples,
                    y: samples.map(s => data[s][idx]),
                    type: 'bar',
                    name: label,
                    customdata: samples.map(() => label),
                    marker: {{
                        color: colors[idx],
                        line: {{
                            color: '#ffffff',
                            width: 1
                        }}
                    }},
                    hovertemplate: '<b>%{{customdata}}</b><br>Value: %{{y}}<extra></extra>'
                }};
            }});

            const layout = {{
                barmode: 'stack',
                xaxis: {{ title: '' }},
                yaxis: {{ title: yLabel }},
                margin: {{ t: 10 }},
                plot_bgcolor: '#f0f0f0',
                paper_bgcolor: '#f0f0f0',
                legend: {{
                    traceorder: 'normal',
                    y: 0.925,
                    yanchor: 'top'
                }}
            }};

            Plotly.newPlot('barplot', traces, layout);
        }}

        function switchBarplot(type) {{
            currentBarplot = type;
            document.querySelectorAll('.tab').forEach(t => t.classList.remove('active'));
            document.getElementById('tab-' + type).classList.add('active');
            drawBarplot(type);
        }}

        function loadSample(sample) {{
            document.getElementById("sample-title").innerText = sample;

            const d = dataBySample[sample];
            if (!d || d.label.length === 0) {{
                document.getElementById("sunburst").innerHTML = "<p>No data for this sample.</p>";
                return;
            }}

            const trace = {{
                type: "sunburst",
                ids: d.id,
                labels: d.label,
                parents: d.parent,
                values: d.value,
                marker: {{ colors: d.color }},
                branchvalues: "total",
                hoverinfo: "label+value+percent entry"
            }};

            const layout = {{
                margin: {{ t: 0, l: 0, r: 0, b: 0 }}
            }};

            Plotly.newPlot("sunburst", [trace], layout);
        }}

        drawBarplot(currentBarplot);
        loadSample("{samples[0]}");
    </script>
</body>
</html>
""")


if __name__ == "__main__":
    args = parse_args()
    df = pd.read_csv(args.input, sep="\t")
    samples, data_by_sample, abs_barplot, rel_barplot = prepare_data(df)
    generate_html(samples, data_by_sample, args.output, abs_barplot, rel_barplot, args.version)
