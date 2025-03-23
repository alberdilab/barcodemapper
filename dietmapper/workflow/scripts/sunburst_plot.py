import pandas as pd
import argparse
import json
import colorsys


def parse_args():
    parser = argparse.ArgumentParser(description="Generate a static HTML sunburst plot with sidebar sample selection.")
    parser.add_argument("--input", required=True, help="Input TSV file with taxonomy and abundance per sample")
    parser.add_argument("--output", required=True, help="Output HTML file")
    return parser.parse_args()


def assign_color(kingdom, depth):
    if kingdom.startswith("k__Animalia"):
        base_hue = 0  # red
    elif kingdom.startswith("k__Fungi"):
        base_hue = 220  # blue
    elif kingdom.startswith("k__Viridiplantae"):
        base_hue = 120  # green
    else:
        base_hue = 0  # fallback gray

    lightness = 0.3 + (depth * 0.1)
    lightness = min(lightness, 0.85)
    rgb = colorsys.hls_to_rgb(base_hue / 360, lightness, 0.6)
    return '#%02x%02x%02x' % tuple(int(x * 255) for x in rgb)


def prepare_data(df):
    samples = df.columns[1:]
    all_data = {}

    for sample in samples:
        subset = df[df[sample] > 0].copy()
        subset = subset.rename(columns={sample: 'Abundance'})
        taxonomy_levels = subset['Taxonomy'].str.split(';', expand=True)
        max_depth = taxonomy_levels.shape[1]
        colnames = [f'Level_{i}' for i in range(max_depth)]
        taxonomy_levels.columns = colnames
        subset = pd.concat([subset, taxonomy_levels], axis=1)

        nodes = {}
        for _, row in subset.iterrows():
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
            # Add abundance only to leaf node
            nodes[path_id]["value"] += row["Abundance"]

        frame = pd.DataFrame(nodes.values())
        frame["color"] = [assign_color(row["kingdom"], row["depth"]) for _, row in frame.iterrows()]
        all_data[sample] = frame.to_dict(orient="list")

    return list(samples), all_data


def generate_html(samples, data_by_sample, output_file):
    menu_html = "".join([f'<div onclick="loadSample(\'{s}\')" class="sample-button">{s}</div>' for s in samples])

    with open(output_file, "w") as f:
        f.write(f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>DietMapper – Sunburst Plot</title>
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
        }}
        #sample-title {{
            font-size: 20px;
            margin-bottom: 10px;
        }}
        #sunburst {{
            width: 100%;
            height: 700px;
        }}
    </style>
</head>
<body>
    <div id="header">DietMapper</div>
    <div id="main">
        <div id="sidebar">
            <h3>Samples</h3>
            {menu_html}
        </div>
        <div id="plot">
            <div id="sample-title">{samples[0]}</div>
            <div id="sunburst"></div>
        </div>
    </div>
    <script>
        const dataBySample = {json.dumps(data_by_sample)};

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

        loadSample("{samples[0]}");
    </script>
</body>
</html>
""")


if __name__ == "__main__":
    args = parse_args()
    df = pd.read_csv(args.input, sep="\t")
    samples, data_by_sample = prepare_data(df)
    generate_html(samples, data_by_sample, args.output)
