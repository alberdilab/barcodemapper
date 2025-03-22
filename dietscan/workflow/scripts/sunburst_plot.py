import pandas as pd
import plotly.graph_objects as go
import argparse

# --- ARGPARSE ---
parser = argparse.ArgumentParser(description="Generate sunburst plot from taxonomic table.")
parser.add_argument("-i", "--input", required=True, help="Input TSV file with taxonomy and sample counts.")
parser.add_argument("-o", "--output", required=True, help="Output HTML file for sunburst plot.")
args = parser.parse_args()

# --- LOAD DATA ---
df = pd.read_csv(args.input, sep="\t")
df_long = df.melt(id_vars=["Taxonomy"], var_name="Sample", value_name="Count")

# Clean counts
df_long["Count"] = df_long["Count"].fillna(0).astype(int)
df_long = df_long[df_long["Count"] > 0]

# --- SUNBURST HELPER ---
def build_sunburst_data(sample_df):
    labels, parents, values = [], [], []
    seen = set()
    for _, row in sample_df.iterrows():
        tax_levels = row["Taxonomy"].split(";")
        count = row["Count"]
        for i, level in enumerate(tax_levels):
            if level not in seen:
                seen.add(level)
                labels.append(level)
                parents.append("" if i == 0 else tax_levels[i - 1])
                values.append(count if i == len(tax_levels) - 1 else 0)
    return labels, parents, values

# --- BUILD FIGURE ---
samples = df_long["Sample"].unique()
sunburst_traces = []
buttons = []

for i, sample in enumerate(samples):
    sample_df = df_long[df_long["Sample"] == sample]
    labels, parents, values = build_sunburst_data(sample_df)
    trace = go.Sunburst(
        labels=labels,
        parents=parents,
        values=values,
        branchvalues="total",
        visible=(i == 0)
    )
    sunburst_traces.append(trace)
    buttons.append({
        "label": sample,
        "method": "update",
        "args": [
            {"visible": [j == i for j in range(len(samples))]},
            {"title": f"Taxonomic Breakdown - {sample}"}
        ]
    })

fig = go.Figure(data=sunburst_traces)
fig.update_layout(
    title=f"Taxonomic Breakdown - {samples[0]}",
    updatemenus=[{
        "buttons": buttons,
        "direction": "down",
        "x": 1.1,
        "y": 1.15,
        "showactive": True
    }],
    margin=dict(t=50, l=0, r=0, b=0)
)

# --- SAVE HTML ---
fig.write_html(args.output)
