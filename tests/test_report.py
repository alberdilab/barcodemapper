"""Checks for the multi-level static HTML report."""

import importlib.util
import json
import re
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd


SCRIPT = Path(__file__).resolve().parents[1] / "barcodemapper/workflow/scripts/html_report.py"


def load_report():
    spec = importlib.util.spec_from_file_location("html_report", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def frame(samples=("alpha", "beta"), empty=False):
    if empty:
        return pd.DataFrame({"Taxonomy": ["Unclassified"], **{sample: [0] for sample in samples}})
    return pd.DataFrame({
        "Taxonomy": ["Unclassified", "k__Fungi", "k__Fungi;p__Ascomycota", "k__Animalia"],
        **{sample: [2, 8, 5, 3] for sample in samples},
    })


def test_cli_single_and_multi_level_tabs(tmp_path):
    for count in (1, 2, 4):
        paths = []
        for index in range(count):
            path = tmp_path / f"input-{index}.tsv"
            frame().to_csv(path, sep="\t", index=False)
            paths.append(path)
        output = tmp_path / f"report-{count}.html"
        args = [sys.executable, str(SCRIPT)]
        for index, path in enumerate(paths):
            args += ["--input", str(path) if count == 1 else f"level {index}={path}"]
            args += ["--limits", f"{'custom' if count == 1 else f'level {index}'}={index + 1}:90"]
        args += ["--output", str(output), "--version", "1.2.3"]
        subprocess.run(args, check=True)
        html = output.read_text()
        assert html.count('id="data-level-') == count
        assert html.count('id="tab-level-') == (count if count > 1 else 0)
        assert ('id="tab-overview"' in html) == (count > 1)
        assert ('id="overview-panel"' in html) == (count > 1)
        for index in range(count):
            assert f'id="data-level-{index}"' in html
            if count > 1:
                assert f'id="tab-level-{index}"' in html
                assert f"level {index}</button>" in html
        assert html.count('id="barplot"') == 1
        assert html.count('id="sunburst"') == 1
        assert "v1.2.3" in html


def test_overview_stats_and_empty_level(tmp_path):
    report = load_report()
    output = tmp_path / "empty.html"
    stats = [{"sample": "alpha", "aligned_fragments": 20, "levels": {
        "ultrastrict": {"fragments": 0}, "strict": {"fragments": 12}
    }}]
    report.generate_report({"ultrastrict": frame(empty=True), "strict": frame()}, stats, output,
                           limits={"ultrastrict": (0.0, 120), "strict": (2.0, 100)})
    html = output.read_text()
    assert "Max mismatch %" in html
    assert "Min aligned bases" in html
    assert "Fragments retained" in html
    assert "% of aligned fragments" in html
    assert "Distinct species taxa" in html
    assert "120" in html
    assert "n/a" in html
    assert "default minimum read length is 100" in html
    assert "Counts are not simple subsets at every rank" in html
    assert "No data for this sample" in html
    assert 'id="classified-chart"' in html
    assert 'id="depth-chart"' in html
    assert 'id="kingdom-chart"' in html
    assert '"retained_fraction": 0.6' in html


def test_cli_custom_limits_and_missing_limits(tmp_path):
    input_path = tmp_path / "input.tsv"
    frame().to_csv(input_path, sep="\t", index=False)
    output = tmp_path / "report.html"
    subprocess.run([sys.executable, str(SCRIPT),
                    "--input", f"custom={input_path}", "--input", f"other={input_path}",
                    "--limits", "custom=1.5:95", "--output", str(output)], check=True)
    html = output.read_text()
    assert "<tr><td>custom</td><td>1.5</td><td>95</td>" in html
    assert "<tr><td>other</td><td>n/a</td><td>n/a</td>" in html


def test_cli_loads_stats_and_reserves_overview_name(tmp_path):
    input_path = tmp_path / "input.tsv"
    frame().to_csv(input_path, sep="\t", index=False)
    stats_path = tmp_path / "stats.json"
    stats_path.write_text(json.dumps({
        "sample": "alpha", "aligned_fragments": 20,
        "levels": {"strict": {"fragments": 12}, "standard": {"fragments": 15}},
    }))
    output = tmp_path / "report.html"
    command = [sys.executable, str(SCRIPT), "--input", f"strict={input_path}",
               "--input", f"standard={input_path}", "--stats", str(stats_path),
               "--limits", "strict=2:100", "--limits", "standard=3:80",
               "--output", str(output)]
    subprocess.run(command, check=True)
    assert '"retained_fraction": 0.6' in output.read_text()
    invalid = subprocess.run([sys.executable, str(SCRIPT), "--input", f"overview={input_path}",
                              "--output", str(output)], capture_output=True, text=True)
    assert invalid.returncode != 0
    assert "reserved for the report tab" in invalid.stderr


def test_names_are_escaped_and_dom_ids_are_indexed(tmp_path):
    report = load_report()
    level = 'overview . "</script><script>alert(1)</script>'
    sample = 'sample . "</script><script>alert(2)</script>'
    table = frame(samples=(sample,))
    output = tmp_path / "escaped.html"
    report.generate_report({level: table, "second.level": table}, [], output, version='"<script>')
    html = output.read_text()
    assert level not in html
    assert sample not in html
    assert "&lt;/script&gt;" in html
    assert "\\u003c/script\\u003e" in html
    assert 'id="tab-level-0"' in html
    assert 'id="tab-level-1"' in html
    assert not re.search(r'id="[^"]*(?:overview \.|second\.level)[^"]*"', html)
    assert "v&quot;&lt;script&gt;" in html


def reference_prepare_data(df, assign_color):
    """Original row-wise report preparation for comparison."""
    samples = df.columns[1:]
    all_data = {}
    all_phylum_rows = []
    for sample in samples:
        subset = df[["Taxonomy", sample]].rename(columns={sample: "Abundance"}).copy()
        taxonomy_levels = subset["Taxonomy"].str.split(";", expand=True)
        colnames = [f"Level_{i}" for i in range(taxonomy_levels.shape[1])]
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
                    nodes[path_id] = {"id": path_id, "label": label, "parent": parent_id,
                                      "depth": depth, "value": 0, "kingdom": levels[0]}
                if depth == len(levels) - 1:
                    nodes[path_id]["value"] += row["Abundance"]
        node_frame = pd.DataFrame(nodes.values(), columns=["id", "label", "parent", "depth", "value", "kingdom"])
        node_frame["color"] = [assign_color(row["kingdom"], row["depth"])
                               for _, row in node_frame.iterrows()]
        all_data[sample] = node_frame.to_dict(orient="list")
        if node_frame.empty:
            phylum_rows = pd.DataFrame({
                "Phylum": pd.Series(dtype="string"), "Kingdom": pd.Series(dtype="string"),
                sample: pd.Series(dtype="int64"), "color": pd.Series(dtype="string")
            })
        else:
            phylum_rows = node_frame[(node_frame["depth"] == 1) | (node_frame["label"] == "Unclassified")][
                ["label", "kingdom", "value", "color"]].copy()
            phylum_rows.columns = ["Phylum", "Kingdom", sample, "color"]
            phylum_rows = phylum_rows.reset_index(drop=True)
            for _, kingdom in node_frame[node_frame["depth"] == 0].iterrows():
                if kingdom["label"] == "Unclassified":
                    continue
                children = node_frame[node_frame["parent"] == kingdom["id"]]
                unresolved = kingdom["value"] - children["value"].sum()
                if unresolved > 0:
                    phylum_rows.loc[len(phylum_rows)] = [
                        "Unresolved within " + kingdom["label"], kingdom["label"],
                        unresolved, kingdom["color"]]
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


def test_prepare_data_matches_original():
    report = load_report()
    table = pd.DataFrame({
        "Taxonomy": ["Unclassified", "k__Fungi", "k__Fungi;p__Ascomycota",
                     "k__Fungi;p__Ascomycota;c__Eurotiomycetes", "k__Animalia", "k__Animalia;p__Chordata"],
        "alpha": [2, 8, 6, 3, 5, 4], "beta": [0, 2, 0, 0, 1, 0], "zero": [0] * 6,
    })
    actual = report.prepare_data(table)
    expected = reference_prepare_data(table, report.assign_color)
    assert actual[0] == expected[0]
    assert actual[1] == expected[1]
    pd.testing.assert_frame_equal(actual[2], expected[2])
    pd.testing.assert_frame_equal(actual[3], expected[3])


def test_zero_sample_report(tmp_path):
    report = load_report()
    table = pd.DataFrame(columns=["Taxonomy"])
    output = tmp_path / "zero-samples.html"
    report.generate_report({"strict": table, "standard": table}, [], output)
    html = output.read_text()
    assert 'id="tab-overview"' in html
    assert '"samples": []' in html
    assert "No data for this sample" in html


def test_report_performance_smoke(tmp_path):
    report = load_report()
    taxonomy = ["Unclassified", "k__Fungi"]
    taxonomy += [f"k__Fungi;p__P{i};c__C{i};o__O{i}" for i in range(350)]
    table = pd.DataFrame({"Taxonomy": taxonomy,
                          **{f"sample {i}": [0, 350] + [1] * 350 for i in range(4)}})
    start = time.monotonic()
    report.generate_report({name: table for name in ("ultrastrict", "strict", "standard", "permissive")},
                           [], tmp_path / "large.html")
    assert time.monotonic() - start < 10
