"""Checks for the multi-level static HTML report."""

import importlib.util
import json
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "barcodemapper/workflow/scripts/html_report.py"


def load_report():
    spec = importlib.util.spec_from_file_location("html_report", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def embedded_json(html, name):
    match = re.search(rf'<script type="application/json" id="{re.escape(name)}">(.*?)</script>', html)
    assert match is not None
    return json.loads(match.group(1))


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
        assert len(re.findall(r'id="data-level-\d+"', html)) == count
        assert html.count('id="tab-level-') == (count if count > 1 else 0)
        assert ('id="tab-overview"' in html) == (count > 1)
        assert ('id="overview-panel"' in html) == (count > 1)
        for index in range(count):
            assert f'id="data-level-{index}"' in html
            if count > 1:
                assert f'id="tab-level-{index}"' in html
                assert f"<strong>level {index}</strong>" in html
                assert f"≤{index + 1}% edits · ≥90 bp" in html
        assert html.count('id="barplot"') == 1
        assert html.count('id="sunburst"') == 1
        assert "v1.2.3" in html
        assert 'id="theme-toggle"' in html
        assert ('id="download-sample-csv"' in html) == (count > 1)
        assert ('id="download-comparison-csv"' in html) == (count > 1)
        assert 'id="download-composition-csv"' in html
        assert 'id="download-taxon-visible-csv"' in html
        assert 'id="download-taxon-all-csv"' in html
        assert 'id="taxon-grid"' in html
        assert 'id="taxon-selection-note"' in html
        assert 'id="taxon-jump-sunburst"' in html
        assert 'id="sunburst-card"' in html
        assert ('id="sample-table"' in html) == (count > 1)
        assert 'role="tabpanel"' in html


def test_overview_stats_and_empty_level(tmp_path):
    report = load_report()
    output = tmp_path / "empty.html"
    stats = [{"sample": "alpha", "aligned_fragments": 20, "levels": {
        "ultrastrict": {"fragments": 0}, "strict": {"fragments": 12}
    }}]
    report.generate_report({"ultrastrict": frame(empty=True), "strict": frame()}, stats, output,
                           limits={"ultrastrict": (0.0, 120), "strict": (2.0, 100)})
    html = output.read_text()
    assert "Thresholds" in html
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
    level = embedded_json(html, "data-level-0")
    assert all(not pairs for pairs in level["taxa"]["values"])
    assert level["composition"] == []


def test_cli_custom_limits_and_missing_limits(tmp_path):
    input_path = tmp_path / "input.tsv"
    frame().to_csv(input_path, sep="\t", index=False)
    output = tmp_path / "report.html"
    subprocess.run([sys.executable, str(SCRIPT),
                    "--input", f"custom={input_path}", "--input", f"other={input_path}",
                    "--limits", "custom=1.5:95", "--output", str(output)], check=True)
    html = output.read_text()
    assert '<strong>custom</strong><small>≤1.5% edits · ≥95 bp</small>' in html
    assert '<th scope="row">custom</th><td>≤1.5% edits · ≥95 bp</td>' in html
    assert '<th scope="row">other</th><td>n/a</td>' in html


def test_named_single_level_uses_its_name_and_limits(tmp_path):
    report = load_report()
    output = tmp_path / "single.html"
    report.generate_report({"strict": frame()}, [], output, limits={"strict": (2.0, 100)})
    html = output.read_text()
    assert 'id="tab-overview"' not in html
    assert '<h1 id="level-title">strict</h1>' in html
    assert embedded_json(html, "data-level-order") == ["strict"]
    assert embedded_json(html, "data-limits") == ["≤2% edits · ≥100 bp"]


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
                              "--output", str(output)], capture_output=True, text=True, check=False)
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
                          **{f"sample {i}": [0, 350] + [1] * 350 for i in range(64)}})
    start = time.monotonic()
    report.generate_report({name: table for name in ("ultrastrict", "strict", "standard", "permissive")},
                           [], tmp_path / "large.html")
    assert time.monotonic() - start < 10


def test_composition_folds_to_top_eight_phyla(tmp_path):
    report = load_report()
    taxonomy = ["Unclassified", "k__Animalia"] + [f"k__Animalia;p__P{i}" for i in range(12)]
    table = pd.DataFrame({"Taxonomy": taxonomy,
                          "sample A": [0, 78] + list(range(12, 0, -1)),
                          "sample B": [0, 12] + [1] * 12})
    output = tmp_path / "folded.html"
    report.generate_report({"strict": table}, [], output)
    html = output.read_text()
    data = embedded_json(html, "data-level-0")
    items = data["composition"]
    assert len(items) == 9
    assert [item["name"] for item in items[:8]] == [f"p__P{i}" for i in range(8)]
    assert items[-1]["name"] == "Other Animalia"
    assert items[-1]["values"] == [10, 4]
    assert sum(item["values"][0] for item in items) == 78
    assert 'id="composition-table-view"' in html
    assert 'title="Download the folded composition categories shown in this table"' in html


def test_theme_level_colours_and_html_composition_legend(tmp_path):
    report = load_report()
    output = tmp_path / "themed.html"
    report.generate_report({name: frame() for name in ("ultrastrict", "strict", "standard", "permissive")},
                           [], output)
    html = output.read_text()
    palettes = embedded_json(html, "data-level-colors")
    assert set(palettes) == {"light", "dark"}
    assert all(len(values) == len(set(values)) == 4 for values in palettes.values())
    assert 'id="composition-legend" class="composition-legend" role="group"' in html


def test_folding_conserves_every_sample_at_every_level(tmp_path):
    report = load_report()
    samples = ("alpha", "beta", "zero")
    taxonomy = ["Unclassified", "k__Animalia"] + [f"k__Animalia;p__P{i}" for i in range(12)]
    levels = {}
    for multiplier, name in enumerate(("ultrastrict", "strict", "standard", "permissive"), start=1):
        columns = {"Taxonomy": taxonomy}
        for sample_index, sample in enumerate(samples):
            counts = [multiplier * (sample_index + 1) * (index + 1) for index in range(12)]
            if sample == "zero":
                counts = [0] * 12
            unresolved = 5 * multiplier if sample != "zero" else 0
            unclassified = 3 * multiplier if sample != "zero" else 0
            columns[sample] = [unclassified, sum(counts) + unresolved, *counts]
        levels[name] = pd.DataFrame(columns)
    output = tmp_path / "folded-all.html"
    report.generate_report(levels, [], output)
    html = output.read_text()
    assert embedded_json(html, "data-level-order") == list(levels)
    for level_index in range(len(levels)):
        data = embedded_json(html, f"data-level-{level_index}")
        assert data["samples"] == list(samples)
        assert any(item["name"] == "Other Animalia" for item in data["composition"])
        assert any(item["name"] == "Unresolved within k__Animalia" for item in data["composition"])
        for sample_index, sample in enumerate(samples):
            depth_zero = levels[list(levels)[level_index]].set_index("Taxonomy")
            expected = int(depth_zero.loc[["Unclassified", "k__Animalia"], sample].sum())
            assert data["totals"][sample_index] == expected
            assert sum(item["values"][sample_index] for item in data["composition"]) == expected


def test_numeric_level_order_and_sample_names(tmp_path):
    report = load_report()
    samples = ("1", "2", "10", "2021-01")
    table = frame(samples=samples)
    levels = {name: table.copy() for name in ("10", "9", "2")}
    output = tmp_path / "numeric-names.html"
    report.generate_report(levels, [], output,
                           limits={"10": (1.0, 120), "9": (2.0, 100), "2": (3.0, 80)})
    html = output.read_text()
    assert embedded_json(html, "data-level-order") == ["10", "9", "2"]
    assert embedded_json(html, "data-limits") == ["≤1% edits · ≥120 bp", "≤2% edits · ≥100 bp",
                                                  "≤3% edits · ≥80 bp"]
    assert [name for name, _ in embedded_json(html, "data-overview")["levels"]] == ["10", "9", "2"]
    for index in range(3):
        data = embedded_json(html, f"data-level-{index}")
        assert data["samples"] == list(samples)
        assert all(isinstance(sample, str) for sample in data["samples"])


def test_hostile_level_and_sample_names_never_appear_raw(tmp_path):
    report = load_report()
    suffix = '"</script><script>alert(1)</script>&'
    levels = [prefix + suffix for prefix in ('=', '+', '-', '@', '"')]
    samples = tuple(prefix + suffix for prefix in ('=', '+', '-', '@', '"'))
    table = frame(samples=samples)
    output = tmp_path / "hostile.html"
    report.generate_report({name: table for name in levels}, [], output)
    html = output.read_text()
    assert embedded_json(html, "data-level-order") == levels
    assert embedded_json(html, "data-level-0")["samples"] == list(samples)
    assert all(name not in html for name in [*levels, *samples])
    assert "\\u003c/script\\u003e" in html


def test_taxonomy_payload_is_dfs_and_preserves_cumulative_counts(tmp_path):
    report = load_report()
    hostile = 's__"</script><script>alert(3)</script>&'
    taxonomy = [
        "Unclassified", "k__Fungi", "k__Fungi;p__Ascomycota",
        "k__Fungi;p__Ascomycota;c__Eurotiomycetes",
        "k__Fungi;p__Ascomycota;c__Eurotiomycetes;o__Eurotiales",
        "k__Fungi;p__Ascomycota;c__Eurotiomycetes;o__Eurotiales;f__Trichocomaceae",
        "k__Fungi;p__Ascomycota;c__Eurotiomycetes;o__Eurotiales;f__Trichocomaceae;g__Aspergillus",
        "k__Fungi;p__Ascomycota;c__Eurotiomycetes;o__Eurotiales;f__Trichocomaceae;g__Aspergillus;" + hostile,
    ]
    table = pd.DataFrame({"Taxonomy": taxonomy, "1": [2, 10, 7, 7, 7, 4, 4, 4],
                          "2021-01": [0, 5, 5, 5, 5, 5, 5, 5]})
    output = tmp_path / "taxonomy.html"
    report.generate_report({"strict": table}, [], output)
    html = output.read_text()
    assert hostile not in html
    data = embedded_json(html, "data-level-0")
    taxa = data["taxa"]
    assert data["samples"] == ["1", "2021-01"]
    assert hostile in taxa["labels"]
    paths = []
    children = [[] for _ in taxa["labels"]]
    for index, (label, parent, depth) in enumerate(zip(taxa["labels"], taxa["parents"], taxa["depths"])):
        assert parent < index
        assert depth == (taxa["depths"][parent] + 1 if parent >= 0 else 0)
        paths.append((paths[parent] + ";" if parent >= 0 else "") + label)
        if parent >= 0:
            children[parent].append(index)
    for path, expected in zip(taxonomy, table["1"]):
        index = paths.index(path)
        assert dict(taxa["values"][index]).get(0, 0) == expected
    for index, descendants in enumerate(children):
        for sample in range(2):
            own = dict(taxa["values"][index]).get(sample, 0)
            summed = sum(dict(taxa["values"][child]).get(sample, 0) for child in descendants)
            assert own >= summed
            if paths[index] not in {"k__Fungi", taxonomy[4]}:
                assert own == summed or not descendants
    assert data["totals"] == [12, 5]


def _script_function(html, name):
    match = re.search(rf"    function {re.escape(name)}\([^\n]*?\) \{{.*?\n    \}}", html, re.DOTALL)
    assert match is not None
    return match.group(0)


def test_embedded_javascript_checks_and_pure_helpers(tmp_path):
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is unavailable")
    report = load_report()
    output = tmp_path / "scripts.html"
    report.generate_report({"strict": frame(), "standard": frame()}, [], output)
    html = output.read_text()
    blocks = re.findall(r'<script([^>]*)>(.*?)</script>', html, re.DOTALL)
    scripts = [body for attributes, body in blocks if 'type="application/json"' not in attributes]
    assert len(scripts) == 2
    for index, script in enumerate(scripts):
        path = tmp_path / f"embedded-{index}.js"
        path.write_text(script)
        subprocess.run([node, "--check", str(path)], check=True)

    helper_script = '''
let explicit = null;
const document = {documentElement:{getAttribute:()=>explicit,createElement:()=>({})},createElement:()=>({})};
const themePreference = {matches:true};
const colors = () => ({neutral:'#888',grid:'#ddd'});
const plotConfig = {};
const Plotly = {purge:()=>{},newPlot:(element)=>element.append({kind:'plot'})};
''' + "\n".join(_script_function(html, name) for name in
                  ("effectiveTheme", "taxonShade", "csvCell", "csvRows", "drawPlot", "showEmpty",
                   "taxonCount", "ensureTree", "taxonSampleOrder", "visibleTaxa", "expandedToDepth",
                   "changeTaxonExpansion", "taxonWindow", "windowCovers", "rovingTaxon", "taxonShare",
                   "taxonCsvData", "clearTaxonSearch", "setTaxonDepth", "drawBarplot", "drawDots")) + '''
if (effectiveTheme() !== 'dark') throw Error('OS dark preference ignored');
explicit='light';
if (effectiveTheme() !== 'light') throw Error('explicit light theme ignored');
for (const [mode,start,end] of [['light',30,75],['dark',40,85]]) {
    explicit=mode;
    const lightness=Array.from({length:10},(_,index)=>Number(taxonShade('k__Animalia',index,10).match(/([0-9.]+)%\\)$/)[1]));
    if (lightness[0] !== start || lightness[9] !== end) throw Error('shade bounds');
    if (lightness.some((value,index)=>index && value-lightness[index-1]<4.9)) throw Error('adjacent shades');
}
for (const prefix of ['=','+','-','@','"','\\t','\\r','\\n']) {
    if (!csvCell(prefix+'sample').startsWith('"\\'')) throw Error('unsafe CSV');
}
if (csvCell(-5)!=='-5' || csvCell(1e-7)!=='1e-7') throw Error('numeric cells changed');
if (csvCell('-5')!==`"'-5"`) throw Error('numeric-looking text was not guarded');
if (csvCell('=kingdom; species')!==`"'=kingdom; species"`) throw Error('lineage was not guarded');
if (csvCell('a,"b')!=='"a,""b"') throw Error('RFC quote escaping');
if (csvRows([['name','count'],['a',-5]])!=='\\ufeff"name","count"\\r\\n"a",-5\\r\\n') throw Error('CSV framing');
if (taxonCount([[0,3],[2,7]],2)!==7 || taxonCount([[0,3],[2,7]],1)!==0) throw Error('sparse taxon lookup');
const taxData={samples:['=sample','-5'],totals:[0,10],composition:[],
    taxa:{labels:['k__Fungi','p__P','c__C','s__S'],parents:[-1,0,1,2],depths:[0,1,2,3],
          values:[[[1,10]],[[1,8]],[[1,5]],[[1,2]]]}};
const tree=ensureTree(taxData);
if (visibleTaxa(taxData,tree,'').visible.join(',')!=='0,1') throw Error('default phylum depth');
const searched=visibleTaxa(taxData,tree,'S');
if (searched.visible.join(',')!=='0,1,2,3' || !searched.expanded.has(2)) throw Error('taxon search ancestry/expansion');
if (tree.expanded.has(1)) throw Error('search changed normal expansion');
tree.visible=searched.visible; tree.order=[1,0];
const ranks=['Kingdom','Phylum','Class','Species'];
const visibleCsv=taxonCsvData(taxData,tree,false,ranks);
if (visibleCsv.length!==5 || visibleCsv[0].join(',')!=='Taxon lineage,Rank,-5,=sample') throw Error('visible CSV/order');
if (visibleCsv[4][2]!==2 || visibleCsv[4][3]!==0) throw Error('CSV sample counts');
tree.visible=[0,1];
if (taxonCsvData(taxData,tree,false,ranks).length!==3 || taxonCsvData(taxData,tree,true,ranks).length!==5) throw Error('visible vs all CSV');
if (taxonSampleOrder(taxData,tree,'name').join(',')!=='1,0') throw Error('name sample sort');
if (taxonSampleOrder(taxData,tree,'total').join(',')!=='1,0') throw Error('total sample sort');
if (taxonShare(5,0)!==0 || taxonShare(5,10)!==.5) throw Error('share zero-total');
const indices=Array.from({length:3000},(_,index)=>index);
const firstWindow=taxonWindow(3000,100,0,540);
if (!firstWindow.virtual || (firstWindow.last-firstWindow.first)*100>6000) throw Error('cell budget');
if (!windowCovers(firstWindow,36,540,3000)) throw Error('programmatic scroll unnecessarily redraws');
if (windowCovers(firstWindow,1000,540,3000)) throw Error('distant scroll not redrawn');
const nextWindow=taxonWindow(3000,100,1000,540);
if (rovingTaxon(indices,nextWindow.first,nextWindow.last,15)!==nextWindow.first) throw Error('roving focus restore');
if (rovingTaxon(indices,firstWindow.first,firstWindow.last,15)!==15) throw Error('visible focus lost');
const searchInput={value:'S'}, depthInput={value:'1'}, grid={scrollTop:100};
document.getElementById=id=>({'taxon-search':searchInput,'taxon-depth':depthInput,'taxon-grid':grid})[id];
let redraws=0;
function drawTaxonGrid(){redraws++;}
let currentLevel=0;
function levelData(){return taxData;}
setTaxonDepth(0);
if (searchInput.value!=='' || tree.expanded.size || depthInput.value!=='0' || grid.scrollTop!==0 || redraws!==1) throw Error('search then collapse');
changeTaxonExpansion(taxData,tree,2,true);
if (visibleTaxa(taxData,tree,'').visible.join(',')!=='0,1,2,3') throw Error('ArrowRight ancestor expansion');
changeTaxonExpansion(taxData,tree,1,false);
if (visibleTaxa(taxData,tree,'').visible.join(',')!=='0,1') throw Error('ArrowLeft collapse');
const element={children:[],replaceChildren(){this.children=[]},append(child){this.children.push(child)}};
showEmpty(element,'empty');
if (element.children[0].textContent !== 'empty') throw Error('missing empty message');
drawPlot(element,[],{});
if (element.children.length !== 1 || element.children[0].kind !== 'plot') throw Error('stale empty message');
const barplot={style:{},children:[],replaceChildren(){this.children=[]},append(child){this.children.push(child)}};
const dotChart={style:{},children:[],replaceChildren(){this.children=[]},append(child){this.children.push(child)}};
const legend={children:['old legend'],replaceChildren(){this.children=[]}};
const compositionTable={children:['old row'],replaceChildren(){this.children=[]}};
document.querySelector=()=>({clientWidth:400});
document.getElementById=id=>({'barplot':barplot,'composition-legend':legend,'composition-table':compositionTable,'classified-chart':dotChart})[id];
function compositionOrder(){return [];}
drawBarplot();
if (compositionTable.children.length || legend.children.length) throw Error('stale empty-level table or legend');
if (barplot.children.length !== 1 || barplot.children[0].textContent !== 'No composition data for this level.') throw Error('empty-level message');
let dotSort=0, dotMeasure='count';
const levelNames=['strict'];
function filteredRows(){return [{name:'10',values:[10],fractions:[1]},{name:'2',values:[2],fractions:[.2]}];}
function levelColor(){return '#123456';}
function esc(value){return String(value);}
function layout(options){return options;}
let dotLayout;
Plotly.newPlot=(element,traces,figureLayout)=>{dotLayout=figureLayout;element.append({kind:'plot'});};
drawDots();
if (dotLayout.yaxis.type !== 'category') throw Error('numeric sample names inferred as numbers');
if (dotLayout.yaxis.categoryarray.join(',') !== '2,10') throw Error('sample sort lost');
'''
    helper_path = tmp_path / "helpers.js"
    helper_path.write_text(helper_script)
    subprocess.run([node, str(helper_path)], check=True)
