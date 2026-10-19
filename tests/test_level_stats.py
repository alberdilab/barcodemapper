"""Checks for the per-level merge and report overview data."""

import importlib.util
import json
import runpy
from collections import defaultdict
from pathlib import Path
from types import SimpleNamespace

import pandas as pd


SCRIPTS = Path(__file__).resolve().parents[1] / "barcodemapper" / "workflow" / "scripts"


def load_level_stats():
    spec = importlib.util.spec_from_file_location("level_stats", SCRIPTS / "level_stats.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def reference_merge(input_files, sample_names, output_file):
    """The original merge_taxonomy.py implementation, kept for byte comparison."""
    taxonomy_dict = defaultdict(lambda: defaultdict(int))

    for file, sample_name in zip(input_files, sample_names):
        df = pd.read_csv(file, sep="\t")
        for _, row in df.iterrows():
            taxonomy = row["Taxonomy"]
            count = row["Count"]
            taxonomy_dict[taxonomy][sample_name] = count

    merged_df = pd.DataFrame.from_dict(taxonomy_dict, orient="index").fillna(0).astype(int)
    for sample_name in sample_names:
        if sample_name not in merged_df:
            merged_df[sample_name] = 0
    merged_df = merged_df[sample_names]
    merged_df.index.name = "Taxonomy"
    merged_df = merged_df.reset_index()
    merged_df = merged_df.sort_values(by="Taxonomy")
    merged_df.to_csv(output_file, sep="\t", index=False)


def test_merge_is_byte_identical_to_original(tmp_path):
    one = tmp_path / "one.tsv"
    two = tmp_path / "two.tsv"
    empty = tmp_path / "empty.tsv"
    one.write_text("Taxonomy\tCount\nk__Fungi;p__Ascomycota\t4\nUnclassified\t2\nk__Fungi\t8\nk__Fungi\t9\n")
    two.write_text("Taxonomy\tCount\nk__Animalia\t3\nk__Fungi\t5\nUnclassified\t0\n")
    empty.write_text("Taxonomy\tCount\n")
    files = [one, two, empty]
    samples = ["one", "two", "empty"]
    expected = tmp_path / "expected.tsv"
    actual = tmp_path / "actual.tsv"
    reference_merge(files, samples, expected)
    mock = SimpleNamespace(input=list(map(str, files)), output=[str(actual)], params=SimpleNamespace(samples=samples))
    runpy.run_path(str(SCRIPTS / "merge_taxonomy.py"), init_globals={"snakemake": mock})
    assert actual.read_bytes() == expected.read_bytes()


def test_merge_all_empty_is_byte_identical_to_original(tmp_path):
    empty = tmp_path / "empty.tsv"
    empty.write_text("Taxonomy\tCount\n")
    expected = tmp_path / "expected.tsv"
    actual = tmp_path / "actual.tsv"
    reference_merge([empty], ["zero"], expected)
    mock = SimpleNamespace(input=[str(empty)], output=[str(actual)], params=SimpleNamespace(samples=["zero"]))
    runpy.run_path(str(SCRIPTS / "merge_taxonomy.py"), init_globals={"snakemake": mock})
    assert actual.read_bytes() == expected.read_bytes() == b"Taxonomy\tzero\n"


def test_overview_stats_uses_cumulative_counts_once_and_preserves_level_order():
    stats = load_level_stats()
    deep = "k__Fungi;p__Ascomycota;c__Saccharomycetes"
    strict = pd.DataFrame({
        "Taxonomy": ["Unclassified", "k__Animalia", "k__Fungi", "k__Fungi;p__Ascomycota", deep],
        "alpha": [2, 1, 8, 6, 3],
        "beta": [0, 0, 4, 2, 0],
    })
    permissive = pd.DataFrame({"Taxonomy": ["Unclassified", "k__Fungi"], "alpha": [3, 10], "beta": [1, 5]})
    extract_stats = [
        {"sample": "alpha", "aligned_fragments": 20, "levels": {
            "strict": {"fragments": 12}, "permissive": {"fragments": 15}
        }},
    ]
    result = stats.overview_stats({"strict": strict, "permissive": permissive}, extract_stats)
    assert list(result["levels"]) == ["strict", "permissive"]
    alpha = result["levels"]["strict"]["samples"]["alpha"]
    assert alpha == {
        "classified_fragments": 9,
        "unclassified_fragments": 2,
        "depth_distribution": [9, 6, 3, 0, 0, 0, 0],
        "distinct_taxa": [2, 1, 1, 0, 0, 0, 0],
        "kingdom_composition": {"k__Animalia": 1, "k__Fungi": 8},
        "fragments_retained": 12,
        "aligned_fragments": 20,
        "retained_fraction": 0.6,
    }
    beta = result["levels"]["strict"]["samples"]["beta"]
    assert beta["classified_fragments"] == 4
    assert beta["depth_distribution"] == [4, 2, 0, 0, 0, 0, 0]
    assert beta["fragments_retained"] is None
    assert beta["aligned_fragments"] is None
    assert beta["retained_fraction"] is None
    totals = result["levels"]["strict"]["totals"]
    assert totals["classified_fragments"] == 13
    assert totals["unclassified_fragments"] == 2
    assert totals["depth_distribution"] == [13, 8, 3, 0, 0, 0, 0]
    assert totals["distinct_taxa"] == [2, 1, 1, 0, 0, 0, 0]
    assert totals["kingdom_composition"] == {"k__Animalia": 1, "k__Fungi": 12}
    assert totals["fragments_retained"] == 12
    assert totals["aligned_fragments"] == 20
    assert totals["retained_fraction"] == 0.6
    assert result["classified_by_sample"] == [
        {"level": "strict", "sample": "alpha", "classified_fragments": 9},
        {"level": "strict", "sample": "beta", "classified_fragments": 4},
        {"level": "permissive", "sample": "alpha", "classified_fragments": 10},
        {"level": "permissive", "sample": "beta", "classified_fragments": 5},
    ]
    json.dumps(result)


def test_overview_stats_handles_empty_levels_and_missing_extract_stats(tmp_path):
    stats = load_level_stats()
    empty = tmp_path / "empty.tsv"
    empty.write_text("Taxonomy\talpha\tbeta\n")
    unclassified = tmp_path / "unclassified.tsv"
    unclassified.write_text("Taxonomy\talpha\tbeta\nUnclassified\t4\t0\n")
    all_zero = tmp_path / "all_zero.tsv"
    all_zero.write_text("Taxonomy\talpha\tbeta\nk__Fungi\t0\t0\nUnclassified\t0\t0\n")
    stats_path = tmp_path / "alpha.stats.json"
    stats_path.write_text(json.dumps({"sample": "alpha", "aligned_fragments": 0, "levels": {
        "empty": {"fragments": 0}, "unclassified": {"fragments": 4}, "all_zero": {"fragments": 0}
    }}))
    result = stats.overview_stats({
        "empty": stats.load_merged(empty),
        "unclassified": stats.load_merged(unclassified),
        "all_zero": stats.load_merged(all_zero),
    }, stats.load_stats([stats_path]))
    for name in ("empty", "all_zero"):
        sample = result["levels"][name]["samples"]["alpha"]
        assert sample["classified_fragments"] == 0
        assert sample["unclassified_fragments"] == 0
        assert sample["depth_distribution"] == [0] * 7
        assert sample["distinct_taxa"] == [0] * 7
        assert sample["kingdom_composition"] == {}
        assert sample["retained_fraction"] is None
    assert result["levels"]["unclassified"]["samples"]["alpha"]["unclassified_fragments"] == 4
    assert result["levels"]["unclassified"]["totals"]["unclassified_fragments"] == 4
    assert result["levels"]["empty"]["samples"]["beta"]["fragments_retained"] is None
    assert stats.overview_stats({}, []) == {"levels": {}, "classified_by_sample": []}
    json.dumps(result)


def test_skipped_rank_still_counts_as_resolved_to_deeper_rank():
    stats = load_level_stats()
    frame = pd.DataFrame({
        "Taxonomy": ["k__Fungi", "k__Fungi;g__Alpha"],
        "sample": [3, 2],
    })
    summary = stats.overview_stats({"strict": frame}, [])["levels"]["strict"]["samples"]["sample"]
    assert summary["depth_distribution"] == [3, 2, 2, 2, 2, 2, 0]
    assert summary["distinct_taxa"] == [1, 0, 0, 0, 0, 1, 0]
