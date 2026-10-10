import csv
import importlib.util
import runpy
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pysam


SCRIPTS = Path(__file__).resolve().parents[1] / "barcodemapper" / "workflow" / "scripts"


def test_unknown_intermediate_rank_keeps_resolved_ancestor(tmp_path):
    assigned = tmp_path / "assigned.tsv"
    assigned.write_text("ReadID\tTaxonomy\nr1\tk__Fungi;p__;c__Agaricomycetes\nr2\t\n")
    output = tmp_path / "aggregate.tsv"
    runpy.run_path(str(SCRIPTS / "aggregate_taxonomy.py"), init_globals={"snakemake": SimpleNamespace(input=[str(assigned)], output=[str(output)])})
    counts = dict(csv.reader(output.open(), delimiter="\t"))
    assert counts["k__Fungi"] == "1"
    assert counts["Unclassified"] == "1"


def test_zero_hit_sample_remains_in_table_and_report(tmp_path):
    positive = tmp_path / "positive.tsv"
    negative = tmp_path / "negative.tsv"
    positive.write_text("Taxonomy\tCount\nk__Fungi\t2\n")
    negative.write_text("Taxonomy\tCount\n")
    merged = tmp_path / "merged.tsv"
    mock = SimpleNamespace(input=[str(positive), str(negative)], output=[str(merged)], params=SimpleNamespace(samples=["positive", "negative"]))
    runpy.run_path(str(SCRIPTS / "merge_taxonomy.py"), init_globals={"snakemake": mock})
    table = pd.read_csv(merged, sep="\t")
    assert list(table.columns) == ["Taxonomy", "positive", "negative"]
    assert table["negative"].tolist() == [0]

    spec = importlib.util.spec_from_file_location("html_report", SCRIPTS / "html_report.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    samples, data, absolute, relative = module.prepare_data(table)
    report = tmp_path / "report.html"
    module.generate_html(samples, data, report, absolute, relative, "1.0.0")
    html = report.read_text()
    assert "plotly-latest" not in html
    assert "No data for this sample" in html


def test_edit_percentage_uses_actual_aligned_length(tmp_path):
    spec = importlib.util.spec_from_file_location("extract_taxonomy", SCRIPTS / "extract_taxonomy.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    bam = tmp_path / "reads.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "ref|ITS|k__Fungi", "LN": 200}]}
    with pysam.AlignmentFile(bam, "wb", header=header) as handle:
        for name, edits in (("accepted", 2), ("rejected", 3)):
            read = pysam.AlignedSegment()
            read.query_name = name
            read.query_sequence = "A" * 100
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 0
            read.mapping_quality = 60
            read.cigar = [(0, 100)]
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            read.set_tag("NM", edits)
            handle.write(read)
    output = tmp_path / "hits.tsv"
    module.process_bam(bam, output, max_edit_percent=2.0, min_covered=80)
    text = output.read_text()
    assert "accepted\t" in text
    assert "rejected\t" not in text


def test_empty_alignment_file_produces_empty_taxonomy_with_header(tmp_path):
    extracted = tmp_path / "empty.tsv"
    extracted.touch()
    assigned = tmp_path / "assigned.tsv"
    runpy.run_path(str(SCRIPTS / "assign_taxonomy.py"), init_globals={
        "snakemake": SimpleNamespace(input=[str(extracted)], output=[str(assigned)])
    })
    assert assigned.read_text() == "ReadID\tTaxonomy\n"

    aggregated = tmp_path / "aggregate.tsv"
    runpy.run_path(str(SCRIPTS / "aggregate_taxonomy.py"), init_globals={
        "snakemake": SimpleNamespace(input=[str(assigned)], output=[str(aggregated)])
    })
    assert aggregated.read_text() == "Taxonomy\tCount\n"


def test_report_barplot_includes_kingdom_only_assignments():
    spec = importlib.util.spec_from_file_location("html_report", SCRIPTS / "html_report.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    table = pd.DataFrame({
        "Taxonomy": ["k__Animalia", "k__Animalia;p__Chordata", "Unclassified"],
        "sample": [5, 3, 2],
    })
    _, _, absolute, relative = module.prepare_data(table)
    assert absolute.set_index("Phylum").loc["Unresolved within k__Animalia", "sample"] == 2
    assert absolute["sample"].sum() == 7
    assert relative["sample"].sum() == 1


def test_report_escapes_sample_names(tmp_path):
    spec = importlib.util.spec_from_file_location("html_report", SCRIPTS / "html_report.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    hostile = 'x\"><script>alert(1)</script>'
    table = pd.DataFrame({"Taxonomy": ["k__Fungi"], hostile: [1]})
    samples, data, absolute, relative = module.prepare_data(table)
    output = tmp_path / "report.html"
    module.generate_html(samples, data, output, absolute, relative, "1.0.0")
    html = output.read_text()
    assert hostile not in html
    assert "&lt;script&gt;" in html
    assert "\\u003cscript\\u003e" in html


def test_report_handles_all_zero_samples(tmp_path):
    spec = importlib.util.spec_from_file_location("html_report", SCRIPTS / "html_report.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    table = pd.DataFrame(columns=["Taxonomy", "empty_one", "empty_two"])
    samples, data, absolute, relative = module.prepare_data(table)
    assert samples == ["empty_one", "empty_two"]
    assert absolute.empty and relative.empty
    output = tmp_path / "empty.html"
    module.generate_html(samples, data, output, absolute, relative, "1.0.0")
    assert "No data for this sample" in output.read_text()
