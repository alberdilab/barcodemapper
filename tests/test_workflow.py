"""Small end-to-end checks using real workflow executables."""

import gzip
import importlib.util
import json
import re
import runpy
import shutil
import subprocess
import sys
from itertools import groupby
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pysam
import pytest
import yaml

from barcodemapper import cli
from tests.make_fixture import TAXON_PREFIX, make_fixture


def test_extracted_rows_group_numeric_looking_read_ids(tmp_path):
    bam = tmp_path / "unsorted.bam"
    unsorted = tmp_path / "hits.tsv.unsorted"
    sorted_tsv = tmp_path / "hits.tsv"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "ref|ITS|k__Fungi", "LN": 100}]}
    with pysam.AlignmentFile(bam, "wb", header=header) as handle:
        for name in ("001", "1", "001"):
            read = pysam.AlignedSegment()
            read.query_name = name
            read.query_sequence = "A" * 100
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 0
            read.mapping_quality = 60
            read.cigar = [(0, 100)]
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            read.set_tag("NM", 0)
            handle.write(read)
    script = Path(__file__).resolve().parents[1] / "barcodemapper/workflow/scripts/extract_taxonomy.py"
    spec = importlib.util.spec_from_file_location("extract_taxonomy", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.process_bam(bam, unsorted)
    assert [line.split("\t", 1)[0] for line in unsorted.read_text().splitlines()] == ["001", "1", "001"]
    subprocess.run(
        ["bash", "-c", "set -euo pipefail; LC_ALL=C sort -s -t $'\\t' -k1,1 \"$1\" > \"$2\"", "bash", str(unsorted), str(sorted_tsv)],
        check=True,
    )
    names = [line.split("\t", 1)[0] for line in sorted_tsv.read_text().splitlines()]
    assert [(name, len(list(rows))) for name, rows in groupby(names)] == [("001", 2), ("1", 1)]


def test_stream_extraction_matches_bam_extraction(tmp_path):
    refs = ["ref1|ITS|k__Fungi;p__A", "ref2|COI|k__Animalia", "bare", "two|parts"]
    header = "@HD\tVN:1.6\n" + "".join(f"@SQ\tSN:{name}\tLN:1000\n" for name in refs)
    seq = "A" * 120
    qual = "I" * 120
    rows = [
        # name, flag, ref, cigar, tags
        ("keep", 0, "ref1|ITS|k__Fungi;p__A", "120M", "NM:i:0"),
        ("mixed", 0, "ref2|COI|k__Animalia", "10S100M5I5M", "NM:i:3"),
        ("eqx", 0, "ref2|COI|k__Animalia", "50=2X68=", "NM:i:2"),
        ("unmapped", 4, "*", "*", "YT:Z:UP"),
        ("nocigar", 0, "ref1|ITS|k__Fungi;p__A", "*", "NM:i:0"),
        ("short", 0, "ref1|ITS|k__Fungi;p__A", "70M50S", "NM:i:0"),
        ("noedit", 0, "ref1|ITS|k__Fungi;p__A", "120M", "AS:i:0"),
        ("toomany", 0, "ref1|ITS|k__Fungi;p__A", "120M", "NM:i:4"),
        ("edge", 0, "ref1|ITS|k__Fungi;p__A", "120M", "NM:i:3"),
        ("bare", 0, "bare", "120M", "NM:i:0"),
        ("two", 16, "two|parts", "120M", "NM:i:1"),
        ("secondary", 256, "ref2|COI|k__Animalia", "120M", "NM:i:0"),
        ("001", 0, "ref1|ITS|k__Fungi;p__A", "120M", "AS:i:0\tNM:i:0"),
        ("1", 0, "ref1|ITS|k__Fungi;p__A", "120M", "NM:i:0"),
        ("001", 0, "ref2|COI|k__Animalia", "120M", "NM:i:0"),
    ]
    lines = []
    for name, flag, ref, cigar, tags in rows:
        pos = 0 if ref == "*" else 1
        lines.append("\t".join([name, str(flag), ref, str(pos), "60" if flag != 4 else "0", cigar,
                               "*", "0", "0", seq, qual, tags]))
    sam = tmp_path / "in.sam"
    sam.write_text(header + "\n".join(lines) + "\n")
    bam = tmp_path / "in.bam"
    with pysam.AlignmentFile(sam, "r") as source, pysam.AlignmentFile(bam, "wb", template=source) as target:
        for read in source:
            target.write(read)
    scripts = Path(__file__).resolve().parents[1] / "barcodemapper/workflow/scripts"
    spec = importlib.util.spec_from_file_location("extract_taxonomy", scripts / "extract_taxonomy.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    for edit, covered in ((3.0, 80), (1.0, 50), (0.0, 30)):
        expected = tmp_path / "bam.tsv"
        module.process_bam(bam, expected, max_edit_percent=edit, min_covered=covered)
        streamed = tmp_path / "stream.tsv"
        with open(sam) as handle:
            subprocess.run(
                [sys.executable, str(scripts / "stream_extract_taxonomy.py"), "--output", str(streamed),
                 "--max-edit-percent", str(edit), "--min-covered", str(covered)],
                stdin=handle, check=True,
            )
        assert streamed.read_text().splitlines() == expected.read_text().splitlines()
        assert expected.read_text()  # not vacuous


def test_stream_multilevel_stats_count_contiguous_fragments(tmp_path):
    script = cli.PACKAGE_DIR / "workflow/scripts/stream_extract_taxonomy.py"

    def record(name, cigar, edits, flag=0):
        return "\t".join((name, str(flag), "ref|ITS|k__Fungi", "1", "60", cigar,
                          "*", "0", "0", "A" * 120, "I" * 120, edits)) + "\n"

    sam = ("@HD\tVN:1.6\n"
           + record("r1", "120M", "NM:i:0")
           + record("r1", "120M", "NM:i:2")
           + record("r1", "90M30S", "NM:i:2")
           + record("r2", "70M50S", "NM:i:3")
           + record("r3", "*", "NM:i:0")
           + record("r3", "120M", "AS:i:0")
           + record("r4", "*", "YT:Z:UP", flag=4))
    levels = {"ultrastrict": (0, 120), "strict": (2, 100),
              "standard": (3, 80), "permissive": (5, 60)}
    command = [sys.executable, str(script)]
    for name, (percent, coverage) in levels.items():
        command += ["--level", f"{name}:{percent}:{coverage}"]
    command += ["--output-pattern", str(tmp_path / "{level}.tsv"),
                "--stats-json", str(tmp_path / "stats.json"), "--sample", "sample"]
    subprocess.run(command, input=sam, text=True, check=True)
    assert json.loads((tmp_path / "stats.json").read_text()) == {
        "sample": "sample", "records": 7, "aligned_fragments": 3,
        "levels": {
            "ultrastrict": {"alignments": 1, "fragments": 1},
            "strict": {"alignments": 2, "fragments": 1},
            "standard": {"alignments": 3, "fragments": 1},
            "permissive": {"alignments": 4, "fragments": 2},
        },
    }


@pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools not installed")
def test_stream_extraction_can_also_write_bam(tmp_path):
    sam = tmp_path / "in.sam"
    sam.write_text("@HD\tVN:1.6\n@SQ\tSN:ref1|ITS|k__Fungi\tLN:1000\n"
                   + "\t".join(["r1", "0", "ref1|ITS|k__Fungi", "1", "60", "120M", "*", "0", "0",
                               "A" * 120, "I" * 120, "NM:i:0"]) + "\n")
    scripts = Path(__file__).resolve().parents[1] / "barcodemapper/workflow/scripts"
    bam = tmp_path / "out.bam"
    with open(sam) as handle:
        subprocess.run([sys.executable, str(scripts / "stream_extract_taxonomy.py"), "--output",
                        str(tmp_path / "out.tsv"), "--bam-out", str(bam)], stdin=handle, check=True)
    with pysam.AlignmentFile(bam, "rb") as reads:
        assert [read.query_name for read in reads] == ["r1"]
    assert (tmp_path / "out.tsv").read_text().startswith("r1\t")


def test_multilevel_stages_match_legacy_standard_and_render_report(tmp_path):
    scripts = cli.PACKAGE_DIR / "workflow/scripts"
    levels = list(cli.PRESETS)
    samples = ("positive", "zero")
    references = {
        "fungus_one": "f1|ITS|k__Fungi;p__Ascomycota;g__Alpha;s__one",
        "fungus_two": "f2|ITS|k__Fungi;p__Ascomycota;g__Alpha;s__two",
        "animal": "a1|COI|k__Animalia;p__Arthropoda",
        "bare": "bare",
    }

    def record(name, reference="fungus_one", cigar="120M", edits="NM:i:0", flag=0):
        target = references.get(reference, "*")
        return "\t".join((name, str(flag), target, "1" if flag != 4 else "0", "60",
                          cigar, "*", "0", "0", "A" * 120, "I" * 120, edits)) + "\n"

    sam_by_sample = {
        "positive": (
            "@HD\tVN:1.6\n"
            + record("perfect")
            + record("shared", "fungus_one")
            + record("shared", "fungus_two")
            + record("strict", "animal", "110M10S", "NM:i:2")
            + record("standard", "fungus_one", "90M30S", "NM:i:2")
            + record("permissive", "animal", "70M50S", "NM:i:3")
            + record("unclassified", "bare")
            + record("missing_nm", edits="AS:i:0")
            + record("no_cigar", cigar="*")
            + record("unmapped", reference="*", cigar="*", flag=4)
        ),
        "zero": "@HD\tVN:1.6\n" + record("unmapped", reference="*", cigar="*", flag=4),
    }

    def run_taxonomy_stages(extracts, destination):
        aggregates = []
        for sample in samples:
            extracted = extracts[sample]
            lines = extracted.read_text().splitlines(keepends=True)
            extracted.write_text("".join(sorted(lines, key=lambda line: line.split("\t", 1)[0])))
            assigned = destination / "taxonomy" / f"{sample}.tsv"
            aggregate = destination / "aggregate" / f"{sample}.tsv"
            assigned.parent.mkdir(parents=True, exist_ok=True)
            aggregate.parent.mkdir(parents=True, exist_ok=True)
            runpy.run_path(str(scripts / "assign_taxonomy.py"), init_globals={
                "snakemake": SimpleNamespace(input=[str(extracted)], output=[str(assigned)])
            })
            runpy.run_path(str(scripts / "aggregate_taxonomy.py"), init_globals={
                "snakemake": SimpleNamespace(input=[str(assigned)], output=[str(aggregate)])
            })
            aggregates.append(str(aggregate))
        merged = destination / "merged.tsv"
        runpy.run_path(str(scripts / "merge_taxonomy.py"), init_globals={
            "snakemake": SimpleNamespace(input=aggregates, output=[str(merged)],
                                         params=SimpleNamespace(samples=list(samples)))
        })
        return merged

    for sample in samples:
        command = [sys.executable, str(scripts / "stream_extract_taxonomy.py")]
        for name, (percent, coverage) in cli.PRESETS.items():
            command += ["--level", f"{name}:{percent}:{coverage}"]
        command += ["--output-pattern", str(tmp_path / "extract/{level}" / f"{sample}.tsv"),
                    "--stats-json", str(tmp_path / "extract" / f"{sample}.stats.json"),
                    "--sample", sample]
        subprocess.run(command, input=sam_by_sample[sample], text=True, check=True)

    merged_by_level = {}
    for level in levels:
        extracts = {sample: tmp_path / "extract" / level / f"{sample}.tsv" for sample in samples}
        assert all(path.is_file() for path in extracts.values())
        merged_by_level[level] = run_taxonomy_stages(extracts, tmp_path / level)
    assert all(path.is_file() for path in merged_by_level.values())

    legacy_extracts = {}
    for sample in samples:
        extracted = tmp_path / "legacy" / f"{sample}.tsv"
        extracted.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run([
            sys.executable, str(scripts / "stream_extract_taxonomy.py"),
            "--output", str(extracted), "--max-edit-percent", "3", "--min-covered", "80",
        ], input=sam_by_sample[sample], text=True, check=True)
        legacy_extracts[sample] = extracted
    legacy_merged = run_taxonomy_stages(legacy_extracts, tmp_path / "legacy")
    assert merged_by_level["standard"].read_bytes() == legacy_merged.read_bytes()

    stats_paths = [tmp_path / "extract" / f"{sample}.stats.json" for sample in samples]
    stats = [json.loads(path.read_text()) for path in stats_paths]
    assert stats[0]["records"] == 10
    assert stats[0]["aligned_fragments"] == 7
    assert [stats[0]["levels"][level]["fragments"] for level in reversed(levels)] == [6, 5, 4, 3]
    assert [stats[0]["levels"][level]["alignments"] for level in reversed(levels)] == [7, 6, 5, 4]
    assert all(stats[1]["levels"][level]["fragments"] == 0 for level in levels)
    for merged in merged_by_level.values():
        table = pd.read_csv(merged, sep="\t")
        assert list(table.columns) == ["Taxonomy", *samples]
        assert table["zero"].eq(0).all()
        assert table.loc[table.Taxonomy == "Unclassified", "positive"].iloc[0] == 1

    report = tmp_path / "report.html"
    report_command = [sys.executable, str(scripts / "html_report.py")]
    for level in levels:
        percent, coverage = cli.PRESETS[level]
        report_command += ["--input", f"{level}={merged_by_level[level]}",
                           "--limits", f"{level}={percent}:{coverage}"]
    for path in stats_paths:
        report_command += ["--stats", str(path)]
    subprocess.run(report_command + ["--output", str(report), "--version", "test"], check=True)
    html = report.read_text()
    assert 'id="tab-overview"' in html
    assert 'id="overview-panel"' in html
    assert [match.group(1) for match in re.finditer(
        r'id="tab-level-\d+"[^>]*><strong>([^<]+)</strong>', html)] == levels


def _require_tools(*names: str) -> None:
    missing = [name for name in names if shutil.which(name) is None]
    if missing:
        pytest.skip("workflow executables unavailable: " + ", ".join(missing))


def _use_local_storage(monkeypatch, root: Path) -> None:
    original = cli._snakemake_command

    def command(tmp_dir, config_path, profile):
        return original(tmp_dir, config_path, profile) + [
            "--local-storage-prefix", str(root / "snakemake-storage"),
            "--shared-fs-usage", "input-output", "persistence", "software-deployment",
            "sources", "storage-local-copies",
        ]

    monkeypatch.setattr(cli, "_snakemake_command", command)


@pytest.mark.integration
@pytest.mark.parametrize(("max_hits", "expected", "excluded", "keep_bam"), [
    (None, "-a", "-k 50", True),
    (50, "-k 50", "-a", True),
    (50, "-k 50", "-a", False),
])
def test_bowtie2_mapping_renders_hit_limit_without_bowtie2(tmp_path, max_hits, expected, excluded, keep_bam):
    _require_tools("snakemake")
    database, reads = make_fixture(tmp_path)
    work = tmp_path / "work"
    (work / "data").mkdir(parents=True)
    (work / "fastp").mkdir()
    index_base = database.with_suffix("")
    for part in cli.INDEX_PARTS:
        Path(f"{index_base}.{part}.bt2l").touch()
    for mate in (1, 2):
        (work / "fastp" / f"positive_{mate}.fq.gz").touch()
        (work / "data" / f"sample_to_reads{mate}.json").write_text(
            json.dumps({"positive": [str(reads / f"positive_{mate}.fq.gz")]})
        )
    config_path = work / "data/config.yaml"
    config_path.write_text(yaml.safe_dump({
        "package_dir": str(cli.PACKAGE_DIR),
        "barcodemapper_db": str(database),
        "barcodemapper_index_base": str(index_base),
        "generate_index": False,
        "generate_headers": False,
        "barcodemapper_db_headers": str(tmp_path / "headers.tsv"),
        "tmp_dir": str(work),
        "output_file": str(tmp_path / "results.tsv"),
        "report_file": str(tmp_path / "results.html"),
        "build_only": "no",
        "max_hits": max_hits,
        "keep_bam": keep_bam,
    }))
    bam = work / ("bowtie/positive.bam" if keep_bam else "extract/custom/positive.tsv")
    result = subprocess.run([
        "snakemake", "-s", str(cli.PACKAGE_DIR / "workflow/Snakefile"),
        "--directory", str(work), "--configfile", str(config_path),
        "--local-storage-prefix", str(tmp_path / "snakemake-storage"),
        "--shared-fs-usage", "input-output", "persistence", "software-deployment",
        "sources", "storage-local-copies",
        "--cores", "1", "-n", "-p", str(bam),
    ], capture_output=True, text=True)
    command = result.stdout + result.stderr
    assert result.returncode == 0, command
    assert "bowtie2" in command and "--seed 1337" in command and "--no-unal" in command
    assert ("--bam-out" in command) is keep_bam
    assert re.search(r"(?<!\S)" + re.escape(expected) + r"(?!\S)", command)
    assert not re.search(r"(?<!\S)" + re.escape(excluded) + r"(?!\S)", command)


@pytest.mark.parametrize("level_names", [list(cli.PRESETS), ["standard"]])
def test_multilevel_dry_run_maps_each_sample_once(tmp_path, level_names):
    pytest.importorskip("snakemake")
    work = tmp_path / "work"
    (work / "data").mkdir(parents=True)
    (work / "fastp").mkdir()
    database = tmp_path / "reference.fa"
    database.write_text(">ref|ITS|k__Fungi\nACGT\n")
    index_base = database.with_suffix("")
    for part in cli.INDEX_PARTS:
        Path(f"{index_base}.{part}.bt2l").touch()
    headers = tmp_path / "headers.tsv"
    headers.touch()
    samples = ["alpha", "beta.v2"]
    for mate in (1, 2):
        reads = {}
        for sample in samples:
            path = tmp_path / f"{sample}_{mate}.fq.gz"
            path.touch()
            reads[sample] = [str(path)]
            (work / "fastp" / f"{sample}_{mate}.fq.gz").touch()
        (work / "data" / f"sample_to_reads{mate}.json").write_text(json.dumps(reads))
    levels = {name: {"max_edit_percent": cli.PRESETS[name][0],
                     "min_coverage": cli.PRESETS[name][1]} for name in level_names}
    outputs = {name: str(tmp_path / f"result.{name}.tsv") for name in level_names}
    config_path = work / "data/config.yaml"
    config_path.write_text(yaml.safe_dump({
        "package_dir": str(cli.PACKAGE_DIR),
        "version": "test",
        "barcodemapper_db": str(database),
        "barcodemapper_index_base": str(index_base),
        "generate_index": False,
        "generate_headers": False,
        "barcodemapper_db_headers": str(headers),
        "tmp_dir": str(work),
        "output_file": str(tmp_path / "result.tsv"),
        "output_files": outputs,
        "report_file": str(tmp_path / "result.html"),
        "build_only": "no",
        "levels": levels,
        "keep_bam": False,
    }, sort_keys=False))
    result = subprocess.run([
        "snakemake", "-s", str(cli.PACKAGE_DIR / "workflow/Snakefile"),
        "--directory", str(work), "--configfile", str(config_path),
        "--local-storage-prefix", str(tmp_path / "snakemake-storage"),
        "--shared-fs-usage", "input-output", "persistence", "software-deployment",
        "sources", "storage-local-copies",
        "--cores", "1", "-n", "-p",
    ], capture_output=True, text=True)
    command = result.stdout + result.stderr
    assert result.returncode == 0, command
    job_stats = command.split("Job stats:", 2)[1]
    for rule, expected in (("bowtie2_mapping", len(samples)),
                           ("assign_taxonomy", len(samples) * len(levels)),
                           ("aggregate_taxonomy", len(samples) * len(levels)),
                           ("merge_taxonomy", len(levels)),
                           ("publish_results", 1), ("html_report", 1)):
        assert re.search(rf"^{rule}\s+{expected}$", job_stats, re.MULTILINE), command
    for name in level_names:
        assert f"--level {name}:{levels[name]['max_edit_percent']}:{levels[name]['min_coverage']}" in command
    report_command = next(line for line in command.splitlines() if "html_report.py" in line)
    input_flags = [f"--input {name}={outputs[name]}" for name in level_names]
    limit_flags = [f"--limits {name}={levels[name]['max_edit_percent']}:{levels[name]['min_coverage']}"
                   for name in level_names]
    assert all(flag in report_command for flag in input_flags)
    assert [report_command.index(flag) for flag in input_flags] == sorted(report_command.index(flag) for flag in input_flags)
    assert all(flag in report_command for flag in limit_flags)
    assert all(f"--stats {work}/extract/{sample}.stats.json" in report_command for sample in samples)


@pytest.mark.integration
def test_synthetic_analysis_keeps_ambiguous_and_zero_hit_samples(tmp_path, monkeypatch):
    _require_tools("snakemake", "fastp", "bowtie2", "bowtie2-build", "samtools")
    _use_local_storage(monkeypatch, tmp_path)
    database, reads = make_fixture(tmp_path)
    output = tmp_path / "results.tsv"
    cli.main([
        "-i", str(reads), "-d", str(database), "-o", str(output),
        "-t", str(tmp_path / "work"), "--adapter-r1", "auto", "--adapter-r2", "auto",
        "--standard",
    ])

    table = pd.read_csv(output, sep="\t")
    assert list(table.columns) == ["Taxonomy", "negative", "positive"]
    assert table.loc[table.Taxonomy == TAXON_PREFIX, "positive"].iloc[0] == 1
    assert table["negative"].eq(0).all()
    assert not table.Taxonomy.str.contains("s__Synthus_").any()
    assert output.with_suffix(".html").is_file()
    manifest = json.loads(output.with_suffix(".manifest.json").read_text())
    assert manifest["samples"]["positive"]["assigned_read_ids"] == {"standard": 1}
    assert manifest["samples"]["negative"]["assigned_read_ids"] == {"standard": 0}
    assert len(manifest["database"]["sha256"]) == 64
    assert (tmp_path / "work/reference/synthetic.1.bt2l").is_file()
    assert (tmp_path / "work/reference/synthetic.tsv").is_file()
    assert not database.with_suffix(".1.bt2l").exists()
    assert not database.with_suffix(".tsv").exists()


@pytest.mark.integration
def test_build_records_source_checksums_and_conflicting_taxa(tmp_path, monkeypatch):
    _require_tools("snakemake", "bowtie2-build")
    _use_local_storage(monkeypatch, tmp_path)
    bold = tmp_path / "bold.fa"
    bold.write_text(
        ">A|COI-5P|x|Animalia,Chordata,Mammalia,Primates,Hominidae,Synthus,one\n"
        + "ACGT" * 40 + "\n"
        ">B|COI-5P|x|Animalia,Chordata,Mammalia,Primates,Hominidae,Synthus,two\n"
        + "ACGT" * 40 + "\n",
        encoding="utf-8",
    )
    unite = tmp_path / "unite.fa"
    unite.write_text(
        ">record|U1|x|x|k__Fungi;p__Basidiomycota;c__Agaricomycetes\n"
        + "TGCA" * 40 + "\n",
        encoding="utf-8",
    )
    database = tmp_path / "combined.fa.gz"
    cli.main([
        "--build", "-b", str(bold), "-u", str(unite), "-d", str(database),
        "-t", str(tmp_path / "work"), "-y", "k__Fungi",
    ])

    with gzip.open(database, "rt") as handle:
        fasta = handle.read()
    assert fasta.count(">") == 3
    assert "s__one" in fasta and "s__two" in fasta
    assert "U1|ITS|k__Fungi" in fasta
    assert cli._database_base(database).with_suffix(".1.bt2l").is_file()
    manifest = json.loads(cli._database_base(database).with_suffix(".manifest.json").read_text())
    assert len(manifest["sources"]["bold"]["sha256"]) == 64
    assert len(manifest["sources"]["unite"]["sha256"]) == 64
