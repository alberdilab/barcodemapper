"""Small end-to-end checks using real workflow executables."""

import gzip
import importlib.util
import json
import re
import shutil
import subprocess
from itertools import groupby
from pathlib import Path

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
@pytest.mark.parametrize(("max_hits", "expected", "excluded"), [
    (None, "-a", "-k 50"),
    (50, "-k 50", "-a"),
])
def test_bowtie2_mapping_renders_hit_limit_without_bowtie2(tmp_path, max_hits, expected, excluded):
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
    }))
    bam = work / "bowtie/positive.bam"
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
    assert "bowtie2" in command and "--seed 1337" in command
    assert re.search(r"(?<!\S)" + re.escape(expected) + r"(?!\S)", command)
    assert not re.search(r"(?<!\S)" + re.escape(excluded) + r"(?!\S)", command)


@pytest.mark.integration
def test_synthetic_analysis_keeps_ambiguous_and_zero_hit_samples(tmp_path, monkeypatch):
    _require_tools("snakemake", "fastp", "bowtie2", "bowtie2-build", "samtools")
    _use_local_storage(monkeypatch, tmp_path)
    database, reads = make_fixture(tmp_path)
    output = tmp_path / "results.tsv"
    cli.main([
        "-i", str(reads), "-d", str(database), "-o", str(output),
        "-t", str(tmp_path / "work"), "--adapter-r1", "auto", "--adapter-r2", "auto",
    ])

    table = pd.read_csv(output, sep="\t")
    assert list(table.columns) == ["Taxonomy", "negative", "positive"]
    assert table.loc[table.Taxonomy == TAXON_PREFIX, "positive"].iloc[0] == 1
    assert table["negative"].eq(0).all()
    assert not table.Taxonomy.str.contains("s__Synthus_").any()
    assert output.with_suffix(".html").is_file()
    manifest = json.loads(output.with_suffix(".manifest.json").read_text())
    assert manifest["samples"]["positive"]["assigned_read_ids"] == 1
    assert manifest["samples"]["negative"]["assigned_read_ids"] == 0
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
