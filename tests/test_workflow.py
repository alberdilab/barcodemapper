"""Small end-to-end checks using real workflow executables."""

import json
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pysam
import pytest

from barcodemapper import cli
from tests.make_fixture import TAXON_PREFIX, make_fixture


@pytest.mark.integration
def test_collation_keeps_numeric_looking_read_ids_together(tmp_path):
    _require_tools("samtools")
    unsorted = tmp_path / "unsorted.bam"
    collated = tmp_path / "collated.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "ref|ITS|k__Fungi", "LN": 100}]}
    with pysam.AlignmentFile(unsorted, "wb", header=header) as handle:
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
            handle.write(read)
    subprocess.run(["samtools", "collate", "-o", str(collated), "-"], input=unsorted.read_bytes(), check=True)
    with pysam.AlignmentFile(collated, "rb") as handle:
        names = [read.query_name for read in handle]
    assert names in (["001", "001", "1"], ["1", "001", "001"])


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
    database = tmp_path / "combined.fa"
    cli.main([
        "--build", "-b", str(bold), "-u", str(unite), "-d", str(database),
        "-t", str(tmp_path / "work"), "-y", "k__Fungi",
    ])

    fasta = database.read_text()
    assert fasta.count(">") == 3
    assert "s__one" in fasta and "s__two" in fasta
    assert "U1|ITS|k__Fungi" in fasta
    assert database.with_suffix(".1.bt2l").is_file()
    manifest = json.loads(database.with_suffix(".manifest.json").read_text())
    assert len(manifest["sources"]["bold"]["sha256"]) == 64
    assert len(manifest["sources"]["unite"]["sha256"]) == 64
