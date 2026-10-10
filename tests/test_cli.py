import json
from pathlib import Path

import pytest
import yaml

from barcodemapper import cli


def _db(root: Path) -> Path:
    path = root / "database.fa"
    path.write_text(">ref|ITS|k__Fungi\nACGT\n")
    return path


def test_unpaired_read_is_rejected(tmp_path):
    reads = tmp_path / "reads"
    reads.mkdir()
    (reads / "sample_1.fq.gz").touch()
    with pytest.raises(ValueError, match="missing R2"):
        cli.inputdir_to_samples(reads, tmp_path / "work")


def test_mate_number_and_duplicate_samples_are_checked(tmp_path):
    first = tmp_path / "a"
    second = tmp_path / "b"
    first.mkdir()
    second.mkdir()
    for path in (first / "sample_1.fq.gz", second / "sample_1.fq.gz", first / "sample_2.fq.gz"):
        path.touch()
    with pytest.raises(ValueError, match="Expected an R2"):
        cli.inputlist_to_samples(str(first / "sample_1.fq.gz"), str(second / "sample_1.fq.gz"), tmp_path / "work")
    with pytest.raises(ValueError, match="multiple files per mate"):
        cli.inputlist_to_samples(
            f"{first / 'sample_1.fq.gz'},{second / 'sample_1.fq.gz'}",
            str(first / "sample_2.fq.gz"),
            tmp_path / "work",
        )


def test_reserved_sample_name_is_rejected(tmp_path):
    reads = tmp_path / "reads"
    reads.mkdir()
    (reads / "Taxonomy_1.fq.gz").touch()
    (reads / "Taxonomy_2.fq.gz").touch()
    with pytest.raises(ValueError, match="reserved for an output column"):
        cli.inputdir_to_samples(reads, tmp_path / "work")


def test_build_and_analysis_require_their_inputs(tmp_path):
    database = _db(tmp_path)
    with pytest.raises(SystemExit) as missing_output:
        cli.main(["-d", str(database), "-i", str(tmp_path)])
    assert missing_output.value.code == 2
    with pytest.raises(SystemExit) as wrong_build_args:
        cli.main(["--build", "-d", str(tmp_path / "new.fa"), "-o", "out.tsv"])
    assert wrong_build_args.value.code == 2


def test_numeric_options_and_paths_with_spaces(monkeypatch, tmp_path):
    root = tmp_path / "space in path"
    root.mkdir()
    database = _db(root)
    reads = root / "reads"
    reads.mkdir()
    (reads / "sample_1.fq.gz").touch()
    (reads / "sample_2.fq.gz").touch()
    work = root / "working files"
    output = root / "result table.tsv"
    captured = []
    monkeypatch.setattr(cli.shutil, "which", lambda name: f"/mock/{name}")
    monkeypatch.setattr(cli.subprocess, "run", lambda command, **kwargs: captured.append(command))
    monkeypatch.setattr(cli, "_write_manifest", lambda *args: None)
    cli.main(["-d", str(database), "-i", str(reads), "-o", str(output), "-t", str(work), "-m", "2", "-c", "100"])
    config = yaml.safe_load((work / "data/config.yaml").read_text())
    assert config["max_edit_percent"] == 2.0
    assert config["min_coverage"] == 100
    assert config["barcodemapper_db"] == str(database)
    assert captured[0][0] == "snakemake"
    assert str(work) in captured[0]
    assert json.loads((work / "data/sample_to_reads1.json").read_text()) == {"sample": [str(reads / "sample_1.fq.gz")]}


def test_nan_edit_limit_is_rejected(tmp_path):
    database = _db(tmp_path)
    with pytest.raises(SystemExit) as error:
        cli.main(["-d", str(database), "-i", str(tmp_path), "-o", str(tmp_path / "out.tsv"), "-m", "nan"])
    assert error.value.code == 2


def test_analysis_rejects_malformed_database_before_mapping(tmp_path, capsys):
    database = tmp_path / "invalid.fa"
    database.write_text(">ref|ITS\nACGT\n")
    with pytest.raises(SystemExit) as error:
        cli.main(["-d", str(database), "-i", str(tmp_path), "-o", str(tmp_path / "out.tsv")])
    assert error.value.code == 2
    assert "invalid database FASTA" in capsys.readouterr().err
