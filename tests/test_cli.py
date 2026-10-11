import gzip
import json
import os
import shlex
import subprocess
from pathlib import Path

import pytest
import yaml

from barcodemapper import cli
from barcodemapper.reference_sources import Release


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
    assert config["barcodemapper_index_base"] == str(work / "reference/database")
    assert config["barcodemapper_db_headers"] == str(work / "reference/database.tsv")
    assert config["generate_index"] and config["generate_headers"]
    assert captured[0][0] == "snakemake"
    assert str(work) in captured[0]
    assert json.loads((work / "data/sample_to_reads1.json").read_text()) == {"sample": [str(reads / "sample_1.fq.gz")]}


@pytest.mark.parametrize("extension", ["bt2l", "bt2"])
def test_analysis_reuses_complete_current_database_artifacts(tmp_path, extension):
    database = _db(tmp_path)
    index_base = database.with_suffix("")
    for part in cli.INDEX_PARTS:
        Path(f"{index_base}.{part}.{extension}").write_text("index")
    headers = database.with_suffix(".tsv")
    headers.write_text("k__Fungi\n")

    assert cli._database_artifacts(database, tmp_path / "work", False) == (
        index_base, extension, False, headers, False
    )

    Path(f"{index_base}.rev.2.{extension}").unlink()
    chosen = cli._database_artifacts(database, tmp_path / "work", False)
    assert chosen[0] == tmp_path / "work/reference/database"
    assert chosen[2] is True
    assert chosen[3] == headers and chosen[4] is False

    Path(f"{index_base}.rev.2.{extension}").write_text("index")
    newer = max(path.stat().st_mtime_ns for path in tmp_path.iterdir()) + 1_000_000_000
    os.utime(database, ns=(newer, newer))
    chosen = cli._database_artifacts(database, tmp_path / "work", False)
    assert chosen[2] is True and chosen[4] is True
    assert chosen[3] == tmp_path / "work/reference/database.tsv"


def test_gzipped_database_uses_clean_artifact_names(tmp_path):
    database = tmp_path / "database.v1.fa.gz"
    with gzip.open(database, "wt") as handle:
        handle.write(">ref|ITS|k__Fungi\nACGT\n")
    base = tmp_path / "database.v1"
    for part in cli.INDEX_PARTS:
        Path(f"{base}.{part}.bt2l").write_text("index")
    headers = Path(f"{base}.tsv")
    headers.write_text("k__Fungi\n")
    assert cli._database_artifact(database, ".manifest.json") == tmp_path / "database.v1.manifest.json"
    assert cli._database_artifacts(database, tmp_path / "work", False) == (
        base, "bt2l", False, headers, False
    )


def test_preindexed_analysis_does_not_require_bowtie2_build(monkeypatch, tmp_path):
    database = _db(tmp_path)
    index_base = database.with_suffix("")
    for part in cli.INDEX_PARTS:
        Path(f"{index_base}.{part}.bt2l").write_text("index")
    database.with_suffix(".tsv").write_text("k__Fungi\n")
    reads = tmp_path / "reads"
    reads.mkdir()
    (reads / "sample_1.fq.gz").touch()
    (reads / "sample_2.fq.gz").touch()
    monkeypatch.setattr(cli.shutil, "which", lambda name: None if name == "bowtie2-build" else f"/mock/{name}")
    monkeypatch.setattr(cli, "run_snakemake", lambda *args: None)
    monkeypatch.setattr(cli, "_write_manifest", lambda *args: None)

    work = tmp_path / "work"
    cli.main(["-d", str(database), "-i", str(reads), "-o", str(tmp_path / "out.tsv"), "-t", str(work)])
    config = yaml.safe_load((work / "data/config.yaml").read_text())
    assert config["generate_index"] is False
    assert config["generate_headers"] is False
    assert config["barcodemapper_index_base"] == str(index_base)


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


def test_latest_build_skips_current_repairs_index_and_replaces_stale(monkeypatch, tmp_path):
    requested_database = tmp_path / "latest.fa"
    database = tmp_path / "latest.fa.gz"
    work = tmp_path / "work"
    release_id = ["BOLD_Public.09-Oct-2026"]
    downloads = []
    runs = []
    fail_next = [False]

    def discover():
        return {
            "bold": Release(release_id[0], "https://doi.org/10.5883/DP-Latest", "https://bold.example/fasta"),
            "unite": Release("2025-02-19", "https://doi.org/10.15156/BIO/3301232", "https://unite.example/data.tgz"),
        }

    def download(releases, directory):
        downloads.append(releases["bold"].identifier)
        directory.mkdir()
        bold = directory / "bold.fasta"
        unite = directory / "unite.fasta"
        bold.write_text(">A|COI-5P|x|Animalia\nACGT\n")
        unite.write_text(">record|U1|x|x|k__Fungi\nTGCA\n")
        return {"bold": bold, "unite": unite}, {
            "bold": {"release": releases["bold"].identifier, "doi": releases["bold"].doi, "sha256": cli._sha256(bold)},
            "unite": {"release": releases["unite"].identifier, "doi": releases["unite"].doi, "sha256": cli._sha256(unite)},
        }

    def run(tmp_dir, config_path, profile):
        config = yaml.safe_load(config_path.read_text())
        runs.append(config)
        if fail_next[0]:
            fail_next[0] = False
            raise RuntimeError("build failed")
        db = Path(config["barcodemapper_db"])
        if "bold_db" in config:
            with gzip.open(db, "wt") as handle:
                handle.write(f">{release_id[0]}|COI-5P|k__Animalia\nACGT\n>U1|ITS|k__Fungi\nTGCA\n")
        base = cli._database_base(db)
        if config["generate_index"]:
            for part in cli.INDEX_PARTS:
                Path(f"{base}.{part}.bt2l").write_text("index")
        if config["generate_headers"]:
            base.with_suffix(".tsv").write_text("k__Animalia\nk__Fungi\n")

    monkeypatch.setattr(cli, "discover_latest", discover)
    monkeypatch.setattr(cli, "download_sources", download)
    monkeypatch.setattr(cli, "run_snakemake", run)
    monkeypatch.setattr(cli.shutil, "which", lambda name: f"/mock/{name}")
    command = ["--build", "-d", str(requested_database), "-t", str(work)]
    cli.main(command)
    assert not requested_database.exists()
    first = database.read_bytes()
    manifest = json.loads(cli._database_base(database).with_suffix(".manifest.json").read_text())
    assert manifest["automatic_sources"] is True
    assert manifest["sources"]["bold"]["release"] == release_id[0]
    assert manifest["parameters"]["barcodemapper_db"] == str(database)
    assert "bold_db" not in manifest["parameters"]
    assert len(downloads) == len(runs) == 1

    cli.main(command)
    assert len(downloads) == len(runs) == 1

    cli._database_base(database).with_suffix(".rev.2.bt2l").unlink()
    cli.main(command)
    assert database.read_bytes() == first
    assert len(downloads) == 1 and len(runs) == 2
    assert cli._database_base(database).with_suffix(".rev.2.bt2l").is_file()

    release_id[0] = "BOLD_Public.16-Oct-2026"
    fail_next[0] = True
    with pytest.raises(RuntimeError, match="build failed"):
        cli.main(command)
    assert database.read_bytes() == first
    assert json.loads(cli._database_base(database).with_suffix(".manifest.json").read_text())["sources"]["bold"]["release"] == "BOLD_Public.09-Oct-2026"
    cli.main(command)
    assert len(downloads) == 3 and len(runs) == 4
    with gzip.open(database, "rt") as handle:
        assert release_id[0] in handle.read()
    assert list(tmp_path.glob(".latest.build-*")) == []


def test_latest_build_will_not_replace_untracked_database(monkeypatch, tmp_path):
    database = _db(tmp_path)
    monkeypatch.setattr(cli, "discover_latest", lambda: {
        "bold": Release("BOLD_Public.09-Oct-2026", "doi", "url"),
        "unite": Release("2025-02-19", "doi", "url"),
    })
    with pytest.raises(SystemExit) as error:
        cli.main(["--build", "-d", str(database)])
    assert error.value.code == 2
    assert database.is_file()


def test_manual_build_resumes_staged_fasta_and_promotes_only_when_complete(monkeypatch, tmp_path):
    bold = tmp_path / "bold.fa"
    unite = tmp_path / "unite.fa"
    bold.write_text(">bold\nACGT\n")
    unite.write_text(">unite\nTGCA\n")
    database = tmp_path / "combined.fa.gz"
    stage = tmp_path / ".combined.fa.gz.build"
    calls = []

    def run(tmp_dir, config_path, profile):
        config = yaml.safe_load(config_path.read_text())
        staged_db = Path(config["barcodemapper_db"])
        assert staged_db == stage / database.name
        assert tmp_dir == stage / "work"
        calls.append(staged_db.exists())
        if len(calls) == 1:
            staged_db.write_bytes(b"prepared")
            raise RuntimeError("index failed")
        for path in cli._build_outputs(staged_db)[1:]:
            path.write_text("artifact")

    monkeypatch.setattr(cli.shutil, "which", lambda name: f"/mock/{name}")
    monkeypatch.setattr(cli, "run_snakemake", run)
    command = ["--build", "-b", str(bold), "-u", str(unite), "-d", str(database)]

    with pytest.raises(RuntimeError, match="index failed"):
        cli.main(command)
    assert not database.exists()
    assert (stage / "stage.json").is_file()
    assert (stage / database.name).read_bytes() == b"prepared"

    cli.main(command)
    assert calls == [False, True]
    assert database.read_bytes() == b"prepared"
    assert not stage.exists()
    manifest = json.loads(cli._database_artifact(database, ".manifest.json").read_text())
    assert manifest["database"]["path"] == str(database)
    assert "automatic_sources" not in manifest
    assert manifest["sources"]["bold"]["sha256"] == cli._sha256(bold)


def _manual_setup(monkeypatch, tmp_path):
    bold = tmp_path / "bold.fa"
    unite = tmp_path / "unite.fa"
    bold.write_text(">bold\nACGT\n")
    unite.write_text(">unite\nTGCA\n")
    database = tmp_path / "combined.fa.gz"
    monkeypatch.setattr(cli.shutil, "which", lambda name: f"/mock/{name}")
    return bold, unite, database, tmp_path / ".combined.fa.gz.build"


def _command(bold, unite, database, *extra):
    return ["--build", "-b", str(bold), "-u", str(unite), "-d", str(database), *extra]


def _fail_with_partial(tmp_dir, config_path, profile):
    staged = Path(yaml.safe_load(config_path.read_text())["barcodemapper_db"])
    staged.write_bytes(b"stale")
    raise RuntimeError("boom")


def _write_all(tmp_dir, config_path, profile):
    staged = Path(yaml.safe_load(config_path.read_text())["barcodemapper_db"])
    staged.write_bytes(b"fresh")
    for path in cli._build_outputs(staged)[1:]:
        path.write_text("artifact")


@pytest.mark.parametrize("change", ["bold", "unite", "bold_retain", "unite_retain"])
def test_manual_build_discards_stale_stage_when_key_changes(monkeypatch, tmp_path, change):
    bold, unite, database, stage = _manual_setup(monkeypatch, tmp_path)
    monkeypatch.setattr(cli, "run_snakemake", _fail_with_partial)
    with pytest.raises(RuntimeError, match="boom"):
        cli.main(_command(bold, unite, database))
    assert (stage / database.name).read_bytes() == b"stale"
    (stage / "marker").write_text("x")

    extra = []
    if change == "bold":
        bold.write_text(">bold\nGGGG\n")
    elif change == "unite":
        unite.write_text(">unite\nCCCC\n")
    elif change == "bold_retain":
        extra = ["-x", "k__Other"]
    else:
        extra = ["-y", "c__Other"]

    seen = []

    def run(tmp_dir, config_path, profile):
        staged = Path(yaml.safe_load(config_path.read_text())["barcodemapper_db"])
        seen.append((staged.exists(), (stage / "marker").exists()))
        _write_all(tmp_dir, config_path, profile)

    monkeypatch.setattr(cli, "run_snakemake", run)
    cli.main(_command(bold, unite, database, *extra))
    assert seen == [(False, False)]
    assert database.read_bytes() == b"fresh"


def test_manual_build_keeps_matching_stage(monkeypatch, tmp_path):
    bold, unite, database, stage = _manual_setup(monkeypatch, tmp_path)
    monkeypatch.setattr(cli, "run_snakemake", _fail_with_partial)
    with pytest.raises(RuntimeError):
        cli.main(_command(bold, unite, database))
    seen = []

    def run(tmp_dir, config_path, profile):
        seen.append((stage / database.name).exists())
        _write_all(tmp_dir, config_path, profile)

    monkeypatch.setattr(cli, "run_snakemake", run)
    cli.main(_command(bold, unite, database))
    assert seen == [True]


def test_manual_build_refuses_complete_database_and_leaves_it_untouched(monkeypatch, tmp_path):
    bold, unite, database, stage = _manual_setup(monkeypatch, tmp_path)
    database.write_bytes(b"existing")
    for path in cli._build_outputs(database)[1:]:
        path.write_text("keep")
    before = {p.name: p.read_bytes() for p in tmp_path.iterdir()}

    def run(*args):
        raise AssertionError("snakemake must not run")

    monkeypatch.setattr(cli, "run_snakemake", run)
    with pytest.raises(SystemExit) as error:
        cli.main(_command(bold, unite, database))
    assert error.value.code == 2
    assert {p.name: p.read_bytes() for p in tmp_path.iterdir()} == before
    assert not stage.exists()


def test_manual_build_refuses_manifest_only_target(monkeypatch, tmp_path, capsys):
    bold, unite, database, stage = _manual_setup(monkeypatch, tmp_path)
    database.write_bytes(b"existing")
    cli._database_artifact(database, ".manifest.json").write_text("{}")
    monkeypatch.setattr(cli, "run_snakemake", lambda *a: pytest.fail("ran"))
    with pytest.raises(SystemExit) as error:
        cli.main(_command(bold, unite, database))
    assert error.value.code == 2
    assert "database artifacts already exist" in capsys.readouterr().err
    assert database.read_bytes() == b"existing"


def test_manual_build_lone_fasta_reports_incomplete(monkeypatch, tmp_path, capsys):
    bold, unite, database, stage = _manual_setup(monkeypatch, tmp_path)
    database.write_bytes(b"lone")
    monkeypatch.setattr(cli, "run_snakemake", lambda *a: pytest.fail("ran"))
    with pytest.raises(SystemExit) as error:
        cli.main(_command(bold, unite, database))
    assert error.value.code == 2
    err = capsys.readouterr().err
    assert f"{database} exists without a manifest or index" in err
    assert "if it is a failed earlier build, remove it, otherwise choose a new -d path" in err
    assert database.read_bytes() == b"lone"
    assert not stage.exists()


def test_manual_build_failure_leaves_no_final_named_artifact(monkeypatch, tmp_path):
    bold, unite, database, stage = _manual_setup(monkeypatch, tmp_path)
    monkeypatch.setattr(cli, "run_snakemake", _fail_with_partial)
    with pytest.raises(RuntimeError, match="boom"):
        cli.main(_command(bold, unite, database))
    for path in cli._build_outputs(database):
        assert not path.exists()
    assert not cli._database_artifact(database, ".manifest.json").exists()
    assert stage.is_dir()


def test_manual_build_success_removes_stage_and_writes_manifest(monkeypatch, tmp_path):
    bold, unite, database, stage = _manual_setup(monkeypatch, tmp_path)
    original_sha256 = cli._sha256
    hashed = []

    def record_hash(path):
        hashed.append(path)
        return original_sha256(path)

    monkeypatch.setattr(cli, "_sha256", record_hash)
    monkeypatch.setattr(cli, "run_snakemake", _write_all)
    cli.main(_command(bold, unite, database))
    assert not stage.exists()
    for path in cli._build_outputs(database):
        assert path.is_file()
    manifest_path = cli._database_artifact(database, ".manifest.json")
    manifest = json.loads(manifest_path.read_text())
    assert manifest["database"]["path"] == str(database)
    assert manifest["sources"]["unite"]["sha256"] == original_sha256(unite)
    assert hashed.count(bold) == hashed.count(unite) == 1
    assert sorted(p.name for p in tmp_path.iterdir() if p.name.startswith(".")) == []


def test_manual_build_missing_staged_artifact_raises_and_keeps_stage(monkeypatch, tmp_path):
    bold, unite, database, stage = _manual_setup(monkeypatch, tmp_path)

    def run(tmp_dir, config_path, profile):
        _write_all(tmp_dir, config_path, profile)
        staged = Path(yaml.safe_load(config_path.read_text())["barcodemapper_db"])
        cli._build_outputs(staged)[-1].unlink()

    monkeypatch.setattr(cli, "run_snakemake", run)
    with pytest.raises(RuntimeError, match="missing staged artifacts"):
        cli.main(_command(bold, unite, database))
    assert stage.is_dir() and (stage / database.name).is_file()
    assert not database.exists()
    assert not cli._database_artifact(database, ".manifest.json").exists()


def test_manual_build_empty_staged_artifact_is_treated_as_missing(monkeypatch, tmp_path):
    bold, unite, database, stage = _manual_setup(monkeypatch, tmp_path)

    def run(tmp_dir, config_path, profile):
        _write_all(tmp_dir, config_path, profile)
        staged = Path(yaml.safe_load(config_path.read_text())["barcodemapper_db"])
        cli._build_outputs(staged)[1].write_text("")

    monkeypatch.setattr(cli, "run_snakemake", run)
    with pytest.raises(RuntimeError, match="missing staged artifacts"):
        cli.main(_command(bold, unite, database))
    assert stage.is_dir()


@pytest.mark.parametrize("requested_gzip", [False, True])
def test_manual_build_reports_uncompressed_lone_fasta(monkeypatch, tmp_path, capsys, requested_gzip):
    bold, unite, _, _ = _manual_setup(monkeypatch, tmp_path)
    uncompressed = tmp_path / "orphan.fa"
    uncompressed.write_bytes(b"keep")
    database = tmp_path / "orphan.fa.gz" if requested_gzip else uncompressed

    with pytest.raises(SystemExit) as error:
        cli.main(_command(bold, unite, database))
    assert error.value.code == 2
    err = capsys.readouterr().err
    assert f"{uncompressed} exists without a manifest or index" in err
    assert "if it is a failed earlier build, remove it, otherwise choose a new -d path" in err
    assert uncompressed.read_bytes() == b"keep"
    assert not (tmp_path / ".orphan.fa.gz.build").exists()


@pytest.mark.parametrize("suffix", [".manifest.json", ".tsv", ".1.bt2l", ".rev.2.bt2"])
def test_manual_build_refuses_stray_final_artifact_without_fasta(monkeypatch, tmp_path, capsys, suffix):
    bold, unite, database, stage = _manual_setup(monkeypatch, tmp_path)
    artifact = cli._database_artifact(database, suffix)
    artifact.write_bytes(b"keep")

    with pytest.raises(SystemExit) as error:
        cli.main(_command(bold, unite, database))
    assert error.value.code == 2
    assert "database artifacts already exist" in capsys.readouterr().err
    assert artifact.read_bytes() == b"keep"
    assert not database.exists() and not stage.exists()


def test_manual_build_different_tmpdir_discards_stage(monkeypatch, tmp_path):
    bold, unite, database, stage = _manual_setup(monkeypatch, tmp_path)
    first_work = tmp_path / "first-work"
    second_work = tmp_path / "second-work"
    monkeypatch.setattr(cli, "run_snakemake", _fail_with_partial)
    with pytest.raises(RuntimeError, match="boom"):
        cli.main(_command(bold, unite, database, "-t", str(first_work)))
    (stage / "marker").write_text("stale")
    seen = []

    def run(tmp_dir, config_path, profile):
        seen.append((tmp_dir, (stage / database.name).exists(), (stage / "marker").exists()))
        _write_all(tmp_dir, config_path, profile)

    monkeypatch.setattr(cli, "run_snakemake", run)
    cli.main(_command(bold, unite, database, "-t", str(second_work)))
    assert seen == [(second_work, False, False)]
    assert not stage.exists()


def test_manual_build_stage_name_includes_full_database_filename(monkeypatch, tmp_path):
    bold, unite, database, first_stage = _manual_setup(monkeypatch, tmp_path)
    other = tmp_path / "combined.fasta.gz"
    other_stage = tmp_path / ".combined.fasta.gz.build"
    monkeypatch.setattr(cli, "run_snakemake", _fail_with_partial)
    for target in (database, other):
        with pytest.raises(RuntimeError, match="boom"):
            cli.main(_command(bold, unite, target))
    assert (first_stage / database.name).read_bytes() == b"stale"
    assert (other_stage / other.name).read_bytes() == b"stale"


def test_manual_build_reused_stage_prints_unlock_hint_for_failed_snakemake(monkeypatch, tmp_path, capsys):
    bold, unite, database, stage = _manual_setup(monkeypatch, tmp_path)
    work = tmp_path / "work files"
    command = _command(bold, unite, database, "-t", str(work), "--slurm")

    def fail(*args):
        raise subprocess.CalledProcessError(1, ["snakemake"])

    monkeypatch.setattr(cli, "run_snakemake", fail)
    with pytest.raises(subprocess.CalledProcessError):
        cli.main(command)
    assert "barcodemapper --unlock" not in capsys.readouterr().err
    with pytest.raises(subprocess.CalledProcessError):
        cli.main(command)
    assert f"barcodemapper --unlock -t {shlex.quote(str(work))} --slurm" in capsys.readouterr().err
    assert stage.is_dir()


@pytest.mark.parametrize("cleanup_fails", [False, True])
def test_manual_build_promotion_failure_keeps_stage_and_rolls_back(monkeypatch, tmp_path, cleanup_fails):
    bold, unite, database, stage = _manual_setup(monkeypatch, tmp_path)
    monkeypatch.setattr(cli, "run_snakemake", _write_all)
    original_replace = Path.replace

    def fail_second_index(source, target):
        if source == stage / "combined.2.bt2l" and target == tmp_path / "combined.2.bt2l":
            raise OSError("promotion failed")
        return original_replace(source, target)

    monkeypatch.setattr(Path, "replace", fail_second_index)
    if cleanup_fails:
        original_rmtree = cli.shutil.rmtree

        def fail_cleanup(path, *args, **kwargs):
            if path == stage / "previous":
                raise OSError("cleanup failed")
            return original_rmtree(path, *args, **kwargs)

        monkeypatch.setattr(cli.shutil, "rmtree", fail_cleanup)
    with pytest.raises(OSError, match="promotion failed") as error:
        cli.main(_command(bold, unite, database))
    assert not database.exists()
    assert not cli._database_artifact(database, ".manifest.json").exists()
    assert (stage / database.name).is_file()
    assert (stage / "previous").exists() is cleanup_fails
    if cleanup_fails:
        assert "cleanup failed" in " ".join(error.value.__notes__)


def test_promote_database_clears_stale_previous_and_moves_manifest_last(monkeypatch, tmp_path):
    stage = tmp_path / "stage"
    stage.mkdir()
    staged = stage / "combined.fa.gz"
    destination = tmp_path / staged.name
    for path in cli._build_outputs(staged) + [cli._database_artifact(staged, ".manifest.json")]:
        path.write_text("new")
    previous = stage / "previous"
    previous.mkdir()
    (previous / "stale").write_text("old")
    moved = []
    original_replace = Path.replace

    def record_replace(source, target):
        moved.append((source, target))
        return original_replace(source, target)

    monkeypatch.setattr(Path, "replace", record_replace)
    cli._promote_database(staged, destination)
    assert moved[-1] == (cli._database_artifact(staged, ".manifest.json"),
                          cli._database_artifact(destination, ".manifest.json"))
    assert not previous.exists()
    assert all(path.is_file() for path in cli._build_outputs(destination))
