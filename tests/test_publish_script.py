"""Offline checks for the maintainer Zenodo publishing script."""

import gzip
import getpass
import hashlib
import importlib.util
import io
import json
import os
import subprocess
import tarfile
from pathlib import Path
from urllib.request import Request

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "publish_zenodo.py"
BUILD_SCRIPT = SCRIPT.with_name("build_and_publish.sh")
STEM = "bmdb_20261021"
TOKEN = "recognisable-secret-token-123"


@pytest.fixture
def publisher():
    spec = importlib.util.spec_from_file_location("publish_zenodo", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def database(tmp_path):
    folder = tmp_path / STEM
    folder.mkdir()
    fasta = gzip.compress(b">B1|COI-5P|k__Animalia\nACGT\n>U1|ITS|k__Fungi\nTGCA\n")
    (folder / f"{STEM}.fa.gz").write_bytes(fasta)
    (folder / f"{STEM}.tsv").write_text("taxon\n", encoding="utf-8")
    manifest = {
        "automatic_sources": True,
        "database": {"path": str(folder / f"{STEM}.fa.gz"), "sha256": hashlib.sha256(fasta).hexdigest()},
        "parameters": {"bold_retain": "k__Animalia", "unite_retain": "k__Fungi"},
        "sources": {
            "bold": {"release": "BOLD_Public.09-Oct-2026", "doi": "https://doi.org/10.5883/DP-Latest",
                     "sha256": "bold-hash"},
            "unite": {"release": "2025-02-19", "doi": "https://doi.org/10.15156/BIO/3301232",
                      "sha256": "unite-hash"},
        },
    }
    (folder / f"{STEM}.manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    for part in ("1", "2", "3", "4", "rev.1", "rev.2"):
        (folder / f"{STEM}.{part}.bt2l").write_bytes(b"index-" + part.encode())
    return folder


def _run(publisher, database, *options):
    publisher.main(["--database", str(database), *options])


def _fake_api(publisher, monkeypatch):
    calls = []
    uploaded = {}

    def request(method, url, token, host, *, payload=None, file_path=None):
        calls.append((method, url, token, host, payload))
        if file_path is not None:
            uploaded[file_path.name] = file_path.read_bytes()
            return {"checksum": "md5:" + hashlib.md5(uploaded[file_path.name]).hexdigest()}
        if method == "POST" and url.endswith("/api/deposit/depositions"):
            return {"id": 123, "links": {"bucket": f"https://{host}/api/files/bucket"}}
        if method == "POST" and url.endswith("/actions/publish"):
            return {"id": 123, "conceptrecid": 100, "doi": "10.5281/zenodo.123",
                    "conceptdoi": "10.5281/zenodo.100"}
        return {}

    monkeypatch.setattr(publisher, "_request", request)
    return calls, uploaded


def test_dry_run_needs_no_token_or_network(publisher, database, monkeypatch, capsys):
    monkeypatch.delenv("ZENODO_TOKEN", raising=False)
    monkeypatch.setattr(publisher, "_request", lambda *args, **kwargs: pytest.fail("unexpected HTTP"))
    _run(publisher, database, "--dry-run")
    output = capsys.readouterr().out
    assert "Planned files:" in output and "Metadata:" in output
    assert f"{STEM}.fa.gz" in output and '"version": "20261021"' in output
    assert "moving alias" in output and "BOLD_Public.09-Oct-2026" in output


@pytest.mark.parametrize("publish", [False, True])
def test_publish_is_opt_in(publisher, database, monkeypatch, capsys, publish):
    monkeypatch.setenv("ZENODO_TOKEN", TOKEN)
    calls, uploaded = _fake_api(publisher, monkeypatch)
    _run(publisher, database, *(["--publish"] if publish else []))
    actions = [url for method, url, *_ in calls if method == "POST" and url.endswith("/actions/publish")]
    assert len(actions) == int(publish)
    assert set(uploaded) == {f"{STEM}.fa.gz", f"{STEM}.tsv", f"{STEM}.manifest.json",
                             "ATTRIBUTION.txt", "SHA256SUMS"}
    output = capsys.readouterr().out
    assert "Draft: https://zenodo.org/deposit/123" in output
    if publish:
        assert "DEFAULT_RECORD = 100" in output
    assert all(token == TOKEN for _, _, token, _, _ in calls)


def test_missing_token_is_clear(publisher, database, monkeypatch, capsys):
    monkeypatch.delenv("ZENODO_TOKEN", raising=False)
    monkeypatch.setattr(publisher, "_request", lambda *args, **kwargs: pytest.fail("unexpected HTTP"))
    with pytest.raises(SystemExit):
        _run(publisher, database)
    assert "ZENODO_TOKEN is required" in capsys.readouterr().err


def test_token_is_redacted_from_errors(publisher, database, monkeypatch, capsys):
    monkeypatch.setenv("ZENODO_TOKEN", TOKEN)

    def fail(*args, **kwargs):
        raise ValueError("request failed: " + TOKEN)

    monkeypatch.setattr(publisher, "_request", fail)
    with pytest.raises(SystemExit) as error:
        _run(publisher, database)
    output = capsys.readouterr()
    assert TOKEN not in output.out + output.err + str(error.value)
    assert "[REDACTED]" in output.err


def test_cross_host_draft_link_is_rejected(publisher, database, monkeypatch, capsys):
    monkeypatch.setenv("ZENODO_TOKEN", TOKEN)
    calls = []

    def request(method, url, token, host, **kwargs):
        calls.append(url)
        if url.endswith("/api/records/12"):
            return {"id": 34}
        return {"links": {"latest_draft": "https://other.example/api/deposit/depositions/12"}}

    monkeypatch.setattr(publisher, "_request", request)
    with pytest.raises(SystemExit):
        _run(publisher, database, "--record", "12")
    assert len(calls) == 2 and "unexpected host" in capsys.readouterr().err


def test_request_rejects_cross_host_before_network(publisher, monkeypatch):
    monkeypatch.setattr(publisher, "build_opener", lambda *args: pytest.fail("unexpected HTTP"))
    with pytest.raises(ValueError, match="unexpected host"):
        publisher._request("GET", "https://other.example/api/deposit/depositions", TOKEN, "zenodo.org")


def test_redirect_keeps_token_on_zenodo_host_only(publisher):
    handler = publisher._SameHostRedirect("zenodo.org")
    request = Request("https://zenodo.org/api/records/12",
                      headers={"Authorization": "Bearer " + TOKEN})
    redirected = handler.redirect_request(
        request, None, 302, "Found", {}, "https://zenodo.org/api/records/34")
    assert redirected.get_header("Authorization") == "Bearer " + TOKEN
    with pytest.raises(ValueError, match="unexpected host"):
        handler.redirect_request(request, None, 302, "Found", {},
                                 "https://other.example/api/records/34")


def test_new_version_resolves_concept_and_removes_old_files(publisher, database, monkeypatch):
    monkeypatch.setenv("ZENODO_TOKEN", TOKEN)
    calls = []

    def request(method, url, token, host, *, payload=None, file_path=None):
        calls.append((method, url))
        if url.endswith("/api/records/12"):
            return {"id": 34}
        if url.endswith("/actions/newversion"):
            return {"links": {"latest_draft": "https://zenodo.org/api/deposit/depositions/56"}}
        if method == "GET" and url.endswith("/depositions/56"):
            return {"id": 56, "files": [{"id": "old-file-id"}],
                    "links": {"bucket": "https://zenodo.org/api/files/new-bucket"}}
        if file_path is not None:
            return {"checksum": "md5:" + hashlib.md5(file_path.read_bytes()).hexdigest()}
        return {}

    monkeypatch.setattr(publisher, "_request", request)
    _run(publisher, database, "--record", "12")
    assert calls[:4] == [
        ("GET", "https://zenodo.org/api/records/12"),
        ("POST", "https://zenodo.org/api/deposit/depositions/34/actions/newversion"),
        ("GET", "https://zenodo.org/api/deposit/depositions/56"),
        ("DELETE", "https://zenodo.org/api/deposit/depositions/56/files/old-file-id"),
    ]


def test_manifest_hash_mismatch_aborts_before_http(publisher, database, monkeypatch, capsys):
    monkeypatch.setenv("ZENODO_TOKEN", TOKEN)
    monkeypatch.setattr(publisher, "_request", lambda *args, **kwargs: pytest.fail("unexpected HTTP"))
    (database / f"{STEM}.fa.gz").write_bytes(b"wrong")
    with pytest.raises(SystemExit):
        _run(publisher, database)
    assert "SHA-256" in capsys.readouterr().err


def test_manual_sources_are_rejected(publisher, database, monkeypatch, capsys):
    monkeypatch.delenv("ZENODO_TOKEN", raising=False)
    path = database / f"{STEM}.manifest.json"
    manifest = json.loads(path.read_text())
    manifest["automatic_sources"] = False
    path.write_text(json.dumps(manifest))
    with pytest.raises(SystemExit):
        _run(publisher, database, "--dry-run")
    assert "automatic source provenance" in capsys.readouterr().err


def test_checksums_and_index_archive_match_fetch_contract(publisher, database, monkeypatch):
    monkeypatch.setenv("ZENODO_TOKEN", TOKEN)
    _, uploaded = _fake_api(publisher, monkeypatch)
    _run(publisher, database, "--include-index")
    expected = {f"{STEM}.fa.gz", f"{STEM}.tsv", f"{STEM}.manifest.json",
                "ATTRIBUTION.txt", f"{STEM}.index.tar", "SHA256SUMS"}
    assert set(uploaded) == expected
    sums = uploaded["SHA256SUMS"].decode().splitlines()
    assert set(sums) == {f"{hashlib.sha256(data).hexdigest()}  {name}"
                         for name, data in uploaded.items() if name != "SHA256SUMS"}
    with tarfile.open(fileobj=io.BytesIO(uploaded[f"{STEM}.index.tar"]), mode="r:") as archive:
        assert [member.name for member in archive] == [f"{STEM}.{part}.bt2l"
                                                      for part in ("1", "2", "3", "4", "rev.1", "rev.2")]
        assert all(member.isfile() for member in archive.getmembers())
    portable = json.loads(uploaded[f"{STEM}.manifest.json"])
    assert portable["database"]["path"] == f"{STEM}.fa.gz"
    assert portable["origin"] == "zenodo"
    assert "BOLD: CC BY-SA 4.0" in uploaded["ATTRIBUTION.txt"].decode()
    assert "UNITE: CC BY 4.0" in uploaded["ATTRIBUTION.txt"].decode()


def test_bold_alias_requires_release(publisher, database, capsys):
    path = database / f"{STEM}.manifest.json"
    manifest = json.loads(path.read_text())
    manifest["sources"]["bold"].pop("release")
    path.write_text(json.dumps(manifest))
    with pytest.raises(SystemExit):
        _run(publisher, database, "--dry-run")
    assert "release" in capsys.readouterr().err


def test_publisher_scrubs_embedded_paths_and_user(publisher, database, capsys):
    path = database / f"{STEM}.manifest.json"
    manifest = json.loads(path.read_text())
    manifest["tools"] = {"bowtie2": f"built by {getpass.getuser()} at /projects/alberdilab/scratch/jpl786/work"}
    path.write_text(json.dumps(manifest))
    _run(publisher, database, "--dry-run")
    assert getpass.getuser() not in capsys.readouterr().out


def test_publisher_rejects_unscrubbed_portable_manifest(publisher, database, monkeypatch, capsys):
    monkeypatch.setattr(publisher, "portable_manifest", lambda manifest: {
        "database": {"sha256": manifest["database"]["sha256"]},
        "tools": {"bad": getpass.getuser()},
    })
    with pytest.raises(SystemExit):
        _run(publisher, database, "--dry-run")
    assert "local user name" in capsys.readouterr().err


def test_slurm_dry_run_uses_submit_dir_without_writable_db_dir(database, tmp_path):
    environment = dict(os.environ, SLURM_JOB_ID="123", SLURM_SUBMIT_DIR=str(SCRIPT.parents[1]))
    environment.pop("BM_REPO", None)
    environment.pop("ZENODO_TOKEN", None)
    result = subprocess.run(["bash", str(BUILD_SCRIPT), "--db-dir", str(tmp_path / "absent"),
                             "--reuse", str(database), "--dry-run"], env=environment,
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "Planned files:" in result.stdout


def test_slurm_reports_missing_publisher_early(database, tmp_path):
    environment = dict(os.environ, SLURM_JOB_ID="123", SLURM_SUBMIT_DIR=str(tmp_path))
    environment.pop("BM_REPO", None)
    result = subprocess.run(["bash", str(BUILD_SCRIPT), "--db-dir", str(tmp_path / "absent"),
                             "--reuse", str(database), "--dry-run"], env=environment,
                            capture_output=True, text=True)
    assert result.returncode == 2
    assert "publisher not found" in result.stderr


@pytest.mark.parametrize("token_in_file", [False, True])
def test_live_build_hides_token_and_uses_reported_folder(database, tmp_path, token_in_file):
    scratch = tmp_path / "scratch"
    db_dir = scratch / "dbs"
    built = db_dir / "bmdb_20261020"
    built.mkdir(parents=True)
    (built / "bmdb_20261020.manifest.json").write_text("{}")
    (built / "bmdb_20261020.fa.gz").write_bytes(b"database")
    newer = db_dir / STEM
    newer.mkdir()
    (newer / f"{STEM}.manifest.json").write_text("{}")
    (newer / f"{STEM}.fa.gz").write_bytes(b"database")
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    (fake_bin / "barcodemapper").write_text(
        "#!/usr/bin/env bash\n[[ -z ${ZENODO_TOKEN:-} ]] || exit 9\n"
        "printf 'Built latest database: %s\\n' \"$TEST_BUILT\"\n")
    (fake_bin / "python").write_text(
        "#!/usr/bin/env bash\n[[ -n ${ZENODO_TOKEN:-} ]] || exit 9\n"
        "printf '%s\\n' \"$*\" > \"$TEST_MARKER\"\n")
    for name in ("snakemake", "bowtie2", "bowtie2-build"):
        (fake_bin / name).write_text("#!/usr/bin/env bash\nexit 0\n")
    for path in fake_bin.iterdir():
        path.chmod(0o755)
    conda_sh = tmp_path / "conda.sh"
    conda_sh.write_text(': "$UNSET_DURING_CONDA_LOAD"\nconda() { : "$UNSET_DURING_CONDA_ACTIVATE"; }\n')
    marker = tmp_path / "publisher-args"
    environment = dict(os.environ, PATH=f"{fake_bin}:{os.environ['PATH']}", BM_REPO=str(SCRIPT.parents[1]),
                       BM_SCRATCH_ROOT=str(scratch), BM_CONDA_SH=str(conda_sh),
                       ZENODO_TOKEN="test-token", TEST_BUILT=str(built), TEST_MARKER=str(marker))
    if token_in_file:
        job_home = tmp_path / "job-home"
        job_home.mkdir()
        (job_home / ".tokens").write_text("ZENODO_TOKEN=test-token\n")
        environment["HOME"] = str(job_home)
        environment.pop("ZENODO_TOKEN")
    result = subprocess.run(["bash", str(BUILD_SCRIPT), "--db-dir", str(db_dir), "--env", "fake"],
                            env=environment, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert f"Database: {built}" in result.stdout
    assert f"--database {built}" in marker.read_text()
    assert "test-token" not in result.stdout + result.stderr


def test_live_build_requires_explicit_environment(tmp_path):
    environment = dict(os.environ)
    environment.pop("BM_ENV", None)
    result = subprocess.run(["bash", str(BUILD_SCRIPT), "--db-dir", str(tmp_path)],
                            env=environment, capture_output=True, text=True)
    assert result.returncode == 2
    assert "--env is required" in result.stderr
