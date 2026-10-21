"""Offline CLI checks for installing published Zenodo databases."""

import gzip
import hashlib
import http.client
import io
import json
import tarfile
from pathlib import Path

import pytest

from barcodemapper import cli, zenodo
from barcodemapper.databases import find_databases, is_outdated
from barcodemapper.reference_sources import Release


VERSION = "20261021"
STEM = f"bmdb_{VERSION}"
DOI = "10.5281/zenodo.123"


def _archive(extra=None, stem=STEM):
    output = io.BytesIO()
    with tarfile.open(fileobj=output, mode="w") as archive:
        for part in cli.INDEX_PARTS:
            name = f"{stem}.{part}.bt2l"
            if extra is not None and extra[0] == name:
                continue
            data = b"index"
            member = tarfile.TarInfo(name)
            member.size = len(data)
            archive.addfile(member, io.BytesIO(data))
        if extra is not None:
            member = tarfile.TarInfo(extra[0])
            member.type = extra[1]
            member.size = 0
            archive.addfile(member)
    return output.getvalue()


def _record(monkeypatch, *, index=True, bad_sha=None, archive=None, edit=None, version=VERSION):
    stem = f"bmdb_{version}"
    fasta = gzip.compress(b">B1|COI-5P|k__Animalia\nACGT\n>U1|ITS|k__Fungi\nTGCA\n")
    manifest = {
        "created_utc": "2026-10-21T00:00:00+00:00",
        "automatic_sources": True,
        "database": {"path": f"/private/build/{stem}.fa.gz", "sha256": hashlib.sha256(fasta).hexdigest()},
        "parameters": {"bold_retain": "k__Animalia", "unite_retain": "k__Viridiplantae,c__Agaricomycetes,c__Pezizomycetes",
                       "tmp_dir": "/private/build/work"},
        "sources": {"bold": {"release": "BOLD_1", "doi": "d1", "path": "/private/bold.fa"},
                    "unite": {"release": "UNITE_1", "doi": "d2", "path": "/private/unite.fa"}},
    }
    payloads = {
        f"{stem}.fa.gz": fasta,
        f"{stem}.tsv": b"taxon\n",
        f"{stem}.manifest.json": json.dumps(manifest).encode(),
        "ATTRIBUTION.txt": b"CC BY-SA 4.0\n",
    }
    if index:
        payloads[f"{stem}.index.tar"] = archive if archive is not None else _archive(stem=stem)
    sums = {name: hashlib.sha256(data).hexdigest() for name, data in payloads.items()}
    if bad_sha is not None:
        sums[bad_sha] = "0" * 64
    payloads["SHA256SUMS"] = "".join(f"{digest}  {name}\n" for name, digest in sums.items()).encode()
    metadata = {
        "id": 123, "doi": DOI, "links": {"html": "https://zenodo.org/records/123"},
        "metadata": {"version": version, "license": "cc-by-sa-4.0"},
        "files": [
            {"key": name, "size": len(data), "checksum": f"md5:{hashlib.md5(data).hexdigest()}",
             "links": {"self": f"https://zenodo.org/api/files/123/{name}"}}
            for name, data in payloads.items()
        ],
    }
    if edit is not None:
        old_sums = payloads["SHA256SUMS"]
        edit(payloads, metadata)
        if payloads["SHA256SUMS"] == old_sums:
            payloads["SHA256SUMS"] = "".join(
                f"{hashlib.sha256(data).hexdigest()}  {name}\n"
                for name, data in payloads.items() if name != "SHA256SUMS"
            ).encode()
        for item in metadata["files"]:
            if item["key"] in payloads:
                data = payloads[item["key"]]
                item["size"] = len(data)
                item["checksum"] = f"md5:{hashlib.md5(data).hexdigest()}"
    calls = []

    def open_zenodo(url, **kwargs):
        calls.append(url)
        if url.endswith("/api/records/1/versions?size=100"):
            return io.BytesIO(json.dumps({"hits": {"hits": [
                {"id": 123, "metadata": {"version": version}},
            ]}}).encode())
        if url.endswith(("/api/records/1", "/api/records/123")):
            return io.BytesIO(json.dumps(metadata).encode())
        return io.BytesIO(payloads[url.rsplit("/", 1)[-1]])

    monkeypatch.setattr(zenodo, "_open", open_zenodo)
    return calls


def _fetch(root, *extra):
    cli.main(["--fetch", "--zenodo-record", "1", "--db-dir", str(root), *extra])


def test_fetch_lists_verifies_and_keeps_portable_provenance(monkeypatch, tmp_path, capsys):
    calls = _record(monkeypatch)
    _fetch(tmp_path)
    folder = tmp_path / STEM
    fasta = folder / f"{STEM}.fa.gz"
    assert find_databases(tmp_path) == [(STEM, fasta)]
    assert {path.name for path in folder.iterdir()} == {f"{STEM}.fa.gz", f"{STEM}.tsv",
                                                         f"{STEM}.manifest.json", "ATTRIBUTION.txt"}
    manifest_text = (folder / f"{STEM}.manifest.json").read_text()
    manifest = json.loads(manifest_text)
    assert manifest["origin"] == "zenodo" and manifest["zenodo"]["doi"] == DOI
    assert manifest["automatic_sources"] is True and manifest["index_included"] is False
    assert "/private/" not in manifest_text
    assert manifest["database"]["sha256"] == cli._sha256(fasta)
    assert calls[1].endswith("/SHA256SUMS")
    assert all(not call.endswith(".index.tar") for call in calls)
    out = capsys.readouterr().out
    assert DOI in out and "CC BY-SA 4.0" in out and "ATTRIBUTION.txt" in out
    assert "built on first use" in out and "bowtie2-build" in out
    cli.main(["--list", "--db-dir", str(tmp_path)])
    assert STEM in capsys.readouterr().out
    cli.main(["--verify", "--db-dir", str(tmp_path)])
    assert "[FAIL]" not in capsys.readouterr().out
    _, _, generate_index, headers, generate_headers = cli._database_artifacts(fasta, tmp_path, False)
    assert generate_index and headers == folder / f"{STEM}.tsv" and not generate_headers
    releases = {"bold": Release("BOLD_1", "d1", "u"), "unite": Release("UNITE_1", "d2", "u")}
    assert not is_outdated(fasta, releases, manifest["parameters"]["bold_retain"],
                           manifest["parameters"]["unite_retain"])


def test_fetch_version_and_record_environment(monkeypatch, tmp_path):
    calls = _record(monkeypatch)
    monkeypatch.setenv("BARCODEMAPPER_ZENODO_RECORD", "1")
    cli.main(["--fetch", "--fetch-version", VERSION, "--db-dir", str(tmp_path)])
    assert calls[:3] == ["https://zenodo.org/api/records/1/versions?size=100",
                         "https://zenodo.org/api/records/123",
                         "https://zenodo.org/api/files/123/SHA256SUMS"]


def test_fetched_database_supports_update_and_check(monkeypatch, tmp_path, capsys):
    _record(monkeypatch)
    _fetch(tmp_path)
    releases = {"bold": Release("BOLD_1", "d1", "u"), "unite": Release("UNITE_1", "d2", "u")}
    monkeypatch.setattr(cli, "discover_latest", lambda: releases)
    built = []
    monkeypatch.setattr(cli, "_build_latest", lambda args, parser, db_path, found=None: built.append(db_path))
    cli.main(["--update", "--db-dir", str(tmp_path)])
    assert not built and "up to date" in capsys.readouterr().out
    releases["bold"] = Release("BOLD_2", "d1", "u")
    with pytest.raises(SystemExit) as error:
        cli.main(["--update", "--check", "--db-dir", str(tmp_path)])
    assert error.value.code == 10 and not built
    cli.main(["--update", "--db-dir", str(tmp_path)])
    assert built == [tmp_path / STEM / f"{STEM}.fa.gz"]


def test_checksum_failure_leaves_no_install_or_stage(monkeypatch, tmp_path, capsys):
    _record(monkeypatch, bad_sha=f"{STEM}.fa.gz")
    with pytest.raises(SystemExit) as error:
        _fetch(tmp_path)
    assert error.value.code == 2
    assert "SHA-256 checksum mismatch" in capsys.readouterr().err
    assert list(tmp_path.iterdir()) == []


def test_existing_folder_needs_force_and_force_replaces_it(monkeypatch, tmp_path, capsys):
    _record(monkeypatch)
    _fetch(tmp_path)
    marker = tmp_path / STEM / "old.txt"
    marker.write_text("old")
    with pytest.raises(SystemExit) as error:
        _fetch(tmp_path)
    assert error.value.code == 2 and "--force" in capsys.readouterr().err
    assert marker.is_file()
    _fetch(tmp_path, "--force")
    assert not marker.exists()
    assert not list(tmp_path.glob(f".{STEM}.fetch-*"))
    assert not list(tmp_path.glob(f".{STEM}.previous-*"))


def test_force_rolls_back_failed_promotion(monkeypatch, tmp_path, capsys):
    _record(monkeypatch)
    _fetch(tmp_path)
    marker = tmp_path / STEM / "old.txt"
    marker.write_text("old")
    replace = Path.replace

    def fail_payload(self, target):
        if self.name == STEM and self.parent.name.startswith(f".{STEM}.fetch-") and target == tmp_path / STEM:
            raise OSError("promotion failed")
        return replace(self, target)

    monkeypatch.setattr(Path, "replace", fail_payload)
    with pytest.raises(SystemExit) as error:
        _fetch(tmp_path, "--force")
    assert error.value.code == 2 and "promotion failed" in capsys.readouterr().err
    assert marker.read_text() == "old"
    assert not list(tmp_path.glob(f".{STEM}.fetch-*"))
    assert not list(tmp_path.glob(f".{STEM}.previous-*"))


@pytest.mark.parametrize("extra", ["--build", "--update", "-d", "-n", "-b", "-u"])
def test_fetch_rejects_other_database_actions(monkeypatch, tmp_path, capsys, extra):
    monkeypatch.setattr(zenodo, "_open", lambda *a, **kw: pytest.fail("unexpected network"))
    command = ["--fetch", "--db-dir", str(tmp_path), extra]
    if extra in ("-d", "-n", "-b", "-u"):
        command.append("value")
    with pytest.raises(SystemExit) as error:
        cli.main(command)
    assert error.value.code == 2
    assert "--fetch" in capsys.readouterr().err


def test_fetch_index_is_current_and_required_by_verify(monkeypatch, tmp_path, capsys):
    _record(monkeypatch)
    _fetch(tmp_path, "--index")
    fasta = tmp_path / STEM / f"{STEM}.fa.gz"
    index, extension, generate, headers, generate_headers = cli._database_artifacts(fasta, tmp_path, False)
    assert index == cli._database_base(fasta) and extension == "bt2l" and not generate
    assert headers.is_file() and not generate_headers
    assert headers.stat().st_mtime_ns >= fasta.stat().st_mtime_ns
    assert all(Path(f"{index}.{part}.bt2l").stat().st_mtime_ns >= fasta.stat().st_mtime_ns
               for part in cli.INDEX_PARTS)
    cli.main(["--verify", "--db-dir", str(tmp_path)])
    assert "[FAIL]" not in capsys.readouterr().out
    Path(f"{index}.1.bt2l").unlink()
    with pytest.raises(SystemExit) as error:
        cli.main(["--verify", "--db-dir", str(tmp_path)])
    assert error.value.code == 1


@pytest.mark.parametrize("extra", [("../escape", tarfile.REGTYPE),
                                     (f"{STEM}.extra.bt2l", tarfile.REGTYPE),
                                     (f"{STEM}.1.bt2l", tarfile.SYMTYPE)])
def test_fetch_rejects_unsafe_index_archive(monkeypatch, tmp_path, capsys, extra):
    _record(monkeypatch, archive=_archive(extra))
    with pytest.raises(SystemExit) as error:
        _fetch(tmp_path, "--index")
    assert error.value.code == 2 and "unsafe or unexpected members" in capsys.readouterr().err
    assert list(tmp_path.iterdir()) == []


def test_missing_record_id_is_clear_without_network(monkeypatch, tmp_path, capsys):
    monkeypatch.delenv("BARCODEMAPPER_ZENODO_RECORD", raising=False)
    monkeypatch.setattr(zenodo, "DEFAULT_RECORD", None)
    monkeypatch.setattr(zenodo, "_open", lambda *a, **kw: pytest.fail("unexpected network"))
    with pytest.raises(SystemExit) as error:
        cli.main(["--fetch", "--db-dir", str(tmp_path)])
    assert error.value.code == 2
    assert "record ID is required" in capsys.readouterr().err


def test_fetch_rejects_other_record_licence(monkeypatch, tmp_path, capsys):
    _record(monkeypatch, edit=lambda payloads, metadata: metadata["metadata"].update(license="cc-by-4.0"))
    with pytest.raises(SystemExit):
        _fetch(tmp_path)
    assert "unexpected Zenodo database licence" in capsys.readouterr().err


@pytest.mark.parametrize("change,expected", [
    (lambda manifest: manifest["database"].update(sha256="0" * 64), "FASTA SHA-256"),
    (lambda manifest: manifest.update(automatic_sources=False), "automatic source provenance"),
])
def test_fetch_rejects_mismatched_or_manual_manifest(monkeypatch, tmp_path, capsys, change, expected):
    def edit(payloads, metadata):
        name = f"{STEM}.manifest.json"
        manifest = json.loads(payloads[name])
        change(manifest)
        payloads[name] = json.dumps(manifest).encode()

    _record(monkeypatch, edit=edit)
    with pytest.raises(SystemExit):
        _fetch(tmp_path)
    assert expected in capsys.readouterr().err
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("sums", [b"malformed\n", (b"0" * 64 + b"  ATTRIBUTION.txt\n") * 2])
def test_fetch_rejects_malformed_or_duplicate_sums(monkeypatch, tmp_path, capsys, sums):
    _record(monkeypatch, edit=lambda payloads, metadata: payloads.update(SHA256SUMS=sums))
    with pytest.raises(SystemExit):
        _fetch(tmp_path)
    assert "invalid SHA256SUMS" in capsys.readouterr().err


def test_fetch_rejects_missing_record_file(monkeypatch, tmp_path, capsys):
    _record(monkeypatch, edit=lambda payloads, metadata: metadata["files"].pop(0))
    with pytest.raises(SystemExit):
        _fetch(tmp_path)
    assert "missing required files" in capsys.readouterr().err


def test_fetch_accepts_labelled_version(monkeypatch, tmp_path):
    _record(monkeypatch, version="20261021_lab")
    _fetch(tmp_path, "--fetch-version", "20261021_lab")
    assert (tmp_path / "bmdb_20261021_lab" / "bmdb_20261021_lab.fa.gz").is_file()


def test_midway_download_failure_removes_stage_and_new_root(monkeypatch, tmp_path, capsys):
    _record(monkeypatch)
    original = zenodo._open

    class Broken(io.BytesIO):
        def read(self, size=-1):
            if self.tell():
                raise http.client.IncompleteRead(b"partial")
            return super().read(4)

    def opened(url, **kwargs):
        if url.endswith(f"/{STEM}.fa.gz"):
            return Broken(b"data")
        return original(url, **kwargs)

    monkeypatch.setattr(zenodo, "_open", opened)
    root = tmp_path / "new-db-dir"
    with pytest.raises(SystemExit):
        _fetch(root)
    assert not root.exists()
    assert "could not fetch" in capsys.readouterr().err


@pytest.mark.parametrize("option", ["--check", "--db-name", "-t", "-s", "-m", "-c"])
def test_fetch_rejects_non_fetch_options(monkeypatch, tmp_path, capsys, option):
    monkeypatch.setattr(zenodo, "_open", lambda *args, **kwargs: pytest.fail("unexpected network"))
    command = ["--fetch", "--db-dir", str(tmp_path), option]
    if option in ("--db-name", "-t", "-m", "-c"):
        command.append("value" if option not in ("-m", "-c") else "1")
    with pytest.raises(SystemExit) as error:
        cli.main(command)
    assert error.value.code == 2
    assert "--fetch" in capsys.readouterr().err


def test_force_index_rolls_back_failed_promotion(monkeypatch, tmp_path, capsys):
    _record(monkeypatch)
    _fetch(tmp_path, "--index")
    old = tmp_path / STEM / f"{STEM}.1.bt2l"
    old.write_bytes(b"old-index")
    replace = Path.replace

    def fail_payload(self, target):
        if self.name == STEM and self.parent.name.startswith(f".{STEM}.fetch-") and target == tmp_path / STEM:
            raise OSError("promotion failed")
        return replace(self, target)

    monkeypatch.setattr(Path, "replace", fail_payload)
    with pytest.raises(SystemExit):
        _fetch(tmp_path, "--force", "--index")
    assert old.read_bytes() == b"old-index"
    assert not list(tmp_path.glob(f".{STEM}.fetch-*"))
    assert "promotion failed" in capsys.readouterr().err


@pytest.mark.parametrize("files", [None, {}, ["bad"], [{"key": "x", "links": []}]])
def test_fetch_malformed_record_reports_one_line_error(monkeypatch, tmp_path, capsys, files):
    monkeypatch.setattr(zenodo, "_open", lambda *args, **kwargs: io.BytesIO(json.dumps({
        "id": 1, "doi": "doi", "metadata": {"version": VERSION, "license": "cc-by-sa-4.0"},
        "links": {"html": "https://zenodo.org/records/1"}, "files": files,
    }).encode()))
    with pytest.raises(SystemExit) as error:
        _fetch(tmp_path)
    assert error.value.code == 2
    assert "could not resolve Zenodo record" in capsys.readouterr().err


def test_fetch_checks_free_space_before_download(monkeypatch, tmp_path, capsys):
    _record(monkeypatch)
    monkeypatch.setattr(cli.shutil, "disk_usage", lambda path: type("Usage", (), {"free": 0})())
    with pytest.raises(SystemExit):
        _fetch(tmp_path)
    assert "insufficient free disk space" in capsys.readouterr().err
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("oversized", [True, False])
def test_fetch_rejects_index_member_size_or_total(monkeypatch, tmp_path, capsys, oversized):
    _record(monkeypatch)

    class FakeArchive:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def getmembers(self):
            members = [tarfile.TarInfo(f"{STEM}.{part}.bt2l") for part in cli.INDEX_PARTS]
            for member in members:
                member.size = zenodo.DATABASE_CAP + 1 if oversized else 2000
            return members

    monkeypatch.setattr(cli.tarfile, "open", lambda *args, **kwargs: FakeArchive())
    with pytest.raises(SystemExit):
        _fetch(tmp_path, "--index")
    assert "index archive exceeds declared size" in capsys.readouterr().err


def test_fetch_reports_coarse_progress(monkeypatch, tmp_path, capsys):
    _record(monkeypatch)
    download = zenodo.download_file

    def report(file, directory, expected_sha256=None, progress=None):
        if file.name.endswith(".fa.gz"):
            progress(100 * 1024 ** 2, 200 * 1024 ** 2)
        return download(file, directory, expected_sha256=expected_sha256, progress=progress)

    monkeypatch.setattr(zenodo, "download_file", report)
    _fetch(tmp_path)
    assert f"Downloading {STEM}.fa.gz: 50%" in capsys.readouterr().err
