"""Offline checks for Zenodo record resolution and verified downloads."""

import hashlib
import io
import json
import os

import pytest

from barcodemapper import zenodo
from urllib.request import Request


class Response(io.BytesIO):
    pass


CONTENT = b"reference database\n"
MD5 = hashlib.md5(CONTENT).hexdigest()
SHA256 = hashlib.sha256(CONTENT).hexdigest()


def record_data(record_id=123, version="20261021", *, base=zenodo.ZENODO_URL, legacy=False):
    item = {
        "checksum": f"md5:{MD5}",
        "links": {"download" if legacy else "self": f"{base}/api/files/{record_id}/database.tar.gz"},
    }
    item.update({"filename": "database.tar.gz", "filesize": len(CONTENT)} if legacy else
                {"key": "database.tar.gz", "size": len(CONTENT)})
    return {
        "id": record_id,
        "doi": f"10.5281/zenodo.{record_id}",
        "conceptdoi": "10.5281/zenodo.1",
        "links": {"html": f"{base}/records/{record_id}"},
        "metadata": {"version": version, "license": "cc-by-sa-4.0"},
        "files": [item],
    }


def test_latest_concept_record_and_download(monkeypatch, tmp_path):
    calls = []

    def open_zenodo(url, **kwargs):
        calls.append(url)
        if url.endswith("/api/records/1"):
            return Response(json.dumps(record_data()).encode())
        assert url == f"{zenodo.ZENODO_URL}/api/files/123/database.tar.gz"
        return Response(CONTENT)

    monkeypatch.setattr(zenodo, "_open", open_zenodo)
    found = zenodo.resolve_record(1)
    assert found.record_id == 123
    assert found.doi == "10.5281/zenodo.123"
    assert found.version == "20261021"
    assert found.url == "https://zenodo.org/records/123"
    assert found.license == "cc-by-sa-4.0"
    assert found.files[0].name == "database.tar.gz"
    assert found.files[0].size == len(CONTENT)
    assert found.files[0].md5 == MD5
    progress = []
    path = zenodo.download_file(found.files[0], tmp_path, expected_sha256=SHA256,
                                progress=lambda current, total: progress.append((current, total)))
    assert path.read_bytes() == CONTENT
    assert progress == [(len(CONTENT), len(CONTENT))]
    assert calls == ["https://zenodo.org/api/records/1", found.files[0].url]


def test_version_pin_uses_versions_endpoint_and_legacy_file_keys(monkeypatch):
    calls = []

    def open_zenodo(url, **kwargs):
        calls.append(url)
        if url.endswith("/versions?size=100"):
            return Response(json.dumps({"hits": {"hits": [
                {"id": 122, "metadata": {"version": "20261020"}},
                {"id": 123, "metadata": {"version": "20261021"}},
            ]}}).encode())
        assert url.endswith("/api/records/122")
        return Response(json.dumps(record_data(122, "20261020", legacy=True)).encode())

    monkeypatch.setattr(zenodo, "_open", open_zenodo)
    found = zenodo.resolve_record(1, version="20261020")
    assert found.record_id == 122
    assert found.files[0].name == "database.tar.gz"
    assert found.files[0].size == len(CONTENT)
    assert found.files[0].url.endswith("/api/files/122/database.tar.gz")
    assert calls == ["https://zenodo.org/api/records/1/versions?size=100", "https://zenodo.org/api/records/122"]


def test_current_zenodo_links_use_self_html_and_license_object(monkeypatch):
    data = record_data(5)
    data["links"] = {"self_html": f"{zenodo.ZENODO_URL}/records/5"}
    data["metadata"]["license"] = {"id": "cc-by-sa-4.0"}
    monkeypatch.setattr(zenodo, "_open", lambda url, **kwargs: Response(json.dumps(data).encode()))
    record = zenodo.resolve_record(5)
    assert record.url == f"{zenodo.ZENODO_URL}/records/5" and record.license == "cc-by-sa-4.0"


def test_unknown_version_is_clear_and_does_not_fetch_a_record(monkeypatch):
    calls = []

    def open_zenodo(url, **kwargs):
        calls.append(url)
        return Response(b'{"hits": {"hits": []}}')

    monkeypatch.setattr(zenodo, "_open", open_zenodo)
    with pytest.raises(ValueError, match="version 'missing' was not found"):
        zenodo.resolve_record(1, version="missing")
    assert calls == ["https://zenodo.org/api/records/1/versions?size=100"]


@pytest.mark.parametrize("expected_sha256,md5,match", [
    (None, "0" * 32, "MD5 checksum mismatch"),
    ("0" * 64, MD5, "SHA-256 checksum mismatch"),
])
def test_checksum_mismatch_removes_partial_file(monkeypatch, tmp_path, expected_sha256, md5, match):
    monkeypatch.setattr(zenodo, "_open", lambda *args, **kwargs: Response(CONTENT))
    file = zenodo.ZenodoFile("database.tar.gz", len(CONTENT), md5,
                             "https://zenodo.org/api/files/database.tar.gz")
    with pytest.raises(ValueError, match=match):
        zenodo.download_file(file, tmp_path, expected_sha256=expected_sha256)
    assert list(tmp_path.iterdir()) == []


def test_non_https_base_and_download_urls_are_rejected(monkeypatch, tmp_path):
    monkeypatch.setattr(zenodo, "_open", lambda *args, **kwargs: pytest.fail("Unexpected network access"))
    with pytest.raises(ValueError, match="HTTPS"):
        zenodo.resolve_record(1, base="http://zenodo.org")
    file = zenodo.ZenodoFile("database.tar.gz", len(CONTENT), MD5,
                             "http://zenodo.org/api/files/database.tar.gz")
    with pytest.raises(ValueError, match="HTTPS"):
        zenodo.download_file(file, tmp_path)


def test_cross_host_download_url_is_rejected(monkeypatch):
    data = record_data()
    data["files"][0]["links"]["self"] = "https://other.example/database.tar.gz"
    monkeypatch.setattr(zenodo, "_open", lambda *args, **kwargs: Response(json.dumps(data).encode()))
    with pytest.raises(ValueError, match="base host"):
        zenodo.resolve_record(1)


@pytest.mark.parametrize("name", ["", "../database", "sub/database", "sub\\database", ".hidden", "a..b", "a\x00b"])
def test_unsafe_file_names_are_rejected(monkeypatch, name):
    data = record_data()
    data["files"][0]["key"] = name
    monkeypatch.setattr(zenodo, "_open", lambda *args, **kwargs: Response(json.dumps(data).encode()))
    with pytest.raises(ValueError, match="Unsafe Zenodo file name"):
        zenodo.resolve_record(1)


def test_missing_record_id_has_clear_error_without_network(monkeypatch):
    monkeypatch.setattr(zenodo, "_open", lambda *args, **kwargs: pytest.fail("Unexpected network access"))
    with pytest.raises(ValueError, match="record ID is required"):
        zenodo.resolve_record()


def test_sandbox_base_from_environment(monkeypatch, tmp_path):
    base = "https://sandbox.zenodo.org"
    monkeypatch.setenv("BARCODEMAPPER_ZENODO_URL", base)
    calls = []

    def open_zenodo(url, **kwargs):
        calls.append(url)
        if url.endswith("/api/records/1"):
            return Response(json.dumps(record_data(base=base)).encode())
        return Response(CONTENT)

    monkeypatch.setattr(zenodo, "_open", open_zenodo)
    found = zenodo.resolve_record(1)
    assert zenodo.download_file(found.files[0], tmp_path).read_bytes() == CONTENT
    assert calls == [f"{base}/api/records/1", found.files[0].url]


def test_explicit_base_overrides_environment(monkeypatch):
    monkeypatch.setenv("BARCODEMAPPER_ZENODO_URL", "https://sandbox.zenodo.org")
    monkeypatch.setattr(zenodo, "_open", lambda url, **kwargs:
                        Response(json.dumps(record_data(base="https://zenodo.org")).encode())
                        if url == "https://zenodo.org/api/records/1" else pytest.fail(url))
    assert zenodo.resolve_record(1, base="https://zenodo.org").record_id == 123


def test_redirect_handler_rejects_cross_host_and_non_https():
    handler = zenodo._SameHostRedirect("zenodo.org")
    request = Request("https://zenodo.org/api/records/1")
    assert handler.redirect_request(request, None, 302, "Found", {},
                                    "https://zenodo.org/api/records/2").full_url.endswith("/2")
    for url in ("https://other.example/file", "http://zenodo.org/file", "file:///tmp/file"):
        with pytest.raises(ValueError):
            handler.redirect_request(request, None, 302, "Found", {}, url)


def test_nonstandard_base_warns(monkeypatch, capsys):
    base = "https://mirror.example"
    monkeypatch.setenv("BARCODEMAPPER_ZENODO_URL", base)
    monkeypatch.setattr(zenodo, "_open", lambda *args, **kwargs:
                        Response(json.dumps(record_data(base=base)).encode()))
    zenodo.resolve_record(1)
    assert "Warning: using nonstandard Zenodo host mirror.example" in capsys.readouterr().err


def test_version_pin_follows_same_host_pages(monkeypatch):
    calls = []

    def opened(url, **kwargs):
        calls.append(url)
        if url.endswith("versions?size=100"):
            return Response(json.dumps({"hits": {"hits": []}, "links": {
                "next": "https://zenodo.org/api/records/1/versions?page=2"}}).encode())
        if url.endswith("versions?page=2"):
            return Response(json.dumps({"hits": {"hits": [{"id": 123, "metadata": {
                "version": "20261021_lab"}}]}}).encode())
        return Response(json.dumps(record_data(version="20261021_lab")).encode())

    monkeypatch.setattr(zenodo, "_open", opened)
    assert zenodo.resolve_record(1, version="20261021_lab").version == "20261021_lab"
    assert len(calls) == 3


@pytest.mark.parametrize("bad", [None, {}, ["not a file"], [{"links": []}],
                                      [{"key": "a", "size": "1", "links": {}}]])
def test_record_rejects_malformed_files(monkeypatch, bad):
    data = record_data()
    data["files"] = bad
    monkeypatch.setattr(zenodo, "_open", lambda *args, **kwargs: Response(json.dumps(data).encode()))
    with pytest.raises((ValueError, TypeError, AttributeError)):
        zenodo.resolve_record(1)


@pytest.mark.parametrize("size", [len(CONTENT) - 1, len(CONTENT) + 1, -1, 1.5,
                                       zenodo.DATABASE_CAP + 1])
def test_download_enforces_declared_size(monkeypatch, tmp_path, size):
    monkeypatch.setattr(zenodo, "_open", lambda *args, **kwargs: Response(CONTENT))
    file = zenodo.ZenodoFile("database.fa.gz", size, MD5, "https://zenodo.org/file")
    with pytest.raises(ValueError):
        zenodo.download_file(file, tmp_path)
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("name,size", [("SHA256SUMS", zenodo.METADATA_CAP + 1),
                                          ("ATTRIBUTION.txt", zenodo.METADATA_CAP + 1),
                                          ("database.fa.gz", zenodo.DATABASE_CAP + 1)])
def test_record_rejects_file_above_cap(monkeypatch, name, size):
    data = record_data()
    data["files"][0].update(key=name, size=size)
    monkeypatch.setattr(zenodo, "_open", lambda *args, **kwargs: Response(json.dumps(data).encode()))
    with pytest.raises(ValueError, match="file size"):
        zenodo.resolve_record(1)


def test_downloaded_file_follows_umask_not_private_temp_mode(monkeypatch, tmp_path):
    monkeypatch.setattr(zenodo, "_open", lambda url, **kwargs: Response(CONTENT))
    record = zenodo.ZenodoFile("database.tar.gz", len(CONTENT), hashlib.md5(CONTENT).hexdigest(),
                               f"{zenodo.ZENODO_URL}/api/files/1/database.tar.gz", "zenodo.org")
    old = os.umask(0o022)
    try:
        path = zenodo.download_file(record, tmp_path)
    finally:
        os.umask(old)
    assert path.stat().st_mode & 0o777 == 0o644
