"""Stage a built BarcodeMapper database and upload it to a Zenodo draft."""

import argparse
import getpass
import hashlib
import json
import os
import re
import shutil
import sys
import tempfile
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlsplit
from urllib.request import HTTPRedirectHandler, Request, build_opener


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from barcodemapper.databases import BMDB_DIR, INDEX_PARTS  # noqa: E402
from barcodemapper.publish_meta import (attribution_text, portable_manifest,
                                        validate_source_manifest, zenodo_metadata)  # noqa: E402


DEFAULT_CREATORS = [
    {"name": "Alberdi, Antton", "affiliation": "University of Copenhagen"},
    {"name": "Aizpurua, Ostaizka", "affiliation": "University of Copenhagen"},
]
SOURCE_LICENCES = {"bold": "CC BY-SA 4.0", "unite": "CC BY 4.0"}
BLOCK_SIZE = 1024 * 1024


def _digest(path: Path, algorithm: str) -> str:
    checksum = hashlib.new(algorithm)
    with path.open("rb") as source:
        for block in iter(lambda: source.read(BLOCK_SIZE), b""):
            checksum.update(block)
    return checksum.hexdigest()


def _preflight(folder: Path) -> dict:
    if not folder.is_dir() or not BMDB_DIR.fullmatch(folder.name):
        raise ValueError("--database must be an existing bmdb_<YYYYMMDD>[_label] folder")
    stem = folder.name
    if ".." in stem:
        raise ValueError("database folder name is unsafe for Zenodo filenames")
    fasta = folder / f"{stem}.fa.gz"
    overview = folder / f"{stem}.tsv"
    manifest_path = folder / f"{stem}.manifest.json"
    for path in (fasta, overview, manifest_path):
        if not path.is_file():
            raise ValueError(f"database file is missing: {path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    validate_source_manifest(manifest, _digest(fasta, "sha256"))
    sources = manifest["sources"]
    for name, expected_licence in SOURCE_LICENCES.items():
        source = sources.get(name)
        if not source.get("sha256"):
            raise ValueError(f"manifest lacks {name.upper()} source SHA-256")
        recorded_licence = source.get("licence", source.get("license"))
        if recorded_licence and recorded_licence != expected_licence:
            raise ValueError(f"{name.upper()} source licence differs from expected {expected_licence}")
    bold = sources["bold"]
    if bold.get("doi") == "https://doi.org/10.5883/DP-Latest" and not bold.get("release"):
        raise ValueError("BOLD DP-Latest DOI requires a recorded snapshot release ID")
    return manifest


def _creators(values: list[str] | None) -> list[dict]:
    if not values:
        return DEFAULT_CREATORS
    creators = []
    for value in values:
        fields = value.split(";")
        if len(fields) != 3 or not all(fields) or "," not in fields[0]:
            raise ValueError("--creator must be 'Family, Given;Affiliation;ORCID'")
        name, affiliation, orcid = (field.strip() for field in fields)
        if not all((name, affiliation, orcid)):
            raise ValueError("--creator has an empty field")
        creators.append({"name": name, "affiliation": affiliation, "orcid": orcid})
    return creators


def _metadata(manifest: dict, version: str, creators: list[dict]) -> dict:
    return zenodo_metadata(manifest, version=version, creators=creators)


def _stage(folder: Path, stage: Path, manifest: dict, include_index: bool) -> list[Path]:
    stem = folder.name
    files = []
    for suffix in (".fa.gz", ".tsv"):
        name = f"{stem}{suffix}"
        with (folder / name).open("rb") as source, (stage / name).open("wb") as target:
            shutil.copyfileobj(source, target, BLOCK_SIZE)
        files.append(stage / name)
    if _digest(files[0], "sha256") != manifest["database"]["sha256"]:
        raise ValueError("staged FASTA SHA-256 does not match the manifest")
    portable = portable_manifest(manifest)
    portable["database"]["path"] = f"{stem}.fa.gz"
    serialized = json.dumps(portable, indent=2) + "\n"
    if getpass.getuser() in serialized or str(Path.home()) in serialized:
        raise ValueError("portable manifest contains a local user name or home path")
    manifest_file = stage / f"{stem}.manifest.json"
    manifest_file.write_text(serialized, encoding="utf-8")
    files.append(manifest_file)
    attribution = attribution_text(manifest)
    attribution += ("\nExpected source licences for redistribution:\n"
                    "  BOLD: CC BY-SA 4.0\n  UNITE: CC BY 4.0\n")
    attribution_file = stage / "ATTRIBUTION.txt"
    attribution_file.write_text(attribution, encoding="utf-8")
    files.append(attribution_file)
    if include_index:
        for part in INDEX_PARTS:
            path = folder / f"{stem}.{part}.bt2l"
            if not path.is_file() or path.stat().st_size == 0:
                raise ValueError(f"missing or empty Bowtie2 index file: {path}")
            files.append(path)
    sums = stage / "SHA256SUMS"
    sums.write_text("".join(f"{_digest(path, 'sha256')}  {path.name}\n" for path in files),
                    encoding="utf-8")
    files.append(sums)
    return files


def _safe_url(url: str, host: str) -> str:
    parts = urlsplit(url)
    if (parts.scheme != "https" or parts.netloc != host or parts.username or parts.password or
            parts.query or parts.fragment or not parts.path.startswith("/api/")):
        raise ValueError("Zenodo API link has an unexpected host or path")
    return url


class _SameHostRedirect(HTTPRedirectHandler):
    def __init__(self, host: str):
        self.host = host

    def redirect_request(self, request, response, code, message, headers, new_url):
        if request.get_method() != "GET":
            raise ValueError("Zenodo API redirected a modifying request")
        _safe_url(new_url, self.host)
        return super().redirect_request(request, response, code, message, headers, new_url)


def _request(method: str, url: str, token: str, host: str, *, payload=None,
             file_path: Path | None = None):
    """Make the one authenticated HTTP call; reject cross-host links and redirects."""
    _safe_url(url, host)
    headers = {"Authorization": f"Bearer {token}"}
    data = None
    if payload is not None:
        data = json.dumps(payload).encode("utf-8")
        headers["Content-Type"] = "application/json"
    if file_path is not None:
        headers["Content-Type"] = "application/octet-stream"
        headers["Content-Length"] = str(file_path.stat().st_size)
        data = _file_blocks(file_path)
    request = Request(url, data=data, headers=headers, method=method)
    try:
        with build_opener(_SameHostRedirect(host)).open(request, timeout=3600) as response:
            body = response.read()
    except HTTPError as exc:
        raise RuntimeError(f"Zenodo API returned HTTP {exc.code}") from None
    except (OSError, URLError) as exc:
        reason = getattr(exc, "reason", None) or exc
        raise RuntimeError(f"Zenodo API request failed ({type(exc).__name__}: {reason})") from None
    return json.loads(body) if body else None


def _file_blocks(path: Path):
    with path.open("rb") as source:
        yield from iter(lambda: source.read(BLOCK_SIZE), b"")


def _record_id(value) -> str:
    if not re.fullmatch(r"[0-9]+", str(value)):
        raise ValueError("Zenodo returned an invalid deposition ID")
    return str(value)


def _file_id(item: dict) -> str:
    file_id = str(item["id"])
    if not re.fullmatch(r"[A-Za-z0-9-]+", file_id):
        raise ValueError("Zenodo returned an invalid file ID")
    return file_id


def _draft_file_name(item: dict) -> str:
    return item.get("filename", item.get("key"))


def _upload(files: list[Path], metadata: dict, record: str | None, draft_id: str | None, publish: bool,
            token: str, host: str) -> None:
    base = f"https://{host}"
    if record:
        latest = _request("GET", f"{base}/api/records/{record}", token, host)
        latest_id = _record_id(latest["id"])
        old = _request("POST", f"{base}/api/deposit/depositions/{latest_id}/actions/newversion", token, host)
        draft_url = _safe_url(old["links"]["latest_draft"], host)
        draft = _request("GET", draft_url, token, host)
    elif draft_id:
        draft = _request("GET", f"{base}/api/deposit/depositions/{draft_id}", token, host)
        if draft.get("state") != "unsubmitted":
            raise ValueError(f"Zenodo deposition {draft_id} is not an unpublished draft")
    else:
        draft = _request("POST", f"{base}/api/deposit/depositions", token, host, payload={})
    draft_id = _record_id(draft["id"])
    draft_base = f"{base}/api/deposit/depositions/{draft_id}"
    print(f"Draft: {base}/deposit/{draft_id}")
    if record:
        for old_file in draft.get("files", []):
            _request("DELETE", f"{draft_base}/files/{_file_id(old_file)}", token, host)
    bucket = _safe_url(draft["links"]["bucket"], host).rstrip("/")
    existing = {_draft_file_name(item): item for item in draft.get("files", [])} if draft_id else {}
    for name in sorted(existing.keys() - {path.name for path in files}):
        print(f"Existing draft file outside this upload plan: {name}")
    for path in files:
        md5 = _digest(path, "md5")
        old_file = existing.get(path.name)
        if old_file is not None:
            old_md5 = str(old_file.get("checksum", "")).removeprefix("md5:")
            if old_md5.lower() == md5:
                print(f"Skipping {path.name}: matching MD5 in draft")
                continue
            _request("DELETE", f"{draft_base}/files/{_file_id(old_file)}", token, host)
        for attempt in range(1, 4):
            print(f"Uploading {path.name}: {path.stat().st_size} bytes (attempt {attempt}/3)")
            try:
                uploaded = _request("PUT", f"{bucket}/{path.name}", token, host, file_path=path)
            except (RuntimeError, OSError) as exc:
                if attempt == 3:
                    raise RuntimeError(f"upload of {path.name} failed after 3 attempts: {exc}") from exc
                current = _request("GET", draft_base, token, host)
                for item in current.get("files", []):
                    if _draft_file_name(item) == path.name:
                        _request("DELETE", f"{draft_base}/files/{_file_id(item)}", token, host)
                time.sleep(0.2 * attempt)
                continue
            checksum = uploaded.get("checksum") if isinstance(uploaded, dict) else None
            if checksum != f"md5:{md5}":
                raise ValueError(f"Zenodo MD5 checksum mismatch for {path.name}")
            break
    _request("PUT", draft_base, token, host, payload={"metadata": metadata})
    if publish:
        published = _request("POST", f"{draft_base}/actions/publish", token, host)
        print(f"DOI: {published['doi']}".replace(token, "[REDACTED]"))
        print(f"Concept DOI: {published['conceptdoi']}".replace(token, "[REDACTED]"))
        print(f"Set barcodemapper.zenodo.DEFAULT_RECORD = {published['conceptrecid']}")


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", required=True, type=Path)
    target = parser.add_mutually_exclusive_group()
    target.add_argument("--record", help="Latest published deposition ID for a new version")
    target.add_argument("--draft", help="Existing unpublished deposition ID to resume")
    parser.add_argument("--include-index", action="store_true")
    parser.add_argument("--sandbox", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--publish", action="store_true")
    parser.add_argument("--creator", action="append")
    parser.add_argument("--staging", type=Path)
    args = parser.parse_args(argv)
    token = os.environ.get("ZENODO_TOKEN", "")
    try:
        if args.record:
            _record_id(args.record)
        if args.draft:
            _record_id(args.draft)
        folder = args.database.expanduser().resolve()
        manifest = _preflight(folder)
        metadata = _metadata(manifest, folder.name.removeprefix("bmdb_"), _creators(args.creator))
        if not args.dry_run and not token:
            raise ValueError("ZENODO_TOKEN is required unless --dry-run is used")
        stage_parent = args.staging.expanduser().resolve() if args.staging else folder.parent
        if not stage_parent.is_dir():
            raise ValueError(f"staging parent directory does not exist: {stage_parent}")
        stage_name = tempfile.mkdtemp(prefix=f".{folder.name}.publish-", dir=stage_parent)
        try:
            files = _stage(folder, Path(stage_name), manifest, args.include_index)
            print("Planned files:")
            for path in files:
                print(f"  {path.name}: {path.stat().st_size} bytes")
            print("Expected source licences: BOLD CC BY-SA 4.0; UNITE CC BY 4.0")
            print("Metadata:")
            metadata_text = json.dumps(metadata, indent=2)
            print(metadata_text.replace(token, "[REDACTED]") if token else metadata_text)
            if not args.dry_run:
                host = "sandbox.zenodo.org" if args.sandbox else "zenodo.org"
                _upload(files, metadata, args.record, args.draft, args.publish, token, host)
        finally:
            # Never let a cleanup failure (e.g. NFS) hide the real error.
            shutil.rmtree(stage_name, ignore_errors=True)
    except (OSError, ValueError, KeyError, TypeError, RuntimeError) as exc:
        # Even a mocked HTTP helper may include credentials in its exception text.
        message = str(exc).replace(token, "[REDACTED]") if token else str(exc)
        parser.error(f"{message}; any created Zenodo draft was left in place")


if __name__ == "__main__":
    main()
