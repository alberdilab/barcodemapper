"""Discover, describe, verify and prune bmdb_<YYYYMMDD> database folders."""

import json
import re
import shutil
from pathlib import Path

from barcodemapper.workflow.scripts.validate_database import validate_database

BMDB_DIR = re.compile(r"bmdb_(?P<date>\d{8})(?:_(?P<label>[A-Za-z0-9][A-Za-z0-9._-]*))?")
INDEX_PARTS = ("1", "2", "3", "4", "rev.1", "rev.2")
ANY = object()
FASTA_SUFFIXES = (".fa.gz", ".fa", ".fasta.gz", ".fasta")


def fasta_in(folder: Path) -> Path | None:
    """Return the FASTA inside a bmdb folder, if it has one."""
    for suffix in FASTA_SUFFIXES:
        candidate = folder / f"{folder.name}{suffix}"
        if candidate.is_file():
            return candidate
    return None


def database_stem(fasta: Path) -> Path:
    name = fasta.name
    for suffix in FASTA_SUFFIXES:
        if name.lower().endswith(suffix):
            return fasta.with_name(name[: -len(suffix)])
    return fasta.with_suffix("")


def find_databases(root: Path, label: object = ANY) -> list[tuple[str, Path]]:
    """Return (folder name, FASTA) pairs under root, oldest first.

    If root is itself a bmdb folder it is returned alone. A label restricts
    the result to databases built with that -n/--name (None: unlabelled only).
    """
    found = []
    folders = [root] if BMDB_DIR.fullmatch(root.name) else (
        sorted(root.iterdir()) if root.is_dir() else []
    )
    for folder in folders:
        match = BMDB_DIR.fullmatch(folder.name)
        if not match or not folder.is_dir() or (fasta := fasta_in(folder)) is None:
            continue
        if label is not ANY and match["label"] != label:
            continue
        found.append((folder.name, fasta))
    return sorted(found)


def select_database(root: Path, pin: str | None = None) -> Path | None:
    """Pick the newest database, or the one matching a folder name or date."""
    found = find_databases(root)
    if pin:
        wanted = pin if pin.startswith("bmdb_") else f"bmdb_{pin}"
        exact = [fasta for name, fasta in found if name == wanted]
        if exact:
            return exact[0]
        dated = [fasta for name, fasta in found if name.startswith(wanted + "_")]
        return dated[-1] if dated else None
    return found[-1][1] if found else None


def read_manifest(fasta: Path) -> dict | None:
    path = Path(f"{database_stem(fasta)}.manifest.json")
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return manifest if isinstance(manifest, dict) else None


def index_status(fasta: Path) -> str:
    """Return 'bt2l', 'bt2' or 'missing' for the Bowtie2 index beside the FASTA."""
    stem = database_stem(fasta)
    for extension in ("bt2l", "bt2"):
        if all((p := Path(f"{stem}.{part}.{extension}")).is_file() and p.stat().st_size > 0
               for part in INDEX_PARTS):
            return extension
    return "missing"


def provenance(fasta: Path) -> dict:
    """Summarise what a database was built from, using its manifest."""
    manifest = read_manifest(fasta) or {}
    sources = manifest.get("sources") if isinstance(manifest.get("sources"), dict) else {}
    params = manifest.get("parameters") if isinstance(manifest.get("parameters"), dict) else {}
    info = {"name": fasta.parent.name, "built_utc": manifest.get("created_utc"),
            "barcodemapper_version": manifest.get("barcodemapper_version")}
    for source in ("bold", "unite"):
        entry = sources.get(source) if isinstance(sources.get(source), dict) else {}
        info[source] = entry.get("release") or ("user-supplied" if entry.get("path") else None)
        info[f"{source}_doi"] = entry.get("doi")
        info[f"{source}_retain"] = params.get(f"{source}_retain")
    return info


def folder_size(folder: Path) -> int:
    return sum(path.stat().st_size for path in folder.iterdir() if path.is_file())


def _human(size: int) -> str:
    value = float(size)
    for unit in ("B", "KB", "MB", "GB"):
        if value < 1024 or unit == "GB":
            return f"{value:.0f} {unit}" if unit == "B" else f"{value:.1f} {unit}"
        value /= 1024
    return f"{size} B"


def describe(root: Path) -> list[str]:
    """Return one report block per database under root, newest first."""
    lines = []
    for name, fasta in reversed(find_databases(root)):
        info = provenance(fasta)
        lines += [
            name,
            f"  path:    {fasta}",
            f"  built:   {info['built_utc'] or 'unknown'} (BarcodeMapper {info['barcodemapper_version'] or 'unknown'})",
            f"  BOLD:    {info['bold'] or 'unknown'}   retain: {info['bold_retain'] or 'unknown'}",
            f"  UNITE:   {info['unite'] or 'unknown'}   retain: {info['unite_retain'] or 'unknown'}",
            f"  index:   {index_status(fasta)}",
            f"  size:    {_human(folder_size(fasta.parent))}",
            "",
        ]
    return lines


def verify(fasta: Path) -> list[tuple[str, bool, str]]:
    """Check that a database is complete and consistent; return (check, ok, detail)."""
    from barcodemapper.cli import _sha256  # local import avoids a cycle

    results = []
    try:
        counts = validate_database(fasta)
        results.append(("FASTA", True, ", ".join(f"{k}: {v}" for k, v in counts.items())))
    except (ValueError, OSError, EOFError) as exc:
        results.append(("FASTA", False, str(exc)))
    status = index_status(fasta)
    results.append(("Bowtie2 index", status != "missing", status))
    stem = database_stem(fasta)
    headers = Path(f"{stem}.tsv")
    results.append(("headers overview", headers.is_file() and headers.stat().st_size > 0, str(headers.name)))
    manifest = read_manifest(fasta)
    if manifest is None:
        results.append(("manifest", False, "missing or unreadable"))
    else:
        recorded = (manifest.get("database") or {}).get("sha256")
        actual = _sha256(fasta)
        results.append(("manifest checksum", recorded == actual,
                        "matches" if recorded == actual else "FASTA differs from the recorded SHA-256"))
    return results


def is_outdated(fasta: Path, releases: dict, bold_retain: str, unite_retain: str) -> bool:
    """True when the latest source releases or filters differ from what fasta was built from."""
    manifest = read_manifest(fasta)
    if not manifest or manifest.get("automatic_sources") is not True:
        return True
    sources = manifest.get("sources") or {}
    params = manifest.get("parameters") or {}
    if params.get("bold_retain") != bold_retain or params.get("unite_retain") != unite_retain:
        return True
    return any(
        (sources.get(name) or {}).get("release") != releases[name].identifier
        or (sources.get(name) or {}).get("doi") != releases[name].doi
        for name in ("bold", "unite")
    )


def prune(root: Path, keep: int, label: str | None) -> list[Path]:
    """Delete all but the newest `keep` databases (same label) under root; return removed folders."""
    if BMDB_DIR.fullmatch(root.name):
        return []
    found = [(name, fasta) for name, fasta in find_databases(root, label)
             if read_manifest(fasta) is not None]
    removed = []
    for _, fasta in found[:-keep] if keep else found:
        shutil.rmtree(fasta.parent)
        removed.append(fasta.parent)
    return removed
