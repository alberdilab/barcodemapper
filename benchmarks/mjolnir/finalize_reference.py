"""Hash and freeze a complete BarcodeMapper reference snapshot."""

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


PARTS = ("1", "2", "3", "4", "rev.1", "rev.2")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta", type=Path)
    args = parser.parse_args()
    fasta = args.fasta.expanduser().resolve()
    if not fasta.name.endswith(".fa.gz"):
        parser.error("Expected a .fa.gz database")
    base = fasta.with_suffix("").with_suffix("")
    artifacts = [fasta, Path(str(base) + ".tsv"),
                 Path(str(base) + ".manifest.json")]
    artifacts += [Path(f"{base}.{part}.bt2l") for part in PARTS]
    for path in artifacts:
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing or empty reference artifact: {path}")
    manifest = json.loads(artifacts[2].read_text(encoding="utf-8"))
    fasta_sha = sha256(fasta)
    if manifest["database"]["sha256"] != fasta_sha:
        raise ValueError("BarcodeMapper database manifest checksum mismatch")
    cleaned_manifest = fasta.parent / "cleaned-sources" / "empty-record-filter.json"
    source_filter = None
    if cleaned_manifest.is_file():
        cleaned = json.loads(cleaned_manifest.read_text(encoding="utf-8"))
        for name in ("bold", "unite"):
            source = manifest["sources"][name]
            cleaned_source = cleaned["sources"][name]
            if (source["path"] != cleaned_source["cleaned_path"] or
                    source["sha256"] != cleaned_source["cleaned_sha256"]):
                raise ValueError(f"Reference source does not match cleaned source: {name}")
        source_filter = {"path": str(cleaned_manifest), "sha256": sha256(cleaned_manifest)}
    checksums = Path(str(base) + ".SHA256SUMS")
    if checksums.exists():
        parser.error("Snapshot checksum file already exists")
    records = [{"path": str(path), "bytes": path.stat().st_size,
                "sha256": sha256(path)} for path in artifacts]
    checksums.write_text("".join(f"{item['sha256']}  {Path(item['path']).name}\n"
                                 for item in records), encoding="utf-8")
    record = Path(str(base) + ".snapshot.json")
    record.write_text(json.dumps({"frozen_utc": datetime.now(timezone.utc).isoformat(),
                                  "source_filter": source_filter,
                                  "artifacts": records,
                                  "checksums_sha256": sha256(checksums)},
                                 indent=2) + "\n", encoding="utf-8")
    for path in artifacts + [checksums, record]:
        path.chmod(path.stat().st_mode & ~0o222)
    print(f"Frozen {len(artifacts)} reference artifacts: {record}")


if __name__ == "__main__":
    main()
