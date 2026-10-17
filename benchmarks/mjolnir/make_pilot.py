"""Create a deterministic, paired 100,000-fragment pilot from a verified run."""

import argparse
import gzip
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def checksum(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_record(handle):
    lines = [handle.readline() for _ in range(4)]
    if not lines[0]:
        return None
    if not all(lines) or not lines[0].startswith(b"@") or not lines[2].startswith(b"+"):
        raise ValueError("Malformed or truncated FASTQ record")
    return lines


def read_id(header):
    name = header.split(None, 1)[0][1:]
    if name.endswith((b"/1", b"/2")):
        name = name[:-2]
    return name


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs_dir", type=Path)
    parser.add_argument("run_accession")
    parser.add_argument("--pairs", type=int, default=100000)
    args = parser.parse_args()
    if args.pairs < 1:
        parser.error("--pairs must be positive")
    inputs = args.inputs_dir.expanduser().resolve()
    raw_manifest = inputs / "raw" / f"{args.run_accession}.manifest.json"
    if not raw_manifest.is_file():
        parser.error("Verified raw FASTQ manifest is missing")
    raw = json.loads(raw_manifest.read_text(encoding="utf-8"))
    if len(raw["files"]) != 2:
        parser.error("Raw manifest must contain exactly two mates")
    source_paths = [Path(item["path"]) for item in raw["files"]]
    pilot_dir = inputs / "pilot"
    pilot_dir.mkdir(parents=True, exist_ok=True)
    targets = [pilot_dir / f"{args.run_accession}_{mate}.fastq.gz" for mate in (1, 2)]
    manifest = pilot_dir / f"{args.run_accession}.manifest.json"
    if manifest.exists() or any(path.exists() for path in targets):
        parser.error("Pilot output already exists; preserve the frozen sample")
    partials = [pilot_dir / (path.name + ".part") for path in targets]
    with gzip.open(str(source_paths[0]), "rb") as r1, gzip.open(str(source_paths[1]), "rb") as r2:
        with partials[0].open("wb") as f1, partials[1].open("wb") as f2:
            with gzip.GzipFile(fileobj=f1, filename="", mode="wb", mtime=0) as w1, gzip.GzipFile(fileobj=f2, filename="", mode="wb", mtime=0) as w2:
                for index in range(args.pairs):
                    mate1 = read_record(r1)
                    mate2 = read_record(r2)
                    if mate1 is None or mate2 is None:
                        raise ValueError(f"Raw run ended after {index} pairs")
                    if read_id(mate1[0]) != read_id(mate2[0]):
                        raise ValueError(f"Unpaired read IDs at pair {index + 1}")
                    w1.writelines(mate1)
                    w2.writelines(mate2)
    for partial, target in zip(partials, targets):
        partial.rename(target)
    manifest.write_text(json.dumps({
        "run_accession": args.run_accession,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "selection": f"first {args.pairs} synchronized pairs from ENA FASTQs",
        "pairs": args.pairs,
        "raw_manifest_sha256": checksum(raw_manifest),
        "mates": [{"path": str(path), "bytes": path.stat().st_size,
                   "sha256": checksum(path)} for path in targets],
    }, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {args.pairs} verified pairs to {pilot_dir}")


if __name__ == "__main__":
    main()
