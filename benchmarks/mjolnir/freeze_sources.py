"""Acquire one immutable BOLD/UNITE source snapshot for the benchmark plans."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import tempfile
from datetime import datetime, timezone
from pathlib import Path

from barcodemapper.__version__ import __version__
from barcodemapper.reference_sources import discover_latest, download_sources


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    destination = args.destination.expanduser().resolve()
    if destination.exists():
        parser.error(f"destination already exists: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix="sources-staging-", dir=destination.parent))
    try:
        releases = discover_latest()
        paths, provenance = download_sources(releases, stage)
        record = {
            "downloaded_utc": datetime.now(timezone.utc).isoformat(),
            "barcodemapper_version": __version__,
            "sources": provenance,
        }
        (stage / "source-releases.json").write_text(
            json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        checksums = "".join(
            f"{provenance[name]['sha256']}  {paths[name].name}\n"
            for name in ("bold", "unite")
        )
        (stage / "SHA256SUMS").write_text(checksums, encoding="utf-8")
        os.replace(stage, destination)
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise


if __name__ == "__main__":
    main()
