"""Index eligible frozen source accessions for taxon-level holdout design."""

import argparse
import hashlib
import json
import os
import sqlite3
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


RANKS = "kpcofgs"


def fasta_records(path):
    with path.open("rb") as handle:
        header = None
        offset = None
        chunks = []
        while True:
            position = handle.tell()
            line = handle.readline()
            if not line:
                if header is not None:
                    yield offset, header, b"".join(chunks)
                return
            if line.startswith(b">"):
                if header is not None:
                    yield offset, header, b"".join(chunks)
                offset = position
                header = line[1:].decode("utf-8").strip()
                chunks = []
            elif header is not None:
                chunks.append(line.strip())


def bold_record(header):
    parts = header.split("|")
    if len(parts) < 4 or parts[1] != "COI-5P":
        return None
    names = parts[3].split(",")
    if len(names) > 7:
        del names[5]
    names = (names + [""] * 7)[:7]
    if names[0] != "Animalia":
        return None
    names = ["" if name.lower() == "none" else name.replace(" ", "_")
             for name in names]
    taxonomy = ";".join(f"{rank}__{name}" for rank, name in zip(RANKS, names))
    return parts[0], "COI-5P", taxonomy, names[4], names[5], names[6]


def unite_record(header):
    parts = header.split("|")
    if len(parts) < 5 or not parts[1]:
        return None
    levels = parts[-1].split(";")
    names = {}
    incertae = {
        "p": "_phy_Incertae_sedis",
        "c": "_cls_Incertae_sedis",
        "o": "_ord_Incertae_sedis",
        "f": "_fam_Incertae_sedis",
        "g": "_gen_Incertae_sedis",
    }
    for level in levels:
        if "__" not in level:
            continue
        rank, name = level.split("__", 1)
        if rank not in RANKS:
            continue
        if (rank == "s" and name.endswith("_sp")) or (rank in incertae and incertae[rank] in name):
            name = ""
        names[rank] = name
    if names.get("k") != "Fungi":
        return None
    taxonomy = ";".join(f"{rank}__{names.get(rank, '')}" for rank in RANKS)
    return parts[1], "ITS", taxonomy, names.get("f", ""), names.get("g", ""), names.get("s", "")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sources_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    sources = args.sources_dir.expanduser().resolve()
    output = args.output_dir.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    database = output / "source-index.sqlite"
    stage = output / "source-index.sqlite.partial"
    summary = output / "source-index.json"
    if database.exists() or summary.exists():
        parser.error("Source index already exists")
    if stage.exists():
        stage.unlink()
    releases = sources / "source-releases.json"
    provenance = json.loads(releases.read_text(encoding="utf-8"))
    conn = sqlite3.connect(str(stage))
    try:
        conn.execute("PRAGMA journal_mode=DELETE")
        conn.execute("PRAGMA synchronous=NORMAL")
        conn.execute("CREATE TABLE records (source TEXT, source_offset INTEGER, accession TEXT, marker TEXT, taxonomy TEXT, family TEXT, genus TEXT, species TEXT, length INTEGER, sequence_sha256 TEXT)")
        counts = Counter()
        eligible = Counter()
        excluded_empty = Counter()
        for source, parser_fn in (("bold", bold_record), ("unite", unite_record)):
            path = sources / f"{source}.fasta"
            batch = []
            for offset, header, sequence in fasta_records(path):
                parsed = parser_fn(header)
                if parsed is None:
                    continue
                accession, marker, taxonomy, family, genus, species = parsed
                sequence = sequence.upper()
                if not accession or not sequence:
                    excluded_empty[marker] += 1
                    continue
                batch.append((source, offset, accession, marker, taxonomy, family,
                              genus, species, len(sequence), hashlib.sha256(sequence).hexdigest()))
                counts[marker] += 1
                if family and genus and species and len(sequence) >= 350:
                    eligible[marker] += 1
                if len(batch) >= 10000:
                    conn.executemany("INSERT INTO records VALUES (?,?,?,?,?,?,?,?,?,?)", batch)
                    conn.commit()
                    batch = []
            if batch:
                conn.executemany("INSERT INTO records VALUES (?,?,?,?,?,?,?,?,?,?)", batch)
                conn.commit()
            print(f"Indexed {counts['COI-5P'] + counts['ITS']} records after {source}", flush=True)
        conn.execute("CREATE INDEX records_species ON records(marker, species)")
        conn.execute("CREATE INDEX records_genus ON records(marker, genus)")
        conn.execute("CREATE INDEX records_family ON records(marker, family)")
        conn.execute("CREATE INDEX records_sequence ON records(marker, sequence_sha256)")
        conn.commit()
    finally:
        conn.close()
    os.replace(str(stage), str(database))
    summary.write_text(json.dumps({
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_release_record": str(releases),
        "source_release_record_sha256": sha256(releases),
        "records_by_marker": dict(sorted(counts.items())),
        "eligible_pe150_by_marker": dict(sorted(eligible.items())),
        "excluded_empty_accession_or_sequence_by_marker": dict(sorted(excluded_empty.items())),
        "eligibility": "family, genus and species named; sequence length >=350 bases",
        "sqlite_path": str(database),
        "sqlite_bytes": database.stat().st_size,
        "sqlite_sha256": sha256(database),
    }, indent=2) + "\n", encoding="utf-8")
    print(f"Source index complete: {database}")


if __name__ == "__main__":
    main()
