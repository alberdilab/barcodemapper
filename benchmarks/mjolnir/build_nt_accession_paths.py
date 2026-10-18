"""Build accession paths from BLAST nt hits and a common taxid path table."""

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
from tempfile import NamedTemporaryFile


RANKS = "kpcofgs"


def open_text(path):
    opener = gzip.open if path.suffix == ".gz" else open
    return opener(path, "rt", encoding="utf-8", newline="")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def resolved_prefix(taxonomy):
    levels = []
    for index, level in enumerate(taxonomy.split(";")):
        if not level or level.endswith("__"):
            break
        if index >= len(RANKS) or not level.startswith(f"{RANKS[index]}__"):
            raise ValueError(f"Invalid taxonomy path {taxonomy!r}")
        levels.append(level)
    return ";".join(levels)


def common_taxonomy(first, second):
    shared = []
    for left, right in zip(first.split(";"), second.split(";")):
        if left != right:
            break
        shared.append(left)
    return ";".join(shared)


def load_paths(path):
    paths = {}
    with open_text(path) as handle:
        for number, row in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if number == 1 and row == ["TaxID", "Taxonomy"]:
                continue
            if len(row) != 2 or not row[0]:
                raise ValueError(f"Invalid taxid path at {path}:{number}")
            taxonomy = resolved_prefix(row[1])
            if row[1] and not taxonomy:
                raise ValueError(f"Unresolved taxid path at {path}:{number}")
            if row[0] in paths and paths[row[0]] != taxonomy:
                raise ValueError(f"Conflicting taxid path for {row[0]!r}")
            paths[row[0]] = taxonomy
    return paths


def build_paths(blast, taxid_paths, output, audit, allow_unmapped=False,
                unmapped_as_empty=False):
    if allow_unmapped and unmapped_as_empty:
        raise ValueError("--allow-unmapped and --unmapped-as-empty are mutually exclusive")
    if output == audit:
        raise ValueError("Output paths must be distinct")
    for path in (output, audit):
        if path.exists():
            raise ValueError(f"Output already exists: {path}")
    paths = load_paths(taxid_paths)
    subjects = {}
    subject_order = {}
    unmapped = {}
    empty_path_subjects = set()
    row_count = 0
    multi_taxid_rows = 0
    multi_taxid_subjects = set()
    no_taxid_rows = 0
    with open_text(blast) as handle:
        for number, row in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 13 or not row[1]:
                raise ValueError(f"Invalid BLAST row at {blast}:{number}; expected staxids in column 13")
            row_count += 1
            subject = row[1]
            subject_order.setdefault(subject, None)
            taxids = [taxid.strip() for taxid in row[12].split(";")]
            if len(taxids) > 1:
                multi_taxid_rows += 1
                multi_taxid_subjects.add(subject)
            if taxids == ["N/A"] or taxids == ["0"] or taxids == [""]:
                no_taxid_rows += 1
            missing = [taxid for taxid in taxids if taxid not in paths]
            empty = any(taxid in paths and not paths[taxid] for taxid in taxids)
            if empty:
                empty_path_subjects.add(subject)
            if missing:
                unmapped.setdefault(subject, set()).update(missing)
                continue
            if empty:
                subjects[subject] = ""
                continue
            taxonomy = paths[taxids[0]]
            for taxid in taxids[1:]:
                taxonomy = common_taxonomy(taxonomy, paths[taxid])
            if subject in empty_path_subjects:
                taxonomy = ""
            elif subject in subjects:
                taxonomy = common_taxonomy(subjects[subject], taxonomy)
            subjects[subject] = taxonomy

    for subject, taxonomy in subjects.items():
        if not taxonomy and subject not in empty_path_subjects:
            unmapped.setdefault(subject, set()).add("no common kingdom")
    missing_taxids = {taxid for taxids in unmapped.values() for taxid in taxids
                      if taxid != "no common kingdom"}
    if unmapped and not (allow_unmapped or unmapped_as_empty):
        raise ValueError(f"{len(unmapped)} BLAST subjects have unmapped taxids or no common kingdom "
                         f"({len(missing_taxids)} distinct unmapped taxids; e.g. {min(unmapped)!r}); "
                         "pass --allow-unmapped to omit them or --unmapped-as-empty to classify them as unclassified")

    kept = (subjects.keys() | unmapped.keys()) if unmapped_as_empty else subjects.keys() - unmapped.keys()
    audit_data = {
        "input_sha256": {"blast": sha256(blast), "taxid_paths": sha256(taxid_paths)},
        "blast_rows": row_count,
        "subject_count": len(subjects.keys() | unmapped.keys() | empty_path_subjects),
        "output_rows": len(kept),
        "multi_taxid_rows": multi_taxid_rows,
        "multi_taxid_subject_count": len(multi_taxid_subjects),
        "no_taxid_rows": no_taxid_rows,
        "unmapped_subject_count": len(unmapped),
        "unmapped_taxid_count": len(missing_taxids),
        "unmapped_subjects": {subject: sorted(taxids) for subject, taxids in sorted(unmapped.items())},
        "empty_path_subjects": len(empty_path_subjects),
        "empty_path_subject_ids": sorted(empty_path_subjects),
    }
    if unmapped_as_empty:
        audit_data["unmapped_as_empty_subject_count"] = len(unmapped)
        audit_data["unmapped_as_empty_subjects"] = {
            subject: sorted(
                "no common kingdom" if taxid == "no common kingdom" else
                f"no taxid ({taxid or 'empty'})" if taxid in ("N/A", "0", "") else
                f"taxid {taxid} absent from taxid paths"
                for taxid in taxids
            )
            for subject, taxids in sorted(unmapped.items())
        }
    temporary = {}
    published = []
    try:
        with NamedTemporaryFile("w", dir=output.parent, prefix=f".{output.name}.",
                                suffix=".tmp", newline="", encoding="utf-8",
                                delete=False) as handle:
            temporary[output] = Path(handle.name)
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(("Accession", "Taxonomy"))
            if unmapped_as_empty:
                for subject in subject_order:
                    writer.writerow((subject, "" if subject in unmapped else subjects[subject]))
            else:
                for subject, taxonomy in subjects.items():
                    if subject in kept:
                        writer.writerow((subject, taxonomy))
        with NamedTemporaryFile("w", dir=audit.parent, prefix=f".{audit.name}.",
                                suffix=".tmp", encoding="utf-8", delete=False) as handle:
            temporary[audit] = Path(handle.name)
            json.dump(audit_data, handle, indent=2)
            handle.write("\n")
        for path, staged in temporary.items():
            if path.exists():
                raise ValueError(f"Output already exists: {path}")
            staged.replace(path)
            published.append(path)
    except Exception:
        for staged in temporary.values():
            staged.unlink(missing_ok=True)
        for path in published:
            path.unlink(missing_ok=True)
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blast", required=True, type=Path)
    parser.add_argument("--taxid-paths", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--audit", type=Path)
    policy = parser.add_mutually_exclusive_group()
    policy.add_argument("--allow-unmapped", action="store_true")
    policy.add_argument("--unmapped-as-empty", action="store_true")
    args = parser.parse_args(argv)
    audit = args.audit or args.output.with_suffix(args.output.suffix + ".audit.json")
    build_paths(args.blast, args.taxid_paths, args.output, audit,
                args.allow_unmapped, args.unmapped_as_empty)


if __name__ == "__main__":
    main()
