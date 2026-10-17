"""Freeze per-run mock taxa and exact reference coverage before viewing calls."""

import argparse
import csv
import hashlib
import json
import sqlite3
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


STUDY_REFERENCE = {
    "PRJEB45431": ("COI-5P", "default"),
    "PRJNA725368": ("ITS", "all-fungi"),
}


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def reference_counts(path):
    ranks = ("s", "g", "f")
    counts = {rank: Counter() for rank in ranks}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            marker = row["Marker"]
            amount = int(row["Records"])
            names = dict(level.split("__", 1) for level in row["Taxonomy"].split(";"))
            for rank in ranks:
                if names.get(rank):
                    counts[rank][(marker, names[rank])] += amount
    return counts


def source_families(connection, marker, species):
    rows = connection.execute(
        "SELECT DISTINCT family FROM records WHERE marker=? AND species=? AND family!=''",
        (marker, species),
    ).fetchall()
    return sorted(row[0] for row in rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bench_root", type=Path)
    args = parser.parse_args()
    root = args.bench_root.expanduser().resolve()
    inputs = root / "projects/01-mock-communities/inputs"
    source_index = root / "projects/02-reference-novelty/inputs/source-index.sqlite"
    members_path = inputs / "expected-members.tsv"
    roles_path = inputs / "sample-roles.tsv"
    output = inputs / "expected-taxa.tsv"
    provenance = inputs / "expected-taxa.provenance.json"
    if output.exists() or provenance.exists():
        parser.error("Expected-taxon annotation already exists")
    grouped = defaultdict(list)
    for row in read_tsv(members_path):
        if row["expected_in_mock"] == "yes":
            grouped[(row["study"], row["species"])].append(row)
    counts = {}
    audits = {}
    for study, (marker, reference) in STUDY_REFERENCE.items():
        audit = root / f"references/audits/{reference}"
        counts[marker] = reference_counts(audit / "taxonomy-paths.tsv")
        audits[reference] = {
            "summary_path": str(audit / "summary.json"),
            "summary_sha256": sha256(audit / "summary.json"),
        }
    source = sqlite3.connect(f"file:{source_index}?mode=ro", uri=True)
    family_by_species = {}
    for study, species in grouped:
        marker = STUDY_REFERENCE[study][0]
        family_by_species[(study, species)] = source_families(
            source, marker, species.replace(" ", "_"))
    source.close()
    fieldnames = (
        "study", "run_accession", "sample_alias", "role", "marker", "reference",
        "expected_species", "paper_rows", "strain_groups", "source_families",
        "exact_species_reference_records", "genus_reference_records",
        "family_reference_records", "coverage_status", "synonym_status",
    )
    rows = []
    for role in read_tsv(roles_path):
        if role["role"] == "extraction_blank":
            continue
        study = role["study"]
        marker, reference = STUDY_REFERENCE[study]
        for (member_study, species), paper_rows in sorted(grouped.items()):
            if member_study != study:
                continue
            label = species.replace(" ", "_")
            genus = species.split()[0]
            families = family_by_species[(study, species)]
            count = counts[marker]
            exact = count["s"][(marker, label)]
            genus_count = count["g"][(marker, genus)]
            family_count = sum(count["f"][(marker, family)] for family in families)
            if exact:
                coverage = "exact_species"
            elif genus_count:
                coverage = "genus_only"
            elif family_count:
                coverage = "family_only"
            elif families:
                coverage = "absent_from_reference"
            else:
                coverage = "unresolved_source_family"
            rows.append({
                "study": study,
                "run_accession": role["run_accession"],
                "sample_alias": role["sample_alias"],
                "role": role["role"],
                "marker": marker,
                "reference": reference,
                "expected_species": species,
                "paper_rows": ",".join(row["paper_row"] for row in paper_rows),
                "strain_groups": ",".join(row["strain_group"] for row in paper_rows if row["strain_group"]),
                "source_families": ",".join(families),
                "exact_species_reference_records": exact,
                "genus_reference_records": genus_count,
                "family_reference_records": family_count,
                "coverage_status": coverage,
                "synonym_status": "unreviewed; exact paper binomial only",
            })
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    stats = Counter(row["coverage_status"] for row in rows)
    provenance.write_text(json.dumps({
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script_sha256": sha256(Path(__file__)),
        "expected_members_sha256": sha256(members_path),
        "sample_roles_sha256": sha256(roles_path),
        "source_index_summary_sha256": sha256(source_index.with_suffix(".json")),
        "reference_audits": audits,
        "expected_taxa_sha256": sha256(output),
        "rows": len(rows),
        "rows_by_coverage_status": dict(sorted(stats.items())),
        "interpretation": "One row per distinct expected binomial per nonblank run. Strain rows with the same binomial are collapsed. Synonyms are not yet adjudicated; exact misses need manual review before computing recall.",
    }, indent=2) + "\n", encoding="utf-8")
    print(f"Annotated {len(rows)} expected run-taxon pairs: {dict(stats)}")


if __name__ == "__main__":
    main()
