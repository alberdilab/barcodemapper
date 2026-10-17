"""Map selected ENA aliases to mock replicates and controls before classification."""

import argparse
import csv
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs_dir", type=Path)
    args = parser.parse_args()
    inputs = args.inputs_dir.expanduser().resolve()
    selection = inputs / "selected-runs.tsv"
    output = inputs / "sample-roles.tsv"
    provenance = inputs / "sample-roles.provenance.json"
    if output.exists() or provenance.exists():
        parser.error("Sample roles already exist")
    with selection.open(newline="", encoding="utf-8") as handle:
        selected = [row for row in csv.DictReader(handle, delimiter="\t")
                    if row["decision"] == "include"]
    records = []
    for row in selected:
        alias = row["sample_alias"]
        study = row["study"]
        if study == "PRJEB45431":
            mock = re.fullmatch(r"MC([1-4])_([1-3])", alias)
            if mock:
                role = "animal_mock"
                community, replicate = mock.groups()
                expected = "28 Table 1 animal species"
                basis = "ENA sample alias matches 4 mock types x 3 replicates in Schmidt et al."
            elif alias == "MC_C1":
                role = "extraction_blank"
                community = replicate = ""
                expected = "none"
                basis = "Inference: sole extra run, control alias, smaller file, and paper reports one extraction blank"
            else:
                raise ValueError(f"Unknown animal alias: {alias}")
        elif study == "PRJNA725368":
            if alias == "Pooled_DNA_illumina_metagenome":
                role = "fungal_pooled_DNA"
            elif alias == "Pooled_biomass_illumina_metagenome":
                role = "fungal_pooled_biomass"
            else:
                raise ValueError(f"Unknown fungal alias: {alias}")
            community = replicate = ""
            expected = "43 Table 1 fungal rows; Purpureocillium flagged as possible contaminant"
            basis = "ENA sample alias and Hu et al. methods"
        else:
            raise ValueError(f"Unexpected study: {study}")
        records.append({
            "study": study,
            "run_accession": row["run_accession"],
            "sample_alias": alias,
            "role": role,
            "community": community,
            "replicate": replicate,
            "expected_member_set": expected,
            "role_basis": basis,
        })
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    provenance.write_text(json.dumps({
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "selected_runs_sha256": sha256(selection),
        "sample_roles_sha256": sha256(output),
        "animal_paper": "https://doi.org/10.1002/ece3.8991",
        "fungal_paper": "https://doi.org/10.1128/mbio.02444-21",
    }, indent=2) + "\n", encoding="utf-8")
    print(f"Assigned roles to {len(records)} selected runs")


if __name__ == "__main__":
    main()
