"""Freeze member names from the two open-access mock-community Table 1s."""

import argparse
import csv
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from xml.etree import ElementTree


SOURCES = (
    ("PRJEB45431", "PMC9170594"),
    ("PRJNA725368", "PMC9040722"),
)


def normalized(element):
    return " ".join("".join(element.itertext()).split())


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def extract(study, paper, path):
    root = ElementTree.parse(str(path)).getroot()
    doi = next((node.text for node in root.findall('.//article-id')
                if node.attrib.get('pub-id-type') == 'doi'), None)
    if not doi:
        raise ValueError(f"No DOI in {path}")
    tables = root.findall('.//table-wrap')
    if not tables:
        raise ValueError(f"No tables in {path}")
    rows = tables[0].findall('.//tr')
    expected_rows = 36 if study == "PRJEB45431" else 45
    if len(rows) != expected_rows:
        raise ValueError(f"Unexpected Table 1 size in {path}: {len(rows)}")
    members = []
    group = ""
    for index, row in enumerate(rows[2:] if study == "PRJEB45431" else rows[1:], 1):
        cells = [item for item in row if item.tag in ("th", "td")]
        if len(cells) == 1 and study == "PRJEB45431":
            group = normalized(cells[0])
            continue
        label = normalized(cells[0])
        name_match = re.match(r"^([A-Z][a-z]+ [a-z]+)\b", label)
        if not name_match:
            raise ValueError(f"Cannot extract binomial from {label!r}")
        species = name_match.group(1)
        variant = re.search(r"\b(VG[IVX]+|VN[IVX]+)\b", label)
        member = {
            "study": study,
            "paper_row": index,
            "paper_label": label,
            "species": species,
            "animal_group": group if study == "PRJEB45431" else "",
            "strain_group": variant.group(1) if variant else "",
            "expected_in_mock": "no" if species == "Purpureocillium lilacinum" else "yes",
            "note": "paper marks possible contamination" if species == "Purpureocillium lilacinum" else "",
            "paper_doi": doi,
            "paper_table": "Table 1",
        }
        members.append(member)
    expected = 28 if study == "PRJEB45431" else 44
    if len(members) != expected:
        raise ValueError(f"Expected {expected} member rows, got {len(members)}")
    return members, doi


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("animal_xml", type=Path)
    parser.add_argument("fungal_xml", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / "expected-members.tsv"
    provenance = output_dir / "expected-members.provenance.json"
    if output.exists() or provenance.exists():
        parser.error("Expected member snapshot already exists")
    paths = (args.animal_xml, args.fungal_xml)
    records = []
    members = []
    for (study, paper), path in zip(SOURCES, paths):
        group, doi = extract(study, paper, path)
        members.extend(group)
        records.append({"study": study, "paper": paper, "doi": doi,
                        "url": f"https://www.ebi.ac.uk/europepmc/webservices/rest/{paper}/fullTextXML",
                        "xml_sha256": sha256(path), "member_rows": len(group)})
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(members[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(members)
    provenance.write_text(json.dumps({
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "sources": records,
        "expected_members_sha256": sha256(output),
        "interpretation": "Purpureocillium lilacinum is a possible contaminant in the paper and is excluded from the primary expected-taxon denominator. Cryptococcus strain variants remain separate rows, but share a species-level name.",
    }, indent=2) + "\n", encoding="utf-8")
    print(f"Extracted {len(members)} paper rows; {sum(row['expected_in_mock'] == 'yes' for row in members)} expected")


if __name__ == "__main__":
    main()
