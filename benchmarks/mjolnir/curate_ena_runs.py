"""Freeze sample metadata and classify the ENA mock-community runs."""

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen
from xml.etree import ElementTree


STUDIES = ("PRJEB45431", "PRJNA725368")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs_dir", type=Path)
    args = parser.parse_args()
    inputs = args.inputs_dir.expanduser().resolve()
    output = inputs / "selected-runs.tsv"
    provenance = inputs / "selected-runs.provenance.json"
    if output.exists() or provenance.exists():
        parser.error("Selection already exists; preserve the frozen decision")
    xml_dir = inputs / "sample-xml"
    xml_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    sample_records = []
    for study in STUDIES:
        report = inputs / f"ena-{study}.tsv"
        with report.open(newline="", encoding="utf-8") as handle:
            runs = list(csv.DictReader(handle, delimiter="\t"))
        for row in runs:
            sample = row["sample_accession"]
            xml_path = xml_dir / f"{sample}.xml"
            url = f"https://www.ebi.ac.uk/ena/browser/api/xml/{sample}"
            if not xml_path.exists():
                with urlopen(Request(url, headers={"User-Agent": "BarcodeMapper-benchmark/1.1.0"}), timeout=60) as response:
                    body = response.read()
                xml_path.write_bytes(body)
            body = xml_path.read_bytes()
            xml = ElementTree.fromstring(body)
            sample_element = xml.find("SAMPLE")
            if sample_element is None or sample_element.attrib.get("accession") != sample:
                raise ValueError(f"Unexpected XML for {sample}")
            alias = sample_element.attrib.get("alias", "")
            description = (sample_element.findtext("DESCRIPTION") or "").strip()
            sample_records.append({"sample": sample, "url": url,
                                   "sha256": hashlib.sha256(body).hexdigest()})
            instrument = row["instrument_model"].casefold()
            include = (row["library_strategy"] == "WGS"
                       and row["library_layout"] == "PAIRED"
                       and "illumina" in instrument)
            if study == "PRJEB45431":
                sample_type = "animal mock or extraction control; verify against paper sample key"
                member_table = "PRJEB45431 paper Table 1"
            else:
                sample_type = alias
                member_table = "Hu et al. 2022 Table 1"
            rows.append({
                "study": study,
                "run_accession": row["run_accession"],
                "sample_accession": sample,
                "sample_alias": alias,
                "sample_description": description,
                "sample_type": sample_type,
                "library_strategy": row["library_strategy"],
                "library_layout": row["library_layout"],
                "instrument_model": row["instrument_model"],
                "reported_member_table": member_table if include else "",
                "r1_url": "https://" + row["fastq_ftp"].split(";")[0] if include else "",
                "r2_url": "https://" + row["fastq_ftp"].split(";")[1] if include else "",
                "r1_ena_md5": row["fastq_md5"].split(";")[0] if include else "",
                "r2_ena_md5": row["fastq_md5"].split(";")[1] if include else "",
                "decision": "include" if include else "exclude",
                "reason": "paired Illumina WGS" if include else "not paired Illumina WGS",
            })
    fields = list(rows[0])
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    provenance.write_text(json.dumps({
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "selection_rule": "library_strategy=WGS, layout=PAIRED, Illumina instrument",
        "sample_xml": sample_records,
        "ena_report_sha256": {study: hashlib.sha256((inputs / f"ena-{study}.tsv").read_bytes()).hexdigest()
                              for study in STUDIES},
        "selected_runs_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
    }, indent=2) + "\n", encoding="utf-8")
    print(f"Selected {sum(row['decision'] == 'include' for row in rows)} of {len(rows)} runs")
    for row in rows:
        print(row["run_accession"], row["sample_alias"], row["decision"])


if __name__ == "__main__":
    main()
