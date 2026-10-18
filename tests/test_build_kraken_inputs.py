"""Small offline fixtures for deterministic Kraken2 reference conversion."""

import gzip
import importlib.util
import json
from pathlib import Path

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "benchmarks/mjolnir/build_kraken_inputs.py"
SPEC = importlib.util.spec_from_file_location("build_kraken_inputs", SCRIPT)
kraken = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(kraken)

RECORDS = (
    ("A1", "COI-5P", "k__Animalia;p__Arthropoda;c__Insecta", "ACGT"),
    ("A2", "COI-5P", "k__Animalia;p__Arthropoda;c__Insecta", "TTAA"),
    ("A3", "COI-5P", "k__Animalia;p__Chordata", "ACGT"),
    ("F1", "ITS", "k__Fungi;p__Ascomycota", "GGCC"),
)


def fixture_fasta(path, records):
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for accession, marker, taxonomy, sequence in records:
            handle.write(f">{accession}|{marker}|{taxonomy}\n{sequence}\n")


def taxid_paths(output):
    rows = (output / "taxid-to-path.tsv").read_text(encoding="utf-8").splitlines()
    return {path: int(taxid) for taxid, path in (row.split("\t") for row in rows[1:])}


def test_stable_taxids_and_preserved_conflicting_records(tmp_path):
    first = tmp_path / "first.fa.gz"
    second = tmp_path / "second.fa.gz"
    fixture_fasta(first, RECORDS)
    fixture_fasta(second, reversed(RECORDS))
    out1, out2 = tmp_path / "out1", tmp_path / "out2"
    kraken.build_inputs(first, out1)
    kraken.build_inputs(second, out2)

    assert taxid_paths(out1) == taxid_paths(out2)
    paths = taxid_paths(out1)
    assert paths["k__Animalia"] < paths["k__Animalia;p__Arthropoda"]
    assert paths["k__Animalia;p__Arthropoda"] < paths["k__Animalia;p__Arthropoda;c__Insecta"]
    assert paths["k__Animalia;p__Chordata"] != paths["k__Animalia;p__Arthropoda"]
    library = (out1 / "library.fa").read_text(encoding="utf-8")
    assert library.count(">kraken:taxid|") == len(RECORDS)
    assert library.count("ACGT\n") == 2  # Same sequence, distinct taxonomy labels.
    assert len((out1 / "seqid2taxid.map").read_text().splitlines()) == len(RECORDS)
    assert (out1 / "taxonomy/nodes.dmp").is_file()
    assert (out1 / "taxonomy/names.dmp").is_file()


def test_count_audit_passes_and_detects_missing_record(tmp_path):
    fasta = tmp_path / "reference.fa.gz"
    fixture_fasta(fasta, RECORDS)
    output = tmp_path / "kraken"
    audit_path = kraken.build_inputs(fasta, output)
    audit = json.loads(audit_path.read_text())
    assert audit["source_records"] == audit["library_records"] == audit["map_records"] == 4
    assert audit["source_records_by_marker"] == {"COI-5P": 3, "ITS": 1}
    assert audit["library_records_by_marker"] == audit["source_records_by_marker"]
    assert audit["map_records_by_marker"] == audit["source_records_by_marker"]
    assert audit["outputs"]["library.fa"] == kraken.sha256(output / "library.fa")

    library = output / "library.fa"
    lines = library.read_text().splitlines(keepends=True)
    library.write_text("".join(lines[:-2]))
    with pytest.raises(ValueError, match="record counts differ"):
        kraken.audit_outputs(kraken.source_inventory(fasta)[0], library,
                             output / "seqid2taxid.map")


def test_build_does_not_publish_audit_on_count_mismatch(tmp_path, monkeypatch):
    fasta = tmp_path / "reference.fa.gz"
    fixture_fasta(fasta, RECORDS)
    original = kraken.write_library

    def omit_last(fasta_path, output, taxids):
        original(fasta_path, output, taxids)
        library = output / "library.fa"
        library.write_text("".join(library.read_text().splitlines(keepends=True)[:-2]))
        mapping = output / "seqid2taxid.map"
        mapping.write_text("".join(mapping.read_text().splitlines(keepends=True)[:-1]))

    monkeypatch.setattr(kraken, "write_library", omit_last)
    output = tmp_path / "kraken"
    with pytest.raises(ValueError, match="Reference record count mismatch"):
        kraken.build_inputs(fasta, output)
    assert not (output / "audit.json").exists()
