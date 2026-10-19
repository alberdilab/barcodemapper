import csv
import importlib.util
import json
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "benchmarks" / "mjolnir"


def load(name):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


build_ncbi_path_map = load("build_ncbi_path_map")
build_nt_accession_paths = load("build_nt_accession_paths")
convert_calls = load("convert_calls")


def rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_nt_builders_feed_both_converters_and_keep_no_hit_pairs(tmp_path):
    ktaxonomy = tmp_path / "ktaxonomy.tsv"
    ktaxonomy.write_text(
        "1\t|\t1\t|\tR\t|\t0\t|\troot\t|\n"
        "2759\t|\t1\t|\tD\t|\t1\t|\tEukaryota\t|\n"
        "33154\t|\t2759\t|\tD1\t|\t2\t|\tOpisthokonta\t|\n"
        "33208\t|\t33154\t|\tD2\t|\t3\t|\tMetazoa\t|\n"
        "3\t|\t33208\t|\tP\t|\t4\t|\tArthropoda\t|\n"
        "4\t|\t3\t|\tC\t|\t4\t|\tInsecta\t|\n")
    taxid_paths = tmp_path / "taxid-paths.tsv"
    kingdom_taxids = tmp_path / "kingdom-taxids.txt"
    kingdom_taxids.write_text("33208\n")
    build_ncbi_path_map.main(["--ktaxonomy", str(ktaxonomy), "--kingdom-taxids",
                              str(kingdom_taxids), "--output", str(taxid_paths)])

    blast = tmp_path / "blast.tsv"
    blast.write_text(
        "r1/1\tacc|nt\t99\t100\t0\t0\t1\t100\t1\t100\t1e-20\t200\t4\n"
        "r2/1\tacc2\t99\t100\t0\t0\t1\t100\t1\t100\t1e-20\t200\t3\n")
    accession_paths = tmp_path / "accession-paths.tsv"
    build_nt_accession_paths.main(["--blast", str(blast), "--taxid-paths",
                                   str(taxid_paths), "--output", str(accession_paths)])

    retained = tmp_path / "read-ids.txt"
    retained.write_text("r1\nr2\nr3\n")
    kraken = tmp_path / "kraken.tsv"
    kraken.write_text("C\tr1/1\t4\t100\t4:20\nU\tr2/1\t0\t100\t0:20\n")
    kraken_output = tmp_path / "kraken-calls.tsv"
    blast_output = tmp_path / "blast-calls.tsv"
    convert_calls.main(["kraken2", "--read-ids", str(retained), "--input", str(kraken),
                        "--taxid-paths", str(taxid_paths), "--output", str(kraken_output)])
    convert_calls.main(["blastn", "--read-ids", str(retained), "--input", str(blast),
                        "--accession-paths", str(accession_paths), "--output", str(blast_output)])

    kraken_rows = rows(kraken_output)
    blast_rows = rows(blast_output)
    assert [row["ReadID"] for row in kraken_rows] == ["r1", "r2", "r3"]
    assert [row["ReadID"] for row in blast_rows] == ["r1", "r2", "r3"]
    assert [(row["Taxonomy"], row["status"]) for row in kraken_rows] == [
        ("k__Animalia;p__Arthropoda;c__Insecta", "classified"),
        ("", "unclassified"), ("", "no_hit")]
    assert [(row["Taxonomy"], row["status"]) for row in blast_rows] == [
        ("k__Animalia;p__Arthropoda;c__Insecta", "classified"),
        ("k__Animalia;p__Arthropoda", "classified"), ("", "no_hit")]


def test_empty_path_blast_hits_are_kept_through_builder_and_converter(tmp_path):
    ktaxonomy = tmp_path / "ktaxonomy.tsv"
    ktaxonomy.write_text(
        "1\t|\t1\t|\tR\t|\t0\t|\troot\t|\n"
        "32630\t|\t1\t|\tR\t|\t1\t|\tsynthetic construct\t|\n"
        "2759\t|\t1\t|\tD\t|\t1\t|\tEukaryota\t|\n"
        "33154\t|\t2759\t|\tD1\t|\t2\t|\tOpisthokonta\t|\n"
        "33208\t|\t33154\t|\tD2\t|\t3\t|\tMetazoa\t|\n")
    taxid_paths = tmp_path / "taxid-paths.tsv"
    kingdom_taxids = tmp_path / "kingdom-taxids.txt"
    kingdom_taxids.write_text("33208\n")
    build_ncbi_path_map.main(["--ktaxonomy", str(ktaxonomy), "--kingdom-taxids",
                              str(kingdom_taxids), "--output", str(taxid_paths)])

    def hit(read, subject, score, taxid):
        return (f"{read}\t{subject}\t99\t100\t0\t0\t1\t100\t1\t100\t"
                f"1e-20\t{score}\t{taxid}\n")

    blast = tmp_path / "blast.tsv"
    blast.write_text(hit("r1/1", "empty", 200, "32630") +
                     hit("r2/1", "classified", 200, "33208") +
                     hit("r2/1", "empty", 200, "32630") +
                     hit("r3/1", "classified", 200, "33208") +
                     hit("r3/1", "empty", 100, "32630"))
    accession_paths = tmp_path / "accession-paths.tsv"
    audit = tmp_path / "accession-audit.json"
    build_nt_accession_paths.main(["--blast", str(blast), "--taxid-paths", str(taxid_paths),
                                   "--output", str(accession_paths), "--audit", str(audit)])
    with accession_paths.open(newline="") as handle:
        paths = {row["Accession"]: row["Taxonomy"] for row in
                 csv.DictReader(handle, delimiter="\t")}
    assert paths == {"empty": "", "classified": "k__Animalia"}
    assert json.loads(audit.read_text())["empty_path_subjects"] == 1

    retained = tmp_path / "read-ids.txt"
    retained.write_text("r1\nr2\nr3\nr4\n")
    converted = tmp_path / "converted.tsv"
    convert_calls.main(["blastn", "--read-ids", str(retained), "--input", str(blast),
                        "--accession-paths", str(accession_paths), "--output", str(converted)])
    assert [(row["Taxonomy"], row["status"]) for row in rows(converted)] == [
        ("", "unclassified"), ("", "unclassified"),
        ("k__Animalia", "classified"), ("", "no_hit")]


def test_unmapped_blast_taxid_becomes_unclassified(tmp_path):
    taxid_paths = tmp_path / "taxid-paths.tsv"
    taxid_paths.write_text("TaxID\tTaxonomy\n10\tk__Animalia\n")
    blast = tmp_path / "blast.tsv"
    blast.write_text("r1/1\tnt-unknown\t99\t100\t0\t0\t1\t100\t1\t100\t1e-20\t200\t999\n")
    accession_paths = tmp_path / "accession-paths.tsv"
    audit = tmp_path / "accession-audit.json"
    build_nt_accession_paths.main([
        "--blast", str(blast), "--taxid-paths", str(taxid_paths),
        "--output", str(accession_paths), "--audit", str(audit),
        "--unmapped-as-empty"])
    assert rows(accession_paths) == [{"Accession": "nt-unknown", "Taxonomy": ""}]
    assert json.loads(audit.read_text())["unmapped_as_empty_subject_count"] == 1

    retained = tmp_path / "read-ids.txt"
    retained.write_text("r1\n")
    converted = tmp_path / "calls.tsv"
    convert_calls.main(["blastn", "--read-ids", str(retained), "--input", str(blast),
                        "--accession-paths", str(accession_paths), "--output", str(converted)])
    assert rows(converted) == [{"ReadID": "r1", "Taxonomy": "", "status": "unclassified"}]
