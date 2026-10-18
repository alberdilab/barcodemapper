import csv
import importlib.util
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
        "2\t|\t2759\t|\tK\t|\t2\t|\tMetazoa\t|\n"
        "3\t|\t2\t|\tP\t|\t3\t|\tArthropoda\t|\n"
        "4\t|\t3\t|\tC\t|\t4\t|\tInsecta\t|\n")
    taxid_paths = tmp_path / "taxid-paths.tsv"
    build_ncbi_path_map.main(["--ktaxonomy", str(ktaxonomy), "--output", str(taxid_paths)])

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
