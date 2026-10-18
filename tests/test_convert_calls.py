import csv
import gzip
import importlib.util
from pathlib import Path

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "benchmarks" / "mjolnir" / "convert_calls.py"
spec = importlib.util.spec_from_file_location("convert_calls", SCRIPT)
convert_calls = importlib.util.module_from_spec(spec)
spec.loader.exec_module(convert_calls)

KRAKEN_SCRIPT = SCRIPT.with_name("build_kraken_inputs.py")
kraken_spec = importlib.util.spec_from_file_location("build_kraken_inputs", KRAKEN_SCRIPT)
build_kraken_inputs = importlib.util.module_from_spec(kraken_spec)
kraken_spec.loader.exec_module(build_kraken_inputs)


def convert(tmp_path, method, calls, extra=(), read_ids="r1\nr2\nr3\n"):
    retained = tmp_path / "read_ids.txt"
    retained.write_text(read_ids)
    source = tmp_path / "calls.tsv"
    source.write_text(calls)
    output = tmp_path / "converted.tsv"
    convert_calls.main([method, "--read-ids", str(retained), "--input", str(source),
                        "--output", str(output), *extra])
    with output.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_barcodemapper_keeps_all_pairs_and_resolves_prefix(tmp_path):
    rows = convert(tmp_path, "barcodemapper",
                   "ReadID\tTaxonomy\nr1/1\tk__Fungi;p__Ascomycota;g__;s__Unknown\n"
                   "r2/2\t\n")
    assert rows == [
        {"ReadID": "r1", "Taxonomy": "k__Fungi;p__Ascomycota", "status": "classified"},
        {"ReadID": "r2", "Taxonomy": "", "status": "unclassified"},
        {"ReadID": "r3", "Taxonomy": "", "status": "no_hit"},
    ]


def test_best_hit_translates_source_statuses(tmp_path):
    rows = convert(tmp_path, "best-hit",
                   "ReadID\tTaxonomy\tstatus\tAccession\n"
                   "r1\tk__Animalia;g__A\taccepted_alignment\ta\n"
                   "r2\t\tno_accepted_alignment\t\n")
    assert [(row["Taxonomy"], row["status"]) for row in rows] == [
        ("k__Animalia;g__A", "classified"), ("", "no_hit"), ("", "no_hit")]


def test_kraken2_maps_prefix_taxids_and_unclassified(tmp_path):
    mapping = tmp_path / "taxids.tsv"
    mapping.write_text("TaxID\tTaxonomy\n10\tk__Fungi\n11\tk__Fungi;p__Ascomycota\n")
    rows = convert(tmp_path, "kraken2",
                   "C\tr1/1\t11\t100\t11:20\nU\tr2/2\t0\t100\t0:20\n"
                   "C\tr3\t0\t100\t0:20\n",
                   ["--taxid-paths", str(mapping)], read_ids="r1\nr2\nr3\nr4\n")
    assert [(row["Taxonomy"], row["status"]) for row in rows] == [
        ("k__Fungi;p__Ascomycota", "classified"), ("", "unclassified"),
        ("", "unclassified"), ("", "no_hit")]


def test_kraken2_accepts_generated_root_taxid_map(tmp_path):
    fasta = tmp_path / "reference.fa.gz"
    with gzip.open(fasta, "wt", encoding="utf-8") as handle:
        handle.write(">a|ITS|k__Fungi;p__Ascomycota\nACGT\n")
    output = tmp_path / "kraken"
    build_kraken_inputs.build_inputs(fasta, output)
    paths = convert_calls.load_paths(output / "taxid-to-path.tsv", "TaxID")
    taxid = next(key for key, value in paths.items() if value == "k__Fungi;p__Ascomycota")
    rows = convert(tmp_path, "kraken2", f"C\tr1\t{taxid}\t4\nC\tr2\t1\t4\n",
                   ["--taxid-paths", str(output / "taxid-to-path.tsv")])
    assert [(row["Taxonomy"], row["status"]) for row in rows] == [
        ("k__Fungi;p__Ascomycota", "classified"), ("", "unclassified"),
        ("", "no_hit")]


@pytest.mark.parametrize("method,hit", [
    ("blastn", "{read}\t{subject}\t99\t100\t0\t0\t1\t100\t1\t100\t1e-20\t100\n"),
    ("vsearch", "{read}\t{subject}\t99\n"),
    ("minimap2", "{read}\t100\t0\t100\t+\t{subject}\t100\t0\t100\t100\t100\t60\n"),
])
def test_alignment_formats_reduce_all_accepted_hits_to_lca(tmp_path, method, hit):
    mapping = tmp_path / "accessions.tsv"
    mapping.write_text("Accession\tTaxonomy\n"
                       "a\tk__Fungi;g__Alpha;s__one\n"
                       "b\tk__Fungi;g__Alpha;s__two\n"
                       "c\tk__Animalia;g__Beta\n")
    calls = (hit.format(read="r1/1", subject="a|ITS|ignored") +
             hit.format(read="r1/2", subject="b") +
             hit.format(read="r1/2", subject="a") +
             hit.format(read="r2", subject="a") +
             hit.format(read="r2", subject="c"))
    rows = convert(tmp_path, method, calls, ["--accession-paths", str(mapping)])
    assert [(row["Taxonomy"], row["status"]) for row in rows] == [
        ("k__Fungi;g__Alpha", "classified"), ("", "unclassified"),
        ("", "no_hit")]


def test_blast_uses_only_top_score_hits_per_mate(tmp_path):
    mapping = tmp_path / "accessions.tsv"
    mapping.write_text("Accession\tTaxonomy\na\tk__Fungi;g__Alpha;s__one\n"
                       "b\tk__Fungi;g__Beta;s__two\n")
    def hit(read, subject, score):
        return f"{read}\t{subject}\t99\t100\t0\t0\t1\t100\t1\t100\t1e-20\t{score}\n"
    rows = convert(tmp_path, "blastn", hit("r1/1", "a", 200) + hit("r1/1", "b", 100),
                   ["--accession-paths", str(mapping)])
    assert (rows[0]["Taxonomy"], rows[0]["status"]) == ("k__Fungi;g__Alpha;s__one", "classified")


def test_blast_capped_top_tie_is_unclassified(tmp_path):
    mapping = tmp_path / "accessions.tsv"
    mapping.write_text("Accession\tTaxonomy\na\tk__Fungi;g__Alpha;s__one\n"
                       "b\tk__Fungi;g__Beta;s__two\n")
    config = tmp_path / "blast_config.yaml"
    config.write_text("max_target_seqs: 2\nlca_tie_policy:\n  score: bit_score\n"
                      "  ties: exact_equal_score\n")
    hit = "r1/1\ta\t99\t100\t0\t0\t1\t100\t1\t100\t1e-20\t200\n"
    rows = convert(tmp_path, "blastn", hit + hit.replace("\ta\t", "\tb\t"),
                   ["--accession-paths", str(mapping), "--blast-config", str(config)])
    assert (rows[0]["Taxonomy"], rows[0]["status"]) == ("", "unclassified")


def test_blast_cap_counts_targets_not_hsp_rows(tmp_path):
    mapping = tmp_path / "accessions.tsv"
    mapping.write_text("Accession\tTaxonomy\na\tk__Fungi;g__Alpha;s__one\n")
    config = tmp_path / "blast_config.yaml"
    config.write_text("max_target_seqs: 2\nlca_tie_policy:\n  score: bit_score\n"
                      "  ties: exact_equal_score\n")
    hit = "r1/1\ta\t99\t100\t0\t0\t1\t100\t1\t100\t1e-20\t200\n"
    rows = convert(tmp_path, "blastn", hit * 2,
                   ["--accession-paths", str(mapping), "--blast-config", str(config)])
    assert (rows[0]["Taxonomy"], rows[0]["status"]) == (
        "k__Fungi;g__Alpha;s__one", "classified")


def test_duplicate_retained_ids_are_rejected(tmp_path):
    with pytest.raises(ValueError, match="Duplicate retained read pair ID"):
        convert(tmp_path, "barcodemapper", "ReadID\tTaxonomy\n", read_ids="r1/1\nr1/2\n")


def test_unmapped_accession_is_rejected(tmp_path):
    mapping = tmp_path / "accessions.tsv"
    mapping.write_text("Accession\tTaxonomy\na\tk__Fungi\n")
    with pytest.raises(ValueError, match="Unmapped accession"):
        convert(tmp_path, "blastn", "r1\tmissing\t99\t100\t0\t0\t1\t100\t1\t100\t1e-20\t100\n",
                ["--accession-paths", str(mapping)])


def test_fastp_r1_fastq_can_supply_retained_ids(tmp_path):
    retained = tmp_path / "retained_1.fq.gz"
    with gzip.open(retained, "wt") as handle:
        handle.write("@r1/1 1:N:0:1\nAC\n+\nII\n@r2/1\nGT\n+\nII\n")
    source = tmp_path / "calls.tsv"
    source.write_text("ReadID\tTaxonomy\nr1\tk__Fungi\n")
    output = tmp_path / "converted.tsv"
    convert_calls.main(["barcodemapper", "--read-ids", str(retained),
                        "--input", str(source), "--output", str(output)])
    with output.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [(row["ReadID"], row["status"]) for row in rows] == [
        ("r1", "classified"), ("r2", "no_hit")]
