import csv
import importlib.util
import json
from pathlib import Path

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "benchmarks" / "mjolnir" / "build_ncbi_path_map.py"
spec = importlib.util.spec_from_file_location("build_ncbi_path_map", SCRIPT)
builder = importlib.util.module_from_spec(spec)
spec.loader.exec_module(builder)


def rows(path):
    with path.open(newline="") as handle:
        return {row["TaxID"]: row["Taxonomy"] for row in csv.DictReader(handle, delimiter="\t")}


def dmp(*fields):
    return "\t|\t".join(fields) + "\t|\n"


def test_taxdump_names_lineage_merges_and_audit(tmp_path):
    taxdump = tmp_path / "taxdump"
    taxdump.mkdir()
    nodes = [("1", "1", "no rank"), ("131567", "1", "no rank"),
             ("2759", "131567", "superkingdom"), ("3", "2759", "kingdom"),
             ("4", "3", "phylum"), ("40", "4", "clade"),
             ("5", "40", "class"), ("6", "5", "order"),
             ("7", "6", "family"), ("8", "7", "genus"),
             ("9", "8", "species"), ("10", "9", "subspecies"),
             ("11", "3", "phylum"), ("12", "11", "genus"),
             ("13", "12", "species"), ("2", "131567", "superkingdom"),
             ("200", "2", "kingdom")]
    (taxdump / "nodes.dmp").write_text("".join(dmp(*row) for row in nodes))
    names = {"1": "root", "131567": "cellular organisms", "2759": "Eukaryota",
             "3": "Metazoa", "4": "Arthropoda", "40": "clade",
             "5": "Insecta", "6": "Diptera", "7": "Testidae", "8": "GenusX",
             "9": "GenusX alpha; beta", "10": "GenusX alpha; beta subsp. test",
             "11": "Mollusca", "12": "Skippus", "13": "Skippus alba",
             "2": "Bacteria", "200": "Bacillati"}
    (taxdump / "names.dmp").write_text(
        "".join(dmp(taxid, name, "", "scientific name") for taxid, name in names.items()) +
        dmp("8", "Oldgenus", "", "synonym"))
    (taxdump / "merged.dmp").write_text(dmp("90", "9"))
    (taxdump / "delnodes.dmp").write_text(dmp("99"))
    name_map = tmp_path / "name-map.tsv"
    name_map.write_text("Rank\tNCBIName\tCommonName\n"
                        "k\tMetazoa\tAnimalia\n"
                        "g\tOldgenus\tCommonGenus\n")
    output = tmp_path / "paths.tsv"
    audit = tmp_path / "audit.json"
    report = tmp_path / "missing.tsv"
    builder.main(["--taxdump", str(taxdump), "--name-map", str(name_map),
                  "--output", str(output), "--audit", str(audit), "--report", str(report)])

    paths = rows(output)
    species = ("k__Animalia;p__Arthropoda;c__Insecta;o__Diptera;f__Testidae;"
               "g__CommonGenus;s__GenusX alpha, beta")
    assert paths["9"] == paths["10"] == paths["90"] == species
    assert paths["12"] == "k__Animalia;p__Mollusca;c__;o__;f__;g__Skippus;s__"
    assert paths["13"].endswith("g__Skippus;s__Skippus alba")
    assert paths["2759"].startswith("k__Eukaryota;p__")
    assert paths["200"].startswith("k__Bacteria;p__")
    assert paths["1"] == paths["131567"] == ""
    assert "99" not in paths
    result = json.loads(audit.read_text())
    assert result["input_sha256"]["nodes"] == builder.sha256(taxdump / "nodes.dmp")
    assert result["name_map_sha256"] == builder.sha256(name_map)
    assert result["merged_count"] == result["merged_rows"] == 1
    assert result["deleted_count"] == 1
    assert result["semicolon_names_replaced"] == 1
    assert result["name_map_changes"] == 2
    assert result["taxa_lost_by_rank"]["c"] > 0
    assert result["resolved_depth_counts"]["0"] == 2
    assert result["resolved_depth_counts"]["2"] >= 3
    assert report.read_text() == "TaxID\n1\n131567\n"


def test_ktaxonomy_plain_ranks_and_default_name_map(tmp_path):
    source = tmp_path / "ktaxonomy.tsv"
    source.write_text("".join(dmp(*row) for row in [
        ("1", "1", "R", "0", "root"),
        ("131567", "1", "R", "1", "cellular organisms"),
        ("2759", "131567", "D", "2", "Eukaryota"),
        ("3", "2759", "K1", "3", "ignored"),
        ("4", "3", "K", "3", "Viridiplantae"),
        ("5", "4", "P", "4", "Tracheophyta"),
        ("6", "5", "P1", "5", "ignored"),
        ("7", "6", "C", "6", "Magnoliopsida"),
        ("20", "131567", "R1", "2", "ignored"),
        ("21", "20", "D1", "3", "ignored"),
        ("2", "21", "D", "4", "Bacteria"),
        ("200", "2", "K", "5", "Bacillati"),
    ]) + "8\t|\t7\t|\tS\t|\t7\t|\tSpecies alpha\n" +
                      "9\t|\t8\t|\tS1\t|\t8\t|\tSpecies alpha strain\n")
    output = tmp_path / "paths.tsv"
    audit = tmp_path / "audit.json"
    builder.main(["--ktaxonomy", str(source), "--output", str(output), "--audit", str(audit)])
    paths = rows(output)
    assert paths["7"].startswith("k__Plantae;p__Tracheophyta;c__Magnoliopsida;")
    assert paths["3"].startswith("k__Eukaryota;p__")
    assert paths["8"].endswith("s__Species alpha")
    assert paths["9"] == paths["8"]
    assert paths["200"].startswith("k__Bacteria;p__")
    assert paths["1"] == paths["131567"] == paths["20"] == paths["21"] == ""
    assert json.loads(audit.read_text())["name_map_changes"] == 1


def test_refuses_to_overwrite_existing_output(tmp_path):
    source = tmp_path / "ktaxonomy.tsv"
    source.write_text(dmp("1", "1", "R", "0", "root"))
    output = tmp_path / "paths.tsv"
    output.write_text("keep")
    with pytest.raises(ValueError, match="already exists"):
        builder.main(["--ktaxonomy", str(source), "--output", str(output)])
    assert output.read_text() == "keep"


def test_merged_chain_and_missing_target_are_counted(tmp_path):
    taxdump = tmp_path / "taxdump"
    taxdump.mkdir()
    (taxdump / "nodes.dmp").write_text(
        dmp("1", "1", "no rank") + dmp("2759", "1", "superkingdom") +
        dmp("3", "2759", "kingdom"))
    (taxdump / "names.dmp").write_text(
        dmp("2759", "Eukaryota", "", "scientific name") +
        dmp("3", "Metazoa", "", "scientific name"))
    (taxdump / "merged.dmp").write_text(
        dmp("90", "91") + dmp("91", "3") + dmp("92", "999") + dmp("93", "92"))
    (taxdump / "delnodes.dmp").write_text(dmp("999"))
    output = tmp_path / "paths.tsv"
    audit = tmp_path / "audit.json"
    builder.main(["--taxdump", str(taxdump), "--output", str(output), "--audit", str(audit)])
    paths = rows(output)
    assert paths["90"] == paths["91"] == paths["3"]
    assert "92" not in paths and "93" not in paths
    result = json.loads(audit.read_text())
    assert result["merged_rows"] == 2
    assert result["merged_missing_target"] == 2


def test_failed_lineage_leaves_no_outputs_or_temporary_files(tmp_path):
    taxdump = tmp_path / "taxdump"
    taxdump.mkdir()
    (taxdump / "nodes.dmp").write_text(
        dmp("1", "1", "no rank") + dmp("2759", "1", "superkingdom") +
        dmp("3", "999", "species"))
    (taxdump / "names.dmp").write_text(
        dmp("2759", "Eukaryota", "", "scientific name") +
        dmp("3", "Species alpha", "", "scientific name"))
    (taxdump / "merged.dmp").write_text("")
    (taxdump / "delnodes.dmp").write_text("")
    output = tmp_path / "paths.tsv"
    audit = tmp_path / "audit.json"
    report = tmp_path / "missing.tsv"
    with pytest.raises(ValueError, match="Missing parent taxid"):
        builder.main(["--taxdump", str(taxdump), "--output", str(output),
                      "--audit", str(audit), "--report", str(report)])
    assert not output.exists() and not audit.exists() and not report.exists()
    assert not list(tmp_path.glob(".*.tmp"))


@pytest.mark.parametrize("source_kind", ["taxdump", "ktaxonomy"])
def test_pre_domain_kraken_assignments_are_unclassified(tmp_path, source_kind):
    if source_kind == "taxdump":
        source = tmp_path / "taxdump"
        source.mkdir()
        (source / "nodes.dmp").write_text(
            dmp("1", "1", "no rank") + dmp("131567", "1", "no rank") +
            dmp("2759", "131567", "superkingdom"))
        (source / "names.dmp").write_text(
            dmp("2759", "Eukaryota", "", "scientific name"))
        (source / "merged.dmp").write_text("")
        (source / "delnodes.dmp").write_text("")
    else:
        source = tmp_path / "ktaxonomy.tsv"
        source.write_text(dmp("1", "1", "R", "0", "root") +
                          dmp("131567", "1", "R", "1", "cellular organisms") +
                          dmp("2759", "131567", "D", "2", "Eukaryota"))
    output = tmp_path / "paths.tsv"
    builder.main([f"--{source_kind}", str(source), "--output", str(output)])

    script = SCRIPT.with_name("convert_calls.py")
    spec = importlib.util.spec_from_file_location("convert_calls", script)
    converter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(converter)
    retained = tmp_path / "reads.txt"
    retained.write_text("r1\nr2\n")
    calls = tmp_path / "kraken.tsv"
    calls.write_text("C\tr1/1\t1\t100\nC\tr2/1\t131567\t100\n")
    converted = tmp_path / "converted.tsv"
    converter.main(["kraken2", "--read-ids", str(retained), "--input", str(calls),
                    "--taxid-paths", str(output), "--output", str(converted)])
    with converted.open(newline="") as handle:
        results = list(csv.DictReader(handle, delimiter="\t"))
    assert [(row["Taxonomy"], row["status"]) for row in results] == [
        ("", "unclassified"), ("", "unclassified")]


def test_ktaxonomy_duplicate_taxid_is_rejected(tmp_path):
    source = tmp_path / "ktaxonomy.tsv"
    source.write_text(dmp("1", "1", "R", "0", "root") +
                      dmp("1", "1", "R", "0", "different root"))
    output = tmp_path / "paths.tsv"
    with pytest.raises(ValueError, match="Duplicate taxid"):
        builder.main(["--ktaxonomy", str(source), "--output", str(output)])
    assert not output.exists()
