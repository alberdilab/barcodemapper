import csv
import gzip
import importlib.util
import json
from pathlib import Path

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "benchmarks" / "mjolnir" / "build_nt_accession_paths.py"
spec = importlib.util.spec_from_file_location("build_nt_accession_paths", SCRIPT)
builder = importlib.util.module_from_spec(spec)
spec.loader.exec_module(builder)


def hit(read, subject, taxids):
    return (f"{read}\t{subject}\t99\t100\t0\t0\t1\t100\t1\t100\t"
            f"1e-20\t200\t{taxids}\n")


def rows(path):
    with path.open(newline="") as handle:
        return {row["Accession"]: row["Taxonomy"] for row in csv.DictReader(handle, delimiter="\t")}


def test_multitaxid_repeated_subject_and_gzip(tmp_path):
    paths = tmp_path / "taxids.tsv.gz"
    with gzip.open(paths, "wt") as handle:
        handle.write("TaxID\tTaxonomy\n"
                     "10\tk__Animalia;p__Arthropoda;c__Insecta;o__;f__;g__;s__\n"
                     "11\tk__Animalia;p__Arthropoda;c__Arachnida;o__;f__;g__;s__\n"
                     "12\tk__Animalia;p__Mollusca;c__;o__;f__;g__;s__\n")
    blast = tmp_path / "blast.tsv.gz"
    with gzip.open(blast, "wt") as handle:
        handle.write(hit("r1/1", "acc|with|pipes", "10;11") +
                     hit("r1/2", "acc|with|pipes", "10") +
                     hit("r2/1", "other", "10") + hit("r2/2", "other", "12"))
    output = tmp_path / "accessions.tsv"
    audit = tmp_path / "audit.json"
    builder.main(["--blast", str(blast), "--taxid-paths", str(paths),
                  "--output", str(output), "--audit", str(audit)])
    assert rows(output) == {"acc|with|pipes": "k__Animalia;p__Arthropoda",
                            "other": "k__Animalia"}
    result = json.loads(audit.read_text())
    assert result["blast_rows"] == 4
    assert result["multi_taxid_rows"] == 1
    assert result["input_sha256"]["blast"] == builder.sha256(blast)


def test_unmapped_taxids_fail_or_omit_subjects(tmp_path):
    paths = tmp_path / "taxids.tsv"
    paths.write_text("TaxID\tTaxonomy\n10\tk__Animalia;p__Arthropoda\n")
    blast = tmp_path / "blast.tsv"
    blast.write_text(hit("r1/1", "good", "10") + hit("r2/1", "bad", "999") +
                     hit("r3/1", "none", "N/A") + hit("r4/1", "zero", "0") +
                     hit("r5/1", "mixed", "10") + hit("r5/2", "mixed", "999"))
    output = tmp_path / "accessions.tsv"
    audit = tmp_path / "audit.json"
    args = ["--blast", str(blast), "--taxid-paths", str(paths),
            "--output", str(output), "--audit", str(audit)]
    with pytest.raises(ValueError, match="4 BLAST subjects have unmapped taxids"):
        builder.main(args)
    assert not output.exists() and not audit.exists()
    builder.main([*args, "--allow-unmapped"])
    assert rows(output) == {"good": "k__Animalia;p__Arthropoda"}
    result = json.loads(audit.read_text())
    assert result["unmapped_subject_count"] == 4
    assert result["unmapped_subjects"] == {"bad": ["999"], "mixed": ["999"],
                                           "none": ["N/A"], "zero": ["0"]}
    assert result["no_taxid_rows"] == 2


def test_unmapped_as_empty_keeps_subjects_with_audited_reasons(tmp_path, capsys):
    paths = tmp_path / "taxids.tsv"
    paths.write_text("TaxID\tTaxonomy\n"
                     "10\tk__Animalia;p__Arthropoda\n"
                     "11\tk__Bacteria\n")
    blast = tmp_path / "blast.tsv"
    blast.write_text(hit("r1/1", "good", "10") +
                     hit("r2/1", "unknown", "999") +
                     hit("r3/1", "none", "N/A") +
                     hit("r4/1", "zero", "0") +
                     hit("r5/1", "mixed", "10") +
                     hit("r5/2", "mixed", "999") +
                     hit("r6/1", "no-common-kingdom", "10;11"))
    output = tmp_path / "accessions.tsv"
    audit = tmp_path / "audit.json"
    args = ["--blast", str(blast), "--taxid-paths", str(paths),
            "--output", str(output), "--audit", str(audit)]
    with pytest.raises(SystemExit) as error:
        builder.main([*args, "--allow-unmapped", "--unmapped-as-empty"])
    assert error.value.code == 2
    assert "not allowed with argument" in capsys.readouterr().err
    builder.main([*args, "--unmapped-as-empty"])
    assert rows(output) == {"good": "k__Animalia;p__Arthropoda", "unknown": "",
                            "none": "", "zero": "", "mixed": "",
                            "no-common-kingdom": ""}
    result = json.loads(audit.read_text())
    assert result["unmapped_as_empty_subject_count"] == 5
    assert result["unmapped_as_empty_subjects"] == {
        "unknown": ["taxid 999 absent from taxid paths"],
        "none": ["no taxid (N/A)"], "zero": ["no taxid (0)"],
        "mixed": ["taxid 999 absent from taxid paths"],
        "no-common-kingdom": ["no common kingdom"],
    }
    assert result["output_rows"] == 6


def test_refuses_to_overwrite_existing_output(tmp_path):
    paths = tmp_path / "taxids.tsv"
    paths.write_text("TaxID\tTaxonomy\n10\tk__Animalia\n")
    blast = tmp_path / "blast.tsv"
    blast.write_text(hit("r1/1", "acc", "10"))
    output = tmp_path / "accessions.tsv"
    output.write_text("keep")
    with pytest.raises(ValueError, match="already exists"):
        builder.main(["--blast", str(blast), "--taxid-paths", str(paths),
                      "--output", str(output)])
    assert output.read_text() == "keep"


def test_empty_taxid_paths_are_written_and_audited(tmp_path):
    paths = tmp_path / "taxids.tsv"
    paths.write_text("TaxID\tTaxonomy\n1\t\n131567\t\n10\tk__Animalia\n")
    blast = tmp_path / "blast.tsv"
    blast.write_text(hit("r1/1", "root", "1") +
                     hit("r2/1", "pre-domain", "131567") +
                     hit("r3/1", "mixed", "10") +
                     hit("r3/2", "mixed", "1") +
                     hit("r4/1", "good", "10") +
                     hit("r5/1", "both", "1;999"))
    output = tmp_path / "accessions.tsv"
    audit = tmp_path / "audit.json"
    builder.main(["--blast", str(blast), "--taxid-paths", str(paths),
                  "--output", str(output), "--audit", str(audit), "--allow-unmapped"])
    assert rows(output) == {"root": "", "pre-domain": "", "mixed": "",
                            "good": "k__Animalia"}
    result = json.loads(audit.read_text())
    assert result["empty_path_subjects"] == 4
    assert result["empty_path_subject_ids"] == ["both", "mixed", "pre-domain", "root"]
    assert result["unmapped_subject_count"] == 1


def test_empty_taxid_path_does_not_require_allow_unmapped(tmp_path):
    paths = tmp_path / "taxids.tsv"
    paths.write_text("TaxID\tTaxonomy\n1\t\n10\tk__Animalia\n")
    blast = tmp_path / "blast.tsv"
    blast.write_text(hit("r1/1", "root", "1") + hit("r2/1", "good", "10"))
    output = tmp_path / "accessions.tsv"
    builder.main(["--blast", str(blast), "--taxid-paths", str(paths),
                  "--output", str(output)])
    assert rows(output) == {"root": "", "good": "k__Animalia"}
    assert json.loads((tmp_path / "accessions.tsv.audit.json").read_text())["empty_path_subjects"] == 1
