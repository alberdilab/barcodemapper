import csv
import hashlib
import importlib.util
from pathlib import Path

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "benchmarks" / "mjolnir" / "evaluate_calls.py"
spec = importlib.util.spec_from_file_location("evaluate_calls", SCRIPT)
evaluate_calls = importlib.util.module_from_spec(spec)
spec.loader.exec_module(evaluate_calls)


def table(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run(tmp_path, calls, option, source, extra=()):
    call_args = []
    for name, content in calls.items():
        path = tmp_path / f"{name}.tsv"
        path.write_text(content)
        call_args += ["--calls", f"{name}={path}"]
    source_path = tmp_path / "source.tsv"
    source_path.write_text(source)
    out = tmp_path / "out"
    evaluate_calls.main([*call_args, option, str(source_path), "--out", str(out), *extra])
    return out


A = "k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Hypocreales;f__Nectriaceae;g__Fusarium;s__Fusarium_oxysporum"
B = A.rsplit(";", 1)[0] + ";s__Fusarium_solani"
C = "k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Agaricales;f__Agaricaceae;g__Agaricus;s__Agaricus_bisporus"


def test_hand_computed_synthetic_metrics_and_depth(tmp_path):
    truth = ("ReadID\tTaxonomy\toff_target\n"
             f"r1\t{A}\tno\n"
             f"r2\t{A}\tno\n"
             f"r3\t{A}\tno\n"
             "r4\t\tyes\n"
             "r5\t\tyes\n"
             f"r6\t{A}\tno\n")
    calls = ("ReadID\tTaxonomy\tstatus\n"
             f"r1\t{A}\tclassified\n"
             f"r2\t{B}\tclassified\n"
             f"r3\t{A.rsplit(';', 1)[0]}\tclassified\n"
             f"r4\t{C}\tclassified\n"
             "r5\t\tunclassified\n"
             "r6\t\tno_hit\n")
    out = run(tmp_path, {"method": calls}, "--truth", truth)
    species = next(row for row in table(out / "rank_metrics.tsv") if row["rank"] == "species")
    # Hand calculation: 1 correct / 3 species calls; 1 correct / 4 target pairs; 1 of 2 off-target pairs called.
    assert (species["TP"], species["FP"], species["FN"], species["precision"], species["recall"], species["wrong_species_overcalls"], species["false_calls_per_million_off_target_pairs"]) == ("1", "2", "3", "0.333333", "0.250000", "1", "500000.000000")
    phylum = next(row for row in table(out / "rank_metrics.tsv") if row["rank"] == "phylum")
    assert (phylum["TP"], phylum["FP"], phylum["FN"]) == ("3", "1", "1")
    depths = {row["depth"]: int(row["pairs"]) for row in table(out / "resolution_depth.tsv")}
    assert depths["species"] == 3
    assert depths["genus"] == depths["unclassified"] == depths["no_hit"] == 1


def test_undefined_denominators_are_blank(tmp_path):
    out = run(tmp_path, {"m": "ReadID\tTaxonomy\tstatus\nr1\t\tno_hit\n"},
              "--truth", f"ReadID\tTaxonomy\toff_target\nr1\t{A}\tno\n")
    species = next(row for row in table(out / "rank_metrics.tsv") if row["rank"] == "species")
    assert species["precision"] == ""
    assert species["recall"] == "0.000000"
    assert species["false_calls_per_million_off_target_pairs"] == ""


def test_species_call_against_genus_truth_is_overcall(tmp_path):
    truth = f"ReadID\tTaxonomy\toff_target\nr1\t{A.rsplit(';', 1)[0]}\tno\n"
    calls = f"ReadID\tTaxonomy\tstatus\nr1\t{A}\tclassified\n"
    out = run(tmp_path, {"m": calls}, "--truth", truth)
    species = next(row for row in table(out / "rank_metrics.tsv") if row["rank"] == "species")
    assert (species["TP"], species["FP"], species["FN"],
            species["wrong_species_overcalls"], species["unscored_calls"]) == (
                "0", "1", "0", "1", "0")


def test_mock_threshold_keeps_reference_absent_member(tmp_path):
    expected = ("run_accession\texpected_species\tsource_families\tcoverage_status\t"
                "exact_species_reference_records\tgenus_reference_records\tfamily_reference_records\n"
                "run1\tFusarium oxysporum\tNectriaceae\texact_species\t3\t5\t8\n"
                "run1\tAgaricus bisporus\tAgaricaceae\tabsent_from_reference\t0\t0\t0\n"
                "run2\tOther species\tOtheraceae\texact_species\t1\t1\t1\n")
    calls = ("ReadID\tTaxonomy\tstatus\n"
             f"r1\t{A}\tclassified\n"
             f"r2\t{A}\tclassified\n"
             f"r3\t{B}\tclassified\n"
             "r4\t\tno_hit\n")
    out = run(tmp_path, {"m": calls}, "--expected-taxa", expected,
              ["--run-accession", "run1", "--min-fragments", "2"])
    detection = {row["expected_species"]: row for row in table(out / "expected_detection.tsv")}
    assert len(detection) == 2
    assert (detection["Fusarium oxysporum"]["fragments"],
            detection["Fusarium oxysporum"]["detected"]) == ("2", "yes")
    assert (detection["Agaricus bisporus"]["coverage_status"],
            detection["Agaricus bisporus"]["fragments"],
            detection["Agaricus bisporus"]["detected"]) == ("absent_from_reference", "0", "no")
    species = next(row for row in table(out / "mock_metrics.tsv") if row["rank"] == "species")
    assert (species["expected_taxa"], species["detected_taxa"], species["unexpected_taxa"]) == ("2", "1", "0")
    assert table(out / "unexpected_calls.tsv") == [{"method": "m", "rank": "species",
                                                  "taxon": "Fusarium_solani", "fragments": "1",
                                                  "passes_threshold": "no"}]


def test_mock_species_calls_support_genus_and_family(tmp_path):
    expected = ("run_accession\texpected_species\tsource_families\tcoverage_status\t"
                "exact_species_reference_records\tgenus_reference_records\tfamily_reference_records\n"
                "run1\tFusarium oxysporum\tNectriaceae\texact_species\t3\t5\t8\n")
    calls = ("ReadID\tTaxonomy\tstatus\n"
             f"r1\t{A}\tclassified\n"
             f"r2\t{A}\tclassified\n")
    out = run(tmp_path, {"m": calls}, "--expected-taxa", expected,
              ["--run-accession", "run1", "--min-fragments", "2"])
    metrics = {row["rank"]: row for row in table(out / "mock_metrics.tsv")}
    assert {rank: metrics[rank]["detected_taxa"] for rank in ("species", "genus", "family")} == {
        "species": "1", "genus": "1", "family": "1"}


@pytest.mark.parametrize("changed_field", range(3))
def test_hash_mismatch_refuses_before_writing_outputs(tmp_path, capsys, changed_field):
    calls = "ReadID\tTaxonomy\tstatus\nr1\t\tno_hit\n"
    truth = f"ReadID\tTaxonomy\toff_target\nr1\t{A}\tno\n"
    hashes = tmp_path / "hashes.tsv"
    first = ["a" * 64, "b" * 64, "c" * 64]
    second = first.copy()
    second[changed_field] = "d" * 64
    hashes.write_text("method\tfastq_r1_sha256\tfastq_r2_sha256\treference_sha256\n"
                      + "one\t" + "\t".join(first) + "\n"
                      + "two\t" + "\t".join(second) + "\n")
    with pytest.raises(SystemExit) as error:
        run(tmp_path, {"one": calls, "two": calls}, "--truth", truth,
            ["--hashes-manifest", str(hashes), "--allow-unverified-hashes"])
    assert error.value.code != 0
    assert "FASTQ/reference sha256 differ" in capsys.readouterr().err
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("missing", ["columns", "cell"])
def test_hash_manifest_requires_file_paths_unless_opted_out(tmp_path, capsys, missing):
    calls = "ReadID\tTaxonomy\tstatus\nr1\t\tno_hit\n"
    truth = f"ReadID\tTaxonomy\toff_target\nr1\t{A}\tno\n"
    hashes = tmp_path / "hashes.tsv"
    digest = hashlib.sha256(b"sample").hexdigest()
    header = "method\tfastq_r1_sha256\tfastq_r2_sha256\treference_sha256"
    row = f"m\t{digest}\t{digest}\t{digest}"
    if missing == "cell":
        header += "\tfastq_r1_path\tfastq_r2_path\treference_path"
        row += "\t  \tr2.fq\treference.fa"
        (tmp_path / "r2.fq").write_bytes(b"sample")
        (tmp_path / "reference.fa").write_bytes(b"sample")
    hashes.write_text(f"{header}\n{row}\n")
    manifest = ["--hashes-manifest", str(hashes)]
    with pytest.raises(SystemExit) as error:
        run(tmp_path, {"m": calls}, "--truth", truth, manifest)
    assert error.value.code != 0
    assert "fastq_r1_path" in capsys.readouterr().err
    assert not (tmp_path / "out").exists()

    run(tmp_path, {"m": calls}, "--truth", truth,
        [*manifest, "--allow-unverified-hashes"])
    assert "unverified" in capsys.readouterr().err.lower()
    assert (tmp_path / "out" / "rank_metrics.tsv").exists()


def test_hash_manifest_verifies_listed_files(tmp_path, capsys):
    calls = "ReadID\tTaxonomy\tstatus\nr1\t\tno_hit\n"
    truth = f"ReadID\tTaxonomy\toff_target\nr1\t{A}\tno\n"
    paths = [tmp_path / name for name in ("r1.fq", "r2.fq", "reference.fa")]
    for path, content in zip(paths, (b"mate 1", b"mate 2", b"reference")):
        path.write_bytes(content)
    hashes = tmp_path / "hashes.tsv"
    hashes.write_text("method\tfastq_r1_sha256\tfastq_r2_sha256\treference_sha256\t"
                      "fastq_r1_path\tfastq_r2_path\treference_path\n"
                      "m\t" + "\t".join(hashlib.sha256(path.read_bytes()).hexdigest()
                                         for path in paths) + "\t" +
                      "\t".join(path.name for path in paths) + "\n")
    run(tmp_path, {"m": calls}, "--truth", truth,
        ["--hashes-manifest", str(hashes)])
    assert (tmp_path / "out" / "rank_metrics.tsv").exists()
    paths[0].write_bytes(b"changed")
    with pytest.raises(SystemExit):
        run(tmp_path, {"m": calls}, "--truth", truth,
            ["--hashes-manifest", str(hashes)])
    assert "sha256 mismatch" in capsys.readouterr().err


def test_differing_read_ids_are_rejected(tmp_path, capsys):
    calls = "ReadID\tTaxonomy\tstatus\nr1\t\tno_hit\n"
    with pytest.raises(SystemExit):
        run(tmp_path, {"one": calls,
                       "two": "ReadID\tTaxonomy\tstatus\nr2\t\tno_hit\n"},
            "--truth", f"ReadID\tTaxonomy\toff_target\nr1\t{A}\tno\n")
    assert "differing ReadID sets" in capsys.readouterr().err
