import gzip
import subprocess
import sys
from pathlib import Path

from barcodemapper.workflow.scripts.prepare_database import dereplicate_bold, rename_unite
from barcodemapper.workflow.scripts.validate_database import validate_database

import pytest


def test_identical_sequence_with_conflicting_taxa_is_preserved(tmp_path):
    source = tmp_path / "bold.fa"
    source.write_text(
        ">A|COI-5P|x|Animalia,Chordata,Mammalia,Primates,Hominidae,Homo,sapiens\nACGTACGT\n"
        ">B|COI-5P|x|Animalia,Chordata,Mammalia,Primates,Hominidae,Pan,troglodytes\nACGTACGT\n"
    )
    records = dereplicate_bold(source)
    assert len(records) == 2
    assert {record.id for record in records} == {"A|COI-5P|x|Animalia,Chordata,Mammalia,Primates,Hominidae,Homo,sapiens", "B|COI-5P|x|Animalia,Chordata,Mammalia,Primates,Hominidae,Pan,troglodytes"}


@pytest.mark.parametrize("header", [
    "ref|ITS",
    "ref|unknown|k__Fungi",
    "ref|ITS|p__Basidiomycota",
    "ref|ITS|k__Fungi;c__Agaricomycetes;p__Basidiomycota",
    "ref|ITS|k__Fungi with spaces",
])
def test_invalid_reference_headers_are_rejected(tmp_path, header):
    database = tmp_path / "invalid.fa"
    database.write_text(f">{header}\nACGT\n")
    with pytest.raises(ValueError):
        validate_database(database)


@pytest.mark.parametrize("compressed", [False, True])
def test_reference_validator_counts_both_markers(tmp_path, compressed):
    database = tmp_path / ("valid.fa.gz" if compressed else "valid.fa")
    fasta = (
        ">coi|COI-5P|k__Animalia;p__Chordata\nACGT\n"
        ">its|ITS|k__Fungi;p__;c__Agaricomycetes\nACGT\n"
    )
    if compressed:
        with gzip.open(database, "wt") as handle:
            handle.write(fasta)
    else:
        database.write_text(fasta)
    assert validate_database(database) == {"COI-5P": 1, "ITS": 1}


def test_prepare_database_writes_gzipped_fasta(tmp_path):
    bold = tmp_path / "bold.fa"
    bold.write_text(">B1|COI-5P|x|Animalia,Chordata\nACGTACGT\n")
    unite = tmp_path / "unite.fa"
    unite.write_text(">record|U1|x|x|k__Fungi;p__Basidiomycota\nTGCA\n")
    database = tmp_path / "combined.fa.gz"
    script = Path(__file__).resolve().parents[1] / "barcodemapper/workflow/scripts/prepare_database.py"
    subprocess.run([sys.executable, str(script), "-b", str(bold), "-u", str(unite),
                    "-y", "k__Fungi", "-o", str(database)], check=True, capture_output=True)

    assert database.read_bytes().startswith(b"\x1f\x8b")
    assert validate_database(database) == {"COI-5P": 1, "ITS": 1}
    assert not (tmp_path / "combined.fa").exists()


def test_malformed_source_headers_fail_instead_of_being_dropped(tmp_path):
    bold = tmp_path / "bold.fa"
    bold.write_text(">broken|COI-5P\nACGT\n")
    with pytest.raises(ValueError, match="Malformed BOLD"):
        dereplicate_bold(bold)

    unite = tmp_path / "unite.fa"
    unite.write_text(">broken|ITS\nACGT\n")
    with pytest.raises(ValueError, match="Malformed UNITE"):
        list(rename_unite(unite))
