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


def test_reference_validator_counts_both_markers(tmp_path):
    database = tmp_path / "valid.fa"
    database.write_text(
        ">coi|COI-5P|k__Animalia;p__Chordata\nACGT\n"
        ">its|ITS|k__Fungi;p__;c__Agaricomycetes\nACGT\n"
    )
    assert validate_database(database) == {"COI-5P": 1, "ITS": 1}


def test_malformed_source_headers_fail_instead_of_being_dropped(tmp_path):
    bold = tmp_path / "bold.fa"
    bold.write_text(">broken|COI-5P\nACGT\n")
    with pytest.raises(ValueError, match="Malformed BOLD"):
        dereplicate_bold(bold)

    unite = tmp_path / "unite.fa"
    unite.write_text(">broken|ITS\nACGT\n")
    with pytest.raises(ValueError, match="Malformed UNITE"):
        list(rename_unite(unite))
