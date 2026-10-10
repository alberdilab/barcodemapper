from barcodemapper.workflow.scripts.prepare_database import dereplicate_bold


def test_identical_sequence_with_conflicting_taxa_is_preserved(tmp_path):
    source = tmp_path / "bold.fa"
    source.write_text(
        ">A|COI-5P|x|Animalia,Chordata,Mammalia,Primates,Hominidae,Homo,sapiens\nACGTACGT\n"
        ">B|COI-5P|x|Animalia,Chordata,Mammalia,Primates,Hominidae,Pan,troglodytes\nACGTACGT\n"
    )
    records = dereplicate_bold(source)
    assert len(records) == 2
    assert {record.id for record in records} == {"A|COI-5P|x|Animalia,Chordata,Mammalia,Primates,Hominidae,Homo,sapiens", "B|COI-5P|x|Animalia,Chordata,Mammalia,Primates,Hominidae,Pan,troglodytes"}
