"""Checks for extracting several stringency levels from one SAM stream."""

import importlib.util
from io import StringIO
import json
import subprocess
import sys
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "barcodemapper/workflow/scripts/stream_extract_taxonomy.py"
PRESETS = {
    "ultrastrict": (0.0, 120),
    "strict": (2.0, 100),
    "standard": (3.0, 80),
    "permissive": (5.0, 60),
}


def sam_input(tmp_path):
    rows = [
        ("r0", 0, "120M", "NM:i:0"),
        ("r1", 0, "110M", "NM:i:0"),
        ("r2", 0, "100M", "NM:i:2"),
        ("r3", 0, "80M", "NM:i:2"),
        ("r4", 0, "60M", "NM:i:3"),
        ("r5", 0, "120M", "NM:i:6"),
        ("r6", 4, "*", "NM:i:0"),
        ("r7", 0, "*", "NM:i:0"),
        ("r8", 0, "120M", "AS:i:0"),
        ("r0", 256, "110M", "NM:i:0"),
        ("r2", 0, "120M", "NM:i:4"),
    ]
    sam = tmp_path / "reads.sam"
    sam.write_text("@HD\tVN:1.6\n@SQ\tSN:ref|ITS|k__Fungi\tLN:1000\n" + "".join(
        "\t".join([name, str(flag), "*" if flag & 4 else "ref|ITS|k__Fungi",
                   "0" if flag & 4 else "1", "60", cigar, "*", "0", "0",
                   "A" * 120, "I" * 120, tag]) + "\n"
        for name, flag, cigar, tag in rows
    ))
    return sam


def run_extract(sam, *args):
    with sam.open() as handle:
        subprocess.run([sys.executable, str(SCRIPT), *map(str, args)], stdin=handle, check=True)


def test_legacy_output_and_stats(tmp_path):
    sam = sam_input(tmp_path)
    output = tmp_path / "legacy.tsv"
    stats = tmp_path / "legacy.stats.json"
    run_extract(sam, "--output", output, "--max-edit-percent", 2, "--min-covered", 100,
                "--stats-json", stats, "--sample", "sample.one")
    assert output.read_text() == "".join(
        f"{name}\tref\tITS\tk__Fungi\t{nm}\t{cov}\n"
        for name, nm, cov in [("r0", 0, 120), ("r1", 0, 110), ("r2", 2, 100), ("r0", 0, 110)]
    )
    assert json.loads(stats.read_text()) == {
        "sample": "sample.one", "records": 11, "aligned_fragments": 7,
        "levels": {"custom": {"alignments": 4, "fragments": 3}},
    }


def test_four_levels_match_separate_runs_and_stats(tmp_path):
    sam = sam_input(tmp_path)
    pattern = tmp_path / "extract/{level}/sample.tsv"
    stats = tmp_path / "extract/sample.stats.json"
    level_args = [arg for name, (edit, cov) in PRESETS.items()
                  for arg in ("--level", f"{name}:{edit}:{cov}")]
    run_extract(sam, *level_args, "--output-pattern", pattern, "--stats-json", stats,
                "--sample", "sample.one")
    for name, (edit, cov) in PRESETS.items():
        separate = tmp_path / f"{name}.tsv"
        run_extract(sam, "--output", separate, "--max-edit-percent", edit, "--min-covered", cov)
        assert Path(str(pattern).replace("{level}", name)).read_bytes() == separate.read_bytes()

    names = {name: [line.split("\t", 1)[0] for line in
                    Path(str(pattern).replace("{level}", name)).read_text().splitlines()]
             for name in PRESETS}
    assert names["ultrastrict"] == ["r0"]
    assert "r1" in names["strict"] and "r1" not in names["ultrastrict"]
    assert all(set(names[earlier]) <= set(names[later]) for earlier, later in
               zip(list(PRESETS)[:-1], list(PRESETS)[1:]))
    assert json.loads(stats.read_text()) == {
        "sample": "sample.one", "records": 11, "aligned_fragments": 7,
        "levels": {
            "ultrastrict": {"alignments": 1, "fragments": 1},
            "strict": {"alignments": 4, "fragments": 3},
            "standard": {"alignments": 5, "fragments": 4},
            "permissive": {"alignments": 8, "fragments": 6},
        },
    }


def test_non_nested_levels_are_independent(tmp_path):
    sam = sam_input(tmp_path)
    pattern = tmp_path / "{level}.tsv"
    run_extract(sam, "--level", "short_clean:0:60", "--level", "long_loose:5:120",
                "--output-pattern", pattern)
    short = (tmp_path / "short_clean.tsv").read_text()
    long = (tmp_path / "long_loose.tsv").read_text()
    assert "r1\t" in short and "r1\t" not in long
    assert "r5\t" in long and "r5\t" not in short


def test_passthrough_receives_headers_and_filtered_records(tmp_path):
    sam = sam_input(tmp_path)
    spec = importlib.util.spec_from_file_location("stream_extract_taxonomy", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    passthrough = StringIO()
    with sam.open() as source:
        module.extract_levels(source, str(tmp_path / "{level}.tsv"), PRESETS,
                              passthrough=passthrough)
    assert passthrough.getvalue() == sam.read_text()
