#!/usr/bin/env python3
"""Score the hit-cap extension runs locally (no network, no HPC).

Reads hit-cap/data/ (result tables and per-read taxonomy pulled from Mjolnir) and writes
extension_compare.tsv (capped vs all-hits assignments per read) and extension_accuracy.tsv
(precision/recall against the expected taxa, same scoring as the pilot).

Rerun:  python3 benchmarks/hit-cap/analyze_extension.py
"""
import glob
import importlib.util
import os

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
BENCH = os.path.dirname(HERE)
DATA = os.path.join(HERE, "data")


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


compare = load("compare_max_hits", os.path.join(BENCH, "mjolnir", "compare_max_hits.py"))
pilot = load("evaluate_pilot", os.path.join(BENCH, "pilot-01", "evaluate_pilot.py"))

# folder -> run accession, marker, reference
SAMPLES = {
    "animal": ("ERR6001161", "COI", "default", "Animalia"),
    "ERR6001164": ("ERR6001164", "COI", "default", "Animalia"),
    "ERR6001167": ("ERR6001167", "COI", "default", "Animalia"),
    "ERR6001170": ("ERR6001170", "COI", "default", "Animalia"),
    "fungal": ("SRR14351851", "ITS", "all-fungi", "Fungi"),
    "SRR14351849": ("SRR14351849", "ITS", "all-fungi", "Fungi"),
}


def tag_key(tag):
    return -1 if tag == "all" else int(tag[1:])


def main():
    expected_taxa = pd.read_csv(os.path.join(BENCH, "pilot-01", "inputs", "expected-taxa.tsv"), sep="\t")
    compare_rows, metric_rows = [], []
    for folder, (run, marker, reference, kingdom) in SAMPLES.items():
        tags = sorted((os.path.basename(p)[:-4] for p in glob.glob(os.path.join(DATA, "results", folder, "*.tsv"))),
                      key=tag_key)
        taxonomy = lambda tag: os.path.join(DATA, "work", folder, tag, "taxonomy", f"{run}.tsv")
        baseline = compare.read_assignments(taxonomy("all"))
        for tag in tags:
            row = {"sample": folder, "run": tag, **compare.compare(baseline, compare.read_assignments(taxonomy(tag)))}
            compare_rows.append(row)
        expected, _ = pilot.build_expected(expected_taxa[expected_taxa.run_accession == run])
        meta = {"marker_label": marker, "reference": reference}
        for tag in tags:
            tax, _ = pilot.parse_taxonomy(os.path.join(DATA, "results", folder, f"{tag}.tsv"), run)
            for metric in pilot.rank_metrics(run, meta, tax, expected, ["all", kingdom]):
                metric.update(sample=folder, setting=tag)
                metric_rows.append(metric)
    pd.DataFrame(compare_rows).to_csv(os.path.join(HERE, "extension_compare.tsv"), sep="\t", index=False)
    pd.DataFrame(metric_rows).to_csv(os.path.join(HERE, "extension_accuracy.tsv"), sep="\t", index=False,
                                     float_format="%.4f")


if __name__ == "__main__":
    main()
