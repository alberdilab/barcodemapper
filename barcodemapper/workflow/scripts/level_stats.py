"""Summary data for the multi-level HTML report.

``overview_stats(merged_by_level, extract_stats)`` returns a JSON-serialisable
dictionary with this exact shape::

    {
      "levels": {
        level_name: {
          "samples": {sample_name: summary, ...},
          "totals": summary
        }, ...
      },
      "classified_by_sample": [
        {"level": level_name, "sample": sample_name,
         "classified_fragments": integer}, ...
      ]
    }

Level order follows ``merged_by_level``. Sample order follows the TSV columns,
then any additional names in ``extract_stats``. The comparison rows are in
level-major, sample-minor order. Each ``summary`` contains:

* ``classified_fragments`` and ``unclassified_fragments``: integers from the
  depth-0 (kingdom and Unclassified) rows, avoiding cumulative double counts.
* ``depth_distribution``: seven integers indexed by rank depth 0..6 (kingdom,
  phylum, class, order, family, genus, species); each counts fragments resolved
  to at least that rank. Skipped intermediate rank names do not remove a
  fragment from shallower thresholds.
* ``distinct_taxa``: seven integers for taxonomy prefixes at each rank whose
  count is positive. Totals count the union across samples.
* ``kingdom_composition``: a mapping of kingdom prefix to its positive depth-0
  count. Totals sum counts across samples.
* ``fragments_retained``, ``aligned_fragments``, ``retained_fraction``: numbers
  from extractor stats. Missing values and fractions with a zero aligned
  denominator are ``None``. Totals sum available stats values; their fraction
  uses only samples with both retained and aligned counts.

All integer values are Python ints. An empty level has zero counts, empty
kingdom composition, and still includes every known sample.
"""

import json

import pandas as pd


RANKS = "kpcofgs"


def load_merged(path):
    """Read one merged level TSV, including headers when it has no rows."""
    return pd.read_csv(path, sep="\t")


def load_stats(paths):
    """Read per-sample extractor stats JSON files in path order."""
    loaded = []
    for path in paths:
        with open(path, encoding="utf-8") as handle:
            loaded.append(json.load(handle))
    return loaded


def _fraction(retained, aligned):
    return retained / aligned if retained is not None and aligned is not None and aligned > 0 else None


def overview_stats(merged_by_level: dict[str, pd.DataFrame], extract_stats: list[dict]) -> dict:
    """Build per-sample and per-level overview data from cumulative TSVs."""
    sample_names = []
    for frame in merged_by_level.values():
        for sample in frame.columns:
            if sample != "Taxonomy" and sample not in sample_names:
                sample_names.append(sample)
    stats_by_sample = {}
    for entry in extract_stats:
        sample = entry.get("sample")
        if sample is not None:
            stats_by_sample[sample] = entry
            if sample not in sample_names:
                sample_names.append(sample)

    levels = {}
    comparison = []
    for level_name, frame in merged_by_level.items():
        taxonomy = frame["Taxonomy"].fillna("").astype(str)
        classified = taxonomy.str.startswith("k__")
        names = taxonomy[classified]
        ranks = names.str.rsplit(";", n=1).str[-1].str[0].map({rank: depth for depth, rank in enumerate(RANKS)})
        parents = names.str.rpartition(";", expand=False).str[0]
        counts = frame.reindex(columns=sample_names, fill_value=0).fillna(0).astype(int)
        samples = {}

        for sample in sample_names:
            values = counts.loc[classified, sample]
            child_counts = values.groupby(parents).sum()
            terminal = values - names.map(child_counts).fillna(0).astype(int)
            depth_distribution = [int(terminal[ranks.ge(depth)].sum()) for depth in range(len(RANKS))]
            distinct_taxa = [int((values.gt(0) & ranks.eq(depth)).sum()) for depth in range(len(RANKS))]
            kingdoms = {name: int(count) for name, count in zip(names[ranks.eq(0)], values[ranks.eq(0)]) if count > 0}
            unclassified = int(counts.loc[taxonomy.eq("Unclassified"), sample].sum())

            entry = stats_by_sample.get(sample, {})
            aligned = entry.get("aligned_fragments")
            retained = entry.get("levels", {}).get(level_name, {}).get("fragments")
            summary = {
                "classified_fragments": int(values[ranks.eq(0)].sum()),
                "unclassified_fragments": unclassified,
                "depth_distribution": depth_distribution,
                "distinct_taxa": distinct_taxa,
                "kingdom_composition": kingdoms,
                "fragments_retained": int(retained) if retained is not None else None,
                "aligned_fragments": int(aligned) if aligned is not None else None,
                "retained_fraction": _fraction(retained, aligned),
            }
            samples[sample] = summary
            comparison.append({"level": level_name, "sample": sample,
                               "classified_fragments": summary["classified_fragments"]})

        retained_pairs = [(item["fragments_retained"], item["aligned_fragments"])
                          for item in samples.values()
                          if item["fragments_retained"] is not None and item["aligned_fragments"] is not None]
        retained_values = [item["fragments_retained"] for item in samples.values()
                           if item["fragments_retained"] is not None]
        aligned_values = [item["aligned_fragments"] for item in samples.values()
                          if item["aligned_fragments"] is not None]
        all_values = counts.loc[classified]
        present = all_values.gt(0).any(axis=1)
        kingdom_composition = {}
        for summary in samples.values():
            for kingdom, value in summary["kingdom_composition"].items():
                kingdom_composition[kingdom] = kingdom_composition.get(kingdom, 0) + value
        totals = {
            "classified_fragments": sum(item["classified_fragments"] for item in samples.values()),
            "unclassified_fragments": sum(item["unclassified_fragments"] for item in samples.values()),
            "depth_distribution": [sum(item["depth_distribution"][depth] for item in samples.values())
                                   for depth in range(len(RANKS))],
            "distinct_taxa": [int((present & ranks.eq(depth)).sum()) for depth in range(len(RANKS))],
            "kingdom_composition": kingdom_composition,
            "fragments_retained": sum(retained_values) if retained_values else None,
            "aligned_fragments": sum(aligned_values) if aligned_values else None,
            "retained_fraction": _fraction(sum(pair[0] for pair in retained_pairs),
                                            sum(pair[1] for pair in retained_pairs)) if retained_pairs else None,
        }
        levels[level_name] = {"samples": samples, "totals": totals}

    return {"levels": levels, "classified_by_sample": comparison}
