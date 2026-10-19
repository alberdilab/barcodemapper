"""Write the Accession<TAB>Taxonomy map convert_calls.py needs, for the
subjects in HITS (column 2 of vsearch rows, column 6 of PAF rows).
Usage: zcat REF.fa.gz | grep '^>' | make_hit_map.py {vsearch|paf} HITS OUT
Header: >accession|marker|taxonomy_path ; key is the full header label."""
import sys

fmt, hits, out = sys.argv[1:4]
column = 1 if fmt == "vsearch" else 5
subjects = set()
with open(hits) as handle:
    for line in handle:
        row = line.rstrip("\n").split("\t")
        if len(row) > column:
            subjects.add(row[column])
found = {}
for line in sys.stdin:
    label = line[1:].strip()
    if label in subjects:
        parts = label.split("|", 2)
        if len(parts) != 3 or not parts[2]:
            raise SystemExit(f"Unexpected reference header {label!r}")
        found[label] = parts[2]
missing = subjects - found.keys()
if missing:
    raise SystemExit(f"{len(missing)} subjects not in reference, e.g. {sorted(missing)[0]!r}")
with open(out, "x") as handle:
    handle.write("Accession\tTaxonomy\n")
    for label in sorted(found):
        handle.write(f"{label}\t{found[label]}\n")
print(f"subjects={len(subjects)} mapped={len(found)}")
