"""Keep only hits that pass BarcodeMapper --standard acceptance (covered >= 80
bases and NM*100 <= 3.0*covered), from vsearch userout or minimap2 PAF.
Usage: accept_hits.py {vsearch|paf} MATE(1|2) IN OUT
vsearch fields: query target id alnlen mism opens gaps ql qlo qhi tl caln qcov
Query IDs lacking a /1 or /2 suffix get /MATE appended. Prints counts."""
import re
import sys

MIN_COVERED = 80
MAX_EDIT_PERCENT = 3.0
OPS = re.compile(r"(\d*)([MIDNSHP=X])")


def ops(text):
    return [(int(n) if n else 1, op) for n, op in OPS.findall(text)]


def vsearch_stats(row):
    parsed = ops(row[11])
    covered = sum(n for n, op in parsed if op == "M")
    # terminal gap operations (first/last I or D) are overhang, not edits
    first_last = {0, len(parsed) - 1}
    gap_bases = sum(n for i, (n, op) in enumerate(parsed)
                    if op in "ID" and i not in first_last)
    return covered, int(row[4]) + gap_bases


def paf_stats(row):
    tags = {t[:2]: t[5:] for t in row[12:] if len(t) > 5 and t[2] == ":"}
    if "NM" not in tags or "cg" not in tags:
        raise ValueError("PAF row lacks NM or cg tag; run minimap2 with -c")
    covered = sum(n for n, op in ops(tags["cg"]) if op in "M=X")
    return covered, int(tags["NM"])


def main():
    fmt, mate, source, target = sys.argv[1:5]
    assert fmt in ("vsearch", "paf") and mate in ("1", "2")
    stats = vsearch_stats if fmt == "vsearch" else paf_stats
    total = kept = 0
    reads_in = set()
    reads_out = set()
    with open(source) as src, open(target, "x") as out:
        for line in src:
            row = line.rstrip("\n").split("\t")
            if not row or not row[0]:
                continue
            total += 1
            covered, nm = stats(row)
            query = row[0]
            reads_in.add(query)
            if covered >= MIN_COVERED and nm * 100 <= MAX_EDIT_PERCENT * covered:
                if not query.endswith(("/1", "/2")):
                    row[0] = f"{query}/{mate}"
                out.write("\t".join(row) + "\n")
                kept += 1
                reads_out.add(row[0])
    print(f"hits_in={total} hits_accepted={kept} mates_with_hit_in={len(reads_in)} mates_with_accepted_hit={len(reads_out)}")


main()
