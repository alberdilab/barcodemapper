"""extract_taxonomy.py logic reading SAM text through pysam (isolates BAM write cost from parsing cost)."""
import argparse
import sys

import pysam

sys.path.insert(0, __import__("os").environ["PKG"])
from extract_taxonomy import parse_cigar  # noqa: E402

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--max-edit-percent", type=float, default=3.0)
parser.add_argument("--min-covered", type=int, default=80)
a = parser.parse_args()

with pysam.AlignmentFile(a.input, "r") as sam, open(a.output, "w") as out:
    for read in sam:
        if read.is_unmapped or not read.has_tag("NM"):
            continue
        nm = read.get_tag("NM")
        cov = parse_cigar(read.cigartuples)
        if cov < a.min_covered or nm * 100 > a.max_edit_percent * cov:
            continue
        p = sam.get_reference_name(read.reference_id).split("|")
        p += [""] * (3 - len(p))
        out.write(f"{read.query_name}\t{p[0]}\t{p[1]}\t{p[2]}\t{nm}\t{cov}\n")
