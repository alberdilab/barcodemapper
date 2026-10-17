"""Stream SAM text on stdin and write the same TSV as extract_taxonomy.py (no BAM, no pysam)."""
import argparse
import re
import sys

CIGAR_OP = re.compile(r"(\d+)([MIDNSHP=X])")
COVERED = {"M", "=", "X"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--max-edit-percent", type=float, default=3.0)
    parser.add_argument("--min-covered", type=int, default=80)
    args = parser.parse_args()
    max_edit, min_cov = args.max_edit_percent, args.min_covered

    with open(args.output, "w", buffering=1 << 20) as out:
        write = out.write
        for line in sys.stdin:
            if line[0] == "@":
                continue
            f = line.split("\t", 11)
            if int(f[1]) & 4:
                continue
            cigar = f[5]
            if cigar == "*":
                continue
            covered = 0
            for length, op in CIGAR_OP.findall(cigar):
                if op in COVERED:
                    covered += int(length)
            if covered < min_cov:
                continue
            nm = None
            for tag in f[11].split("\t") if len(f) > 11 else ():
                if tag.startswith("NM:i:"):
                    nm = int(tag[5:])
                    break
            if nm is None or nm * 100 > max_edit * covered:
                continue
            parts = f[2].split("|")
            ref = parts[0]
            marker = parts[1] if len(parts) > 1 else ""
            tax = parts[2] if len(parts) > 2 else ""
            write(f"{f[0]}\t{ref}\t{marker}\t{tax}\t{nm}\t{covered}\n")


if __name__ == "__main__":
    main()
