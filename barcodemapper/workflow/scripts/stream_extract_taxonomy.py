"""Stream SAM text on stdin and write the same TSV as extract_taxonomy.py (no BAM, no pysam)."""
import argparse
import re
import subprocess
import sys

CIGAR_OP = re.compile(r"(\d+)([MIDNSHP=X])")
COVERED = {"M", "=", "X"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--max-edit-percent", type=float, default=3.0)
    parser.add_argument("--min-covered", type=int, default=80)
    parser.add_argument("--bam-out", help="Also write the incoming SAM to this BAM (needs samtools)")
    parser.add_argument("--threads", type=int, default=1, help="samtools compression threads")
    args = parser.parse_args()
    max_edit, min_cov = args.max_edit_percent, args.min_covered

    bam = None
    if args.bam_out:
        bam = subprocess.Popen(
            ["samtools", "view", "-b", "-@", str(args.threads), "-o", args.bam_out, "-"],
            stdin=subprocess.PIPE, text=True)
    try:
        extract(sys.stdin, args.output, max_edit, min_cov, bam.stdin if bam else None)
    except BaseException:
        if bam:
            bam.kill()
        raise
    if bam:
        bam.stdin.close()
        if bam.wait() != 0:
            sys.exit(f"samtools failed writing {args.bam_out}")


def extract(lines, output, max_edit, min_cov, passthrough=None):
    with open(output, "w", buffering=1 << 20) as out:
        write = out.write
        for line in lines:
            if passthrough is not None:
                passthrough.write(line)
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
