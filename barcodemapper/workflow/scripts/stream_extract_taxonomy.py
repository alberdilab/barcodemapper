"""Stream SAM text on stdin and write the same TSV as extract_taxonomy.py (no BAM, no pysam)."""
import argparse
from contextlib import ExitStack
import json
from pathlib import Path
import re
import subprocess
import sys

CIGAR_OP = re.compile(r"(\d+)([MIDNSHP=X])")
COVERED = {"M", "=", "X"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output")
    parser.add_argument("--max-edit-percent", type=float, default=3.0)
    parser.add_argument("--min-covered", type=int, default=80)
    parser.add_argument("--level", action="append", default=[], metavar="NAME:MAXPCT:MINCOV")
    parser.add_argument("--output-pattern", help="Output TSV path containing {level}")
    parser.add_argument("--stats-json", help="Write per-sample extraction counts as JSON")
    parser.add_argument("--sample", help="Sample name recorded in stats JSON")
    parser.add_argument("--bam-out", help="Also write the incoming SAM to this BAM (needs samtools)")
    parser.add_argument("--threads", type=int, default=1, help="samtools compression threads")
    args = parser.parse_args()
    max_edit, min_cov = args.max_edit_percent, args.min_covered
    levels = {}
    if args.level:
        if args.output or not args.output_pattern or "{level}" not in args.output_pattern:
            parser.error("--level requires --output-pattern containing {level}, without --output")
        for value in args.level:
            parts = value.split(":")
            if len(parts) != 3 or not parts[0] or parts[0] == "overview" or parts[0] in levels:
                parser.error(f"invalid or duplicate --level: {value}")
            try:
                levels[parts[0]] = (float(parts[1]), int(parts[2]))
            except ValueError:
                parser.error(f"invalid --level thresholds: {value}")
    elif not args.output or args.output_pattern:
        parser.error("--output is required without --level; --output-pattern requires --level")

    bam = None
    if args.bam_out:
        bam = subprocess.Popen(
            ["samtools", "view", "-b", "-@", str(args.threads), "-o", args.bam_out, "-"],
            stdin=subprocess.PIPE, text=True)
    try:
        if levels:
            stats = extract_levels(sys.stdin, args.output_pattern, levels,
                                   bam.stdin if bam else None, args.stats_json is not None, args.sample)
        else:
            stats = extract(sys.stdin, args.output, max_edit, min_cov,
                            bam.stdin if bam else None, args.stats_json is not None, args.sample)
    except BaseException:
        if bam:
            bam.kill()
        raise
    if bam:
        bam.stdin.close()
        if bam.wait() != 0:
            sys.exit(f"samtools failed writing {args.bam_out}")
    if args.stats_json:
        path = Path(args.stats_json)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w") as out:
            json.dump(stats, out)
            out.write("\n")


def extract(lines, output, max_edit, min_cov, passthrough=None, collect_stats=False, sample=None):
    with open(output, "w", buffering=1 << 20) as out:
        return _extract(lines, {"custom": out}, {"custom": (max_edit, min_cov)},
                        passthrough, collect_stats, sample)


def extract_levels(lines, output_pattern, levels, passthrough=None, collect_stats=False, sample=None):
    with ExitStack() as stack:
        outputs = {}
        for name in levels:
            path = Path(output_pattern.replace("{level}", name))
            path.parent.mkdir(parents=True, exist_ok=True)
            outputs[name] = stack.enter_context(path.open("w", buffering=1 << 20))
        return _extract(lines, outputs, levels, passthrough, collect_stats, sample)


def _extract(lines, outputs, levels, passthrough, collect_stats, sample):
    records = 0
    aligned_fragments = 0
    kept_fragments = {name: 0 for name in levels} if collect_stats else None
    # Bowtie2 emits a fragment's alignments contiguously, so only the last read ID
    # is needed to count distinct fragments without storing whole-sample ID sets.
    last_aligned = None
    last_kept = {name: None for name in levels} if collect_stats else None
    alignments = {name: 0 for name in levels} if collect_stats else None
    for line in lines:
        if passthrough is not None:
            passthrough.write(line)
        if line[0] == "@":
            continue
        if collect_stats:
            records += 1
        f = line.split("\t", 11)
        if int(f[1]) & 4:
            continue
        cigar = f[5]
        if cigar == "*":
            continue
        if collect_stats:
            if f[0] != last_aligned:
                aligned_fragments += 1
                last_aligned = f[0]
        covered = 0
        for length, op in CIGAR_OP.findall(cigar):
            if op in COVERED:
                covered += int(length)
        if all(covered < min_cov for _, min_cov in levels.values()):
            continue
        nm = None
        for tag in f[11].split("\t") if len(f) > 11 else ():
            if tag.startswith("NM:i:"):
                nm = int(tag[5:])
                break
        if nm is None:
            continue
        kept = [name for name, (max_edit, min_cov) in levels.items()
                if covered >= min_cov and nm * 100 <= max_edit * covered]
        if not kept:
            continue
        parts = f[2].split("|")
        ref = parts[0]
        marker = parts[1] if len(parts) > 1 else ""
        tax = parts[2] if len(parts) > 2 else ""
        row = f"{f[0]}\t{ref}\t{marker}\t{tax}\t{nm}\t{covered}\n"
        for name in kept:
            outputs[name].write(row)
            if collect_stats:
                alignments[name] += 1
                if f[0] != last_kept[name]:
                    kept_fragments[name] += 1
                    last_kept[name] = f[0]
    if collect_stats:
        return {"sample": sample, "records": records, "aligned_fragments": aligned_fragments,
                "levels": {name: {"alignments": alignments[name], "fragments": kept_fragments[name]}
                           for name in levels}}


if __name__ == "__main__":
    main()
