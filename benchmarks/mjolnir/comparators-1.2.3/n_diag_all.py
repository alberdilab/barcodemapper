"""Read-only N-run diagnostic for several methods in ONE pass over a reference.
Usage: zcat REF.fa.gz | python -I n_diag_all.py OUT.tsv LABEL:FORMAT:FILE [...]
FORMAT: vsearch  (accepted vsearch userout, caln column 12)
        paf      (accepted PAF: target col6, start col8, end col9)
        sam      (`samtools view | cut -f1-6,12-` text; the BarcodeMapper rule covered>=80 and NM*100<=3*covered is applied here)
For every accepted alignment the aligned reference stretch is cut from the reference and its N content counted.
Output columns: label, accepted_alignments, with_any_N, with_ge50pct_N, frac_ge50pct_N."""
import re, sys
from collections import defaultdict

CIG = re.compile(r"(\d+)([MIDNSHP=X])")
VOPS = re.compile(r"(\d*)([MID])")
out_path, specs = sys.argv[1], sys.argv[2:]
regions = defaultdict(list)   # label -> [(subject, start, span)]
for spec in specs:
    label, fmt, path = spec.split(":", 2)
    for line in open(path):
        r = line.rstrip("\n").split("\t")
        if not r or not r[0]:
            continue
        if fmt == "vsearch":
            ops = [(int(n) if n else 1, o) for n, o in VOPS.findall(r[11])]
            start = ops[0][0] if ops[0][1] == "I" else 0
            inner = ops[1:-1] if ops[-1][1] == "I" and len(ops) > 1 else ops[1:] if ops[0][1] == "I" else ops
            regions[label].append((r[1], start, sum(n for n, o in inner if o in "MI")))
        elif fmt == "paf":
            regions[label].append((r[5], int(r[7]), int(r[8]) - int(r[7])))
        elif fmt == "sam":
            if int(r[1]) & 4 or r[5] == "*":
                continue
            ops = [(int(n), o) for n, o in CIG.findall(r[5])]
            covered = sum(n for n, o in ops if o in "M=X")
            nm = next((int(t[5:]) for t in r[6:] if t.startswith("NM:i:")), None)
            if nm is None or covered < 80 or nm * 100 > 3.0 * covered:
                continue
            regions[label].append((r[2], int(r[3]) - 1, sum(n for n, o in ops if o in "MD=X")))
        else:
            raise SystemExit(f"bad format {fmt}")
subjects = {s for v in regions.values() for s, _, _ in v}
seqs, cur = {}, None
for line in sys.stdin:
    if line[0] == ">":
        label = line[1:].strip(); cur = label if label in subjects else None
        if cur: seqs[cur] = []
    elif cur:
        seqs[cur].append(line.strip())
with open(out_path, "x") as o:
    o.write("label\taccepted_alignments\twith_any_N\twith_ge50pct_N\tfrac_ge50pct_N\n")
    for label, regs in regions.items():
        any_n = ge50 = 0
        full = {}
        for s, a, span in regs:
            if s not in full: full[s] = "".join(seqs[s]).upper()
            reg = full[s][a:a + span]
            k = reg.count("N")
            any_n += k > 0
            ge50 += bool(reg) and k / len(reg) >= 0.5
        n = len(regs)
        o.write(f"{label}\t{n}\t{any_n}\t{ge50}\t{ge50 / n if n else float('nan'):.4f}\n")
