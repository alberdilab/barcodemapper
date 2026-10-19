"""Read-only diagnostic: share of accepted vsearch hits whose aligned target region contains N.
Usage: zcat REF.fa.gz | python -I vsearch_n_diagnostic.py ACCEPTED.tsv
vsearch caln: leading 'I' consumes target only (terminal), M consumes both."""
import re, sys
OPS = re.compile(r"(\d*)([MID])")
hits = []
subjects = set()
for line in open(sys.argv[1]):
    r = line.rstrip("\n").split("\t")
    ops = [(int(n) if n else 1, o) for n, o in OPS.findall(r[11])]
    start = ops[0][0] if ops[0][1] == "I" else 0
    inner = ops[1:-1] if ops[-1][1] == "I" and len(ops) > 1 else ops[1:] if ops[0][1] == "I" else ops
    span = sum(n for n, o in inner if o in "MI")
    hits.append((r[1], start, span)); subjects.add(r[1])
seqs = {}
cur = None
for line in sys.stdin:
    if line[0] == ">":
        label = line[1:].strip(); cur = label if label in subjects else None
        if cur: seqs[cur] = []
    elif cur:
        seqs[cur].append(line.strip())
n_hits = len(hits); with_n = 0; mostly_n = 0
for s, start, span in hits:
    reg = "".join(seqs[s])[start:start + span].upper()
    k = reg.count("N")
    if k: with_n += 1
    if reg and k / len(reg) >= 0.5: mostly_n += 1
print(f"accepted_hits={n_hits} hits_with_any_N_in_aligned_target={with_n} hits_with_>=50%_N={mostly_n}")
