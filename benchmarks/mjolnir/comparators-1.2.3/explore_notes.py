"""Exploratory, read-only: (a) 10 fungal pairs hit by vsearch (n_mismatch) but not by BarcodeMapper standard;
(b) Kraken2 k-mer hit statistics from output.txt. Writes plain text to stdout."""
import csv, re, sys
from collections import defaultdict
bench = sys.argv[1]
P = f"{bench}/projects/03-method-comparison"; R = f"{P}/results/pilots-1.2.3"; N = f"{P}/work/comparators-1.2.3/ndiag"
OPS = re.compile(r"(\d*)([MIDNSHP=X])")
run = "SRR14351851"

def calls(path):
    return {r["ReadID"]: r for r in csv.DictReader(open(path), delimiter="\t")}
v = calls(f"{R}/evaluation/calls/{run}.vsearch_n_mismatch.tsv"); b = calls(f"{R}/evaluation/calls/{run}.barcodemapper_standard.tsv")
only = sorted(k for k, r in v.items() if r["status"] != "no_hit" and b[k]["status"] == "no_hit")
print(f"(a) {run}: vsearch(n_mismatch) hit pairs {sum(r['status']!='no_hit' for r in v.values())}, BarcodeMapper standard {sum(r['status']!='no_hit' for r in b.values())}, vsearch-only {len(only)}")
step = len(only) / 10
pick = [only[int(i * step)] for i in range(10)]
hits = defaultdict(list)
for line in open(f"{R}/vsearch-nmismatch/{run}.accepted.tsv"):
    r = line.rstrip("\n").split("\t"); hits[r[0].rsplit("/", 1)[0]].append(r)
bam = defaultdict(list)
for line in open(f"{N}/{run}.bm.sam.txt"):
    r = line.rstrip("\n").split("\t"); bam[r[0]].append(r)
cat = defaultdict(int)
for rid in pick:
    print(f"\n{rid}: vsearch call {v[rid]['Taxonomy'].split(';')[-1]}")
    for r in hits[rid][:3]:
        ops = [(int(n) if n else 1, o) for n, o in OPS.findall(r[11])]
        cov = sum(n for n, o in ops if o == "M")
        print(f"  vsearch {r[0][-2:]} target={r[1].split('|')[0]}|{r[1].split('|')[1]} id={r[2]} covered={cov} mism={r[4]} gaps={r[6]} NM_approx={int(r[4])+0} ql={r[7]} qlo-qhi={r[8]}-{r[9]} tl={r[10]} qcov={r[12]}")
    al = bam.get(rid, [])
    if not al:
        print("  Bowtie2: no alignment for either mate in the BAM (read not aligned end-to-end)"); cat["no Bowtie2 alignment at all"] += 1
    for r in al[:4]:
        cig = [(int(n), o) for n, o in re.findall(r"(\d+)([MIDNSHP=X])", r[5])]
        covb = sum(n for n, o in cig if o in "M=X"); nm = next((int(t[5:]) for t in r[6:] if t.startswith("NM:i:")), None)
        why = "accepted" if covb >= 80 and nm * 100 <= 3 * covb else ("covered<80" if covb < 80 else f"NM {nm} > 3% of covered {covb}")
        print(f"  Bowtie2 flag={r[1]} ref={r[2].split('|')[0]} cigar={r[5]} NM={nm} covered={covb} -> {why}")
        if why != "accepted": cat[why.split(" ")[0]] += 1
print("\nbowtie2 mode: end-to-end (no --local), default sensitivity, paired -1/-2; reports unpaired alignments when no concordant pair")
print("category counts over the 10 pairs (first listed alignment):", dict(cat))

print("\n(b) Kraken2 k-mer hits (output.txt; pair has >=1 k-mer with a nonzero taxid in the DB):")
for d in ("kraken2", "kraken2-conf0.1"):
    for rr in ("SRR14351851", "ERR6001161", "ERR6001173"):
        n = anyk = C = totk = hitk = 0
        for line in open(f"{R}/{d}/{rr}.output.txt"):
            f = line.rstrip("\n").split("\t"); n += 1; C += f[0] == "C"
            toks = [t for t in f[4].split() if t != "|:|"]
            tot = hit = 0
            for t in toks:
                tx, c = t.split(":"); c = int(c); tot += c
                if tx not in ("0", "A"): hit += c
            totk += tot; hitk += hit; anyk += hit > 0
        print(f"  {d} {rr}: pairs {n}, classified(C) {C} ({C/n:.4f}), pairs with >=1 DB k-mer {anyk} ({anyk/n:.4f}), k-mers in DB {hitk} of {totk} ({hitk/totk:.5f})")
