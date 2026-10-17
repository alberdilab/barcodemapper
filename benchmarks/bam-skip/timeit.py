#!/usr/bin/env python3
"""Run a shell command and append `variant rep elapsed user sys max_rss_kb` to a TSV (no /usr/bin/time needed)."""
import resource
import subprocess
import sys
import time

out, variant, rep, command = sys.argv[1:5]
start = time.time()
code = subprocess.call(["bash", "-c", command])
elapsed = time.time() - start
usage = resource.getrusage(resource.RUSAGE_CHILDREN)
with open(out, "a") as handle:
    handle.write(f"{variant}\t{rep}\t{elapsed:.2f}\t{usage.ru_utime:.2f}\t{usage.ru_stime:.2f}\t{usage.ru_maxrss}\n")
sys.exit(code)
