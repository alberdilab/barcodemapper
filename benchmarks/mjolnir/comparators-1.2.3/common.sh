# Sourced by the comparator scripts. Needs BENCH_ROOT and RUN.
set -euo pipefail
P="$BENCH_ROOT/projects/03-method-comparison"
TASKS="$BENCH_ROOT/tasks/1.2.3-comparators"
WORK="$P/work/comparators-1.2.3"
RES="$P/results/pilots-1.2.3"
READS="$P/work/pilots-1.2.3/$RUN-pilot-standard/fastp"
case "$RUN" in
  SRR14351851) REFNAME=all-fungi ;;
  ERR6001161|ERR6001173) REFNAME=default ;;
  *) echo "unknown run $RUN" >&2; exit 2 ;;
esac
REF="$BENCH_ROOT/references-1.2.3/$REFNAME.fa.gz"
CONFIG="$WORK/comparator_config.yaml"
QDIR="$WORK/queries"
