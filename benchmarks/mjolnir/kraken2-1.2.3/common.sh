# Sourced by the Kraken2 scripts. Needs BENCH_ROOT, RUN, SETTING (kraken2|kraken2-conf0.1).
set -euo pipefail
P="$BENCH_ROOT/projects/03-method-comparison"
TASKS="$BENCH_ROOT/tasks/1.2.3-kraken2"
WORK="$P/work/comparators-1.2.3"
RES="$P/results/pilots-1.2.3"
READS="$P/work/pilots-1.2.3/$RUN-pilot-standard/fastp"
case "$RUN" in
  SRR14351851) DBNAME=all-fungi ;;
  ERR6001161|ERR6001173) DBNAME=default ;;
  *) echo "unknown run $RUN" >&2; exit 2 ;;
esac
case "$SETTING" in
  kraken2) CONF=0.0 ;;
  kraken2-conf0.1) CONF=0.1 ;;
  *) echo "unknown setting $SETTING" >&2; exit 2 ;;
esac
DB="$P/work/kraken2/$DBNAME"
CONFIG="$WORK/comparator_config.kraken2.yaml"
OUT="$RES/$SETTING"
