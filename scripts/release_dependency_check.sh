#!/usr/bin/env bash
# Release-time dependency check: are there newer versions of the pinned
# dependencies, and does the pipeline still give identical results with them?
#
# Looks up the latest versions, and if any are newer, runs the regression
# check in a throwaway environment on Mjolnir. Run before tagging a release.
#
# BM_HPC_SSH is split on whitespace, so it cannot hold quoted arguments.
#
# usage: scripts/release_dependency_check.sh [--force] [--keep]
#   --force  run the HPC check even when every pin is already the latest
#   --keep   keep the throwaway environment on the HPC for inspection
#
# Environment:
#   BM_HPC_SSH      ssh command prefix (default: ssh mjolnirgate.unicph.domain)
#   BM_HPC_SCRATCH  remote scratch directory (default: /projects/alberdilab/scratch/jpl786)
set -uo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
force=
keep=
for arg in "$@"; do
    case "$arg" in
        --force) force=1 ;;
        --keep) keep=--keep ;;
        *) echo "unknown option: $arg" >&2; exit 64 ;;
    esac
done

# shellcheck disable=SC2206
ssh_cmd=(${BM_HPC_SSH:-ssh mjolnirgate.unicph.domain})
scratch=${BM_HPC_SCRATCH:-/projects/alberdilab/scratch/jpl786}
candidate=$(mktemp -t bm_candidate.XXXXXX)
mv "$candidate" "$candidate.yml"
candidate=$candidate.yml
remote=
trap 'rm -f "$candidate"; [ -n "$remote" ] && "${ssh_cmd[@]}" "rm -rf \"$remote\"" >/dev/null 2>&1' EXIT

echo "== dependency versions (pinned vs latest)"
python3 "$root/scripts/check_dependencies.py" --write-candidate "$candidate" --fail-if-outdated
code=$?
case "$code" in
    0) echo "All pins are already the latest versions." ;;
    10) ;;
    *) echo "Version lookup failed (exit $code); cannot decide." >&2; exit "$code" ;;
esac
if [ "$code" -eq 0 ] && [ -z "$force" ]; then
    echo "Nothing newer to test (use --force to run the regression check anyway)."
    exit 0
fi

stamp=$(date +%Y%m%d%H%M%S)
remote=$scratch/tmp_bm_release_check_$stamp
echo "== copying the checkout to the HPC ($remote)"
(cd "$root" && COPYFILE_DISABLE=1 tar --exclude='._*' --exclude=.git --exclude=build --exclude=dist --exclude=.codeome \
    --exclude=.claude --exclude='*.egg-info' --exclude=benchmarks --exclude=__pycache__ \
    --exclude=.DS_Store -cf - .) |
    "${ssh_cmd[@]}" "mkdir -p '$remote/repo' && tar -xf - -C '$remote/repo'" || exit 4
"${ssh_cmd[@]}" "cat > '$remote/candidate.yml'" < "$candidate" || exit 4

echo "== regression check on the HPC with the newer versions"
"${ssh_cmd[@]}" "BM_SCRATCH='$scratch' bash '$remote/repo/scripts/mjolnir_dependency_check.sh' '$remote/repo' '$remote/candidate.yml' $keep"
result=$?

echo
if [ "$result" -eq 0 ]; then
    echo "VERDICT: the pipeline gives identical results with the newer dependencies."
    echo "         Consider adopting the versions above in barcodemapper_env.yml."
else
    echo "VERDICT: do NOT adopt the newer dependencies yet (see the output above)."
fi
exit "$result"
