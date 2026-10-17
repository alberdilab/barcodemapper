#!/usr/bin/env bash
# Use the login-node Slurm configuration for accounting queries on compute nodes.
set -euo pipefail
script_path=$(readlink -f -- "$0")
export SLURM_CONF="$(dirname -- "$script_path")/slurm-accounting.conf"
exec /usr/bin/sacct "$@"
