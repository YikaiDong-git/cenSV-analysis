#!/usr/bin/env bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail
D=${CENSV_ROOT}/analyses/cdr_calls
LIST=${1:-$D/code/polish_haps.txt}
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
set -- $line
bash $D/code/run_cdr_one.sh "$1" "$2"
