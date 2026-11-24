#!/bin/bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
#SBATCH --job-name=gao_xcheck
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --output=${CENSV_ROOT}/analyses/cdr_definition/logs/gao_xcheck_%j.log
set -euo pipefail
export OMP_NUM_THREADS=1
U=${CENSV_ROOT}/analyses/cdr_definition
SRC=${CENSV_ROOT}/analyses/cdr_calls/output/calls
TAGS=$(ls "$SRC" | awk 'NR % 16 == 1' | tr '\n' ' ')
echo "[$(date)] tags: $TAGS"
python3 "$U/code/pipeline/02_merge_crosscheck.py" $TAGS
echo "[$(date)] XCHECK_DONE"
