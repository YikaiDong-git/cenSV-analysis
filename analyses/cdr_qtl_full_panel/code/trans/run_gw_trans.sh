#!/usr/bin/env bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail
P1=${CENSV_ROOT}
WU=$P1/analyses/cdr_qtl_full_panel
NPERM=${NPERM:-100000}
SEED=${SEED:-20260805}

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$WU/output/trans/gw_tasks.txt")
[ -z "$LINE" ] && { echo "no task for index ${SLURM_ARRAY_TASK_ID}"; exit 1; }
read -r PHENO CHROM <<< "$LINE"
echo "[$(date)] task ${SLURM_ARRAY_TASK_ID}: $PHENO / $CHROM  nperm=$NPERM"
cd "$P1"
~/.pixi/bin/pixi run python -u "$WU/code/trans/genomewide_trans_scan.py" "$PHENO" "$CHROM" "$NPERM" "$SEED"
echo "[$(date)] DONE $PHENO / $CHROM"
