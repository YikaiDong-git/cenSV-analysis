#!/bin/bash
#SBATCH -J ssv_fanout -c 4 --mem=16G -t 6:00:00
#SBATCH -o ${CENSV_SCRATCH}/ssv_fanout/logs/%A_%a.log
#SBATCH -e ${CENSV_SCRATCH}/ssv_fanout/logs/%A_%a.log
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_SCRATCH:?set CENSV_SCRATCH (see README)}"
set -e
FAN=${CENSV_SCRATCH}/ssv_fanout
NEW=${CENSV_ROOT}/analyses/finemap_phewide
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $FAN/chunks.txt)
[ -z "$LINE" ] && { echo "no chunk for task $SLURM_ARRAY_TASK_ID"; exit 0; }
REGION=$(echo "$LINE" | cut -f1); CS=$(echo "$LINE" | cut -f2); CE=$(echo "$LINE" | cut -f3); CID=$(echo "$LINE" | cut -f4)
mkdir -p $FAN/out $FAN/logs
source $NEW/code/pixi_env_hook.sh
export RETICULATE_PYTHON=$(which python3)
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMEXPR_NUM_THREADS=4
echo "=== chunk ${SLURM_ARRAY_TASK_ID} $REGION [$CS-$CE] node=$(hostname) START $(date +%H:%M:%S) ==="
Rscript $FAN/fanout_run.R "$REGION" "$FAN/manifest_fanout.tsv" "$FAN/out" "$FAN" "$CS" "$CE" "$CID"
echo "=== $REGION c$CID DONE $(date +%H:%M:%S) ==="
