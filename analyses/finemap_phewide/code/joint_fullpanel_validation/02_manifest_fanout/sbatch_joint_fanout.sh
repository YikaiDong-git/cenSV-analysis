#!/bin/bash
#SBATCH --job-name=jfan -c 4 --mem=24G -t 12:00:00
#SBATCH --output=${CENSV_SCRATCH}/ssv_joint/logs/jfan.%A_%a.out
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_SCRATCH:?set CENSV_SCRATCH (see README)}"
set -e
FAN=${CENSV_SCRATCH}/ssv_joint
NEW=${CENSV_ROOT}/analyses/finemap_phewide
: "${VER:?set VER=0p3 or 0p1}"
CHUNKS=$FAN/chunks_${VER}.txt
MANIFEST=$FAN/manifest_joint_${VER}.tsv
OUTDIR=$FAN/fanout_out/${VER}
mkdir -p $OUTDIR $FAN/logs
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $CHUNKS)
[ -z "$LINE" ] && { echo "no chunk for task $SLURM_ARRAY_TASK_ID"; exit 0; }
REGION=$(echo "$LINE" | cut -f1); CS=$(echo "$LINE" | cut -f2); CE=$(echo "$LINE" | cut -f3); CID=$(echo "$LINE" | cut -f4)
OUTF=$OUTDIR/${REGION}_c${CID}_fanout_diag.tsv
if [ -s "$OUTF" ]; then echo "exists, skip: $OUTF"; exit 0; fi
source $NEW/code/pixi_env_hook.sh
export RETICULATE_PYTHON=$(which python3)
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMEXPR_NUM_THREADS=4
echo "=== chunk ${SLURM_ARRAY_TASK_ID} $REGION [$CS-$CE] c$CID node=$(hostname) START $(date +%H:%M:%S) ==="
Rscript $FAN/fanout_run.R "$REGION" "$MANIFEST" "$OUTDIR" "$FAN" "$CS" "$CE" "$CID"
echo "=== $REGION c$CID DONE $(date +%H:%M:%S) ==="
