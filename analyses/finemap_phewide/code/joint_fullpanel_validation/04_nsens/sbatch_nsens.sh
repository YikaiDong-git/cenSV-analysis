#!/bin/bash
#SBATCH --job-name=nsens -c 4 --mem=24G -t 6:00:00
#SBATCH --output=${CENSV_SCRATCH}/ssv_joint/logs/nsens.%A_%a.out
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_SCRATCH:?set CENSV_SCRATCH (see README)}"
set -e
FAN=${CENSV_SCRATCH}/ssv_joint
NEW=${CENSV_ROOT}/analyses/finemap_phewide
CHUNKS=$FAN/chunks_nsens.txt
MANIFEST=$FAN/manifest_nsens.tsv
OUTDIR=$FAN/nsens_out
mkdir -p $OUTDIR $FAN/logs
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $CHUNKS)
[ -z "$LINE" ] && { echo "no chunk"; exit 0; }
REGION=$(echo "$LINE" | cut -f1); CS=$(echo "$LINE" | cut -f2); CE=$(echo "$LINE" | cut -f3); CID=$(echo "$LINE" | cut -f4)
OUTF=$OUTDIR/${REGION}_c${CID}_nsens_diag.tsv
[ -s "$OUTF" ] && { echo "exists"; exit 0; }
source $NEW/code/pixi_env_hook.sh
export RETICULATE_PYTHON=$(which python3)
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMEXPR_NUM_THREADS=4
echo "=== nsens chunk $SLURM_ARRAY_TASK_ID $REGION [$CS-$CE] node=$(hostname) ==="
Rscript $FAN/fanout_run_n.R "$REGION" "$MANIFEST" "$OUTDIR" "$FAN" "$CS" "$CE" "$CID"
echo "=== nsens $REGION c$CID DONE ==="
