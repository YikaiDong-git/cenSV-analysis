#!/bin/bash
#SBATCH --job-name=fg_ear_dl
#SBATCH --time=2:00:00
#SBATCH --mem=4G
#SBATCH -c 2
#SBATCH --output=${CENSV_ROOT}/analyses/disease_replication/finngen_replication/logs/dl_%A_%a.log
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"
set -euo pipefail
BASE=${CENSV_WORKSPACE}
PROJ=$BASE/Project1_Centromere/centromere_twas/finngen
STAGE=$BASE/_shared/datasets/GWAS/finngen/v_R12_20260517
MANIFEST=${MANIFEST:?set MANIFEST via --export}
ROW=$((SLURM_ARRAY_TASK_ID + 1))
PHENOCODE=$(awk -v r="$ROW" -F'\t' 'NR==r{print $1}' "$MANIFEST")
[ -z "$PHENOCODE" ] && { echo "no phenocode at row $ROW"; exit 2; }
echo "[task $SLURM_ARRAY_TASK_ID] downloading $PHENOCODE"
bash $PROJ/scripts/download_one_finngen.sh "$PHENOCODE" "$STAGE"
