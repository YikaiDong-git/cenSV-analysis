#!/bin/bash
#SBATCH --job-name=fg_ear_harm
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH -c 2
#SBATCH --output=${CENSV_ROOT}/analyses/disease_replication/finngen_replication/logs/harm_%A_%a.log
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"
set -euo pipefail
BASE=${CENSV_WORKSPACE}
PROJ=$BASE/Project1_Centromere/centromere_twas/finngen
PIXI_PY=$BASE/_shared/tools/bayesian_finemap/v1_20260515/.pixi/envs/default/bin/python
MANIFEST=${MANIFEST:?set MANIFEST via --export}
ROW=$((SLURM_ARRAY_TASK_ID + 1))
PHENOCODE=$(awk -v r="$ROW" -F'\t' 'NR==r{print $1}' "$MANIFEST")
N_CASES=$(awk -v r="$ROW" -F'\t' 'NR==r{print $2}' "$MANIFEST")
N_CTRLS=$(awk -v r="$ROW" -F'\t' 'NR==r{print $3}' "$MANIFEST")
[ -z "$PHENOCODE" ] && { echo "no phenocode at row $ROW"; exit 2; }
echo "[task $SLURM_ARRAY_TASK_ID] harmonize $PHENOCODE ($N_CASES ca / $N_CTRLS co)"
$PIXI_PY $PROJ/scripts/harmonize_finngen_one.py "$PHENOCODE" "$N_CASES" "$N_CTRLS"
