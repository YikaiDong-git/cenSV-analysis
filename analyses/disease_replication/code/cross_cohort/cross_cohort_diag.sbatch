#!/bin/bash
#SBATCH --job-name=xcohort_diag
#SBATCH --mem=16G
#SBATCH -c 1
#SBATCH --nodes=1
#SBATCH --time=0:30:00
#SBATCH --output=${CENSV_ROOT}/analyses/disease_replication/logs/xcohort_diag_%j.log
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"
set -uo pipefail
P=${CENSV_WORKSPACE}
BT=$P/_shared/envs/workspace/v1/.pixi/envs/biotools/bin
D=$P/Project1_Centromere/analyses/disease_replication
$BT/python $D/code/cross_cohort_diag.py
echo "[sbatch] exit=$?"
