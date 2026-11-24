#!/bin/bash
#SBATCH --job-name=aou_extract
#SBATCH --mem=12G
#SBATCH -c 1
#SBATCH --nodes=1
#SBATCH --time=0:20:00
#SBATCH --output=${CENSV_ROOT}/analyses/disease_replication/logs/aou_extract_%j.log
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"
set -uo pipefail
BT=${CENSV_WORKSPACE}/_shared/envs/workspace/v1/.pixi/envs/biotools/bin
DD=${CENSV_ROOT}/analyses/disease_replication
$BT/python $DD/code/extract_aou_weights.py
echo "[sbatch] exit=$?"
