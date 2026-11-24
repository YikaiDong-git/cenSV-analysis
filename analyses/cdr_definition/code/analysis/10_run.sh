#!/bin/bash
#SBATCH --job-name=gao_percall
#SBATCH --mem=24G
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --output=${CENSV_ROOT}/analyses/cdr_definition/logs/gao_percall_%j.log
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail
export OMP_NUM_THREADS=1
cd ${CENSV_ROOT}
~/.pixi/bin/pixi run python analyses/cdr_definition/code/analysis/10_per_unit_positions.py
echo PERCALL_DONE
