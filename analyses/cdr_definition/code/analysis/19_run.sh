#!/bin/bash
#SBATCH --job-name=gao_jack
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --output=${CENSV_ROOT}/analyses/cdr_definition/logs/gao_jack_%j.log
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail
export OMP_NUM_THREADS=1
cd ${CENSV_ROOT}
~/.pixi/bin/pixi run python analyses/cdr_definition/code/analysis/19_dip_jackknife.py
echo JACK_DONE
