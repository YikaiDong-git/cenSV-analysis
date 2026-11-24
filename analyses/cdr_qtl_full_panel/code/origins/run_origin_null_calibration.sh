#!/bin/bash
#SBATCH --job-name=origin_null
#SBATCH --nodes=1
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 1:00:00
#SBATCH -o ${CENSV_ROOT}/analyses/cdr_qtl_full_panel/output/origin_null_%j.log
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4
P1=${CENSV_ROOT}
cd "$P1"
~/.pixi/bin/pixi run python \
  "$P1/analyses/cdr_qtl_full_panel/code/origins/origin_null_calibration.py"
echo "===ORIGIN_NULL_DONE==="
