#!/bin/bash
#SBATCH --job-name=winavail
#SBATCH --output=${CENSV_ROOT}/analyses/meqtl_replication/logs/window_availability/winavail_%j.log
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -uo pipefail
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2
P=${CENSV_ROOT}
U=$P/analyses/meqtl_replication
echo "HOST $(hostname)  START $(date -u +%Y-%m-%dT%H:%M:%SZ)"
cd "$P" || exit 1
~/.pixi/bin/pixi run python "$U/code/11_availability/window_availability_control.py"
rc=$?
echo "EXIT $rc  END $(date -u +%Y-%m-%dT%H:%M:%SZ)"
[ $rc -eq 0 ] && echo PASS || echo "FAIL rc=$rc"
exit $rc
