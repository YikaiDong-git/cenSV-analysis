#!/bin/bash
#SBATCH -J svfu_light
#SBATCH -c 1
#SBATCH --mem 16G
#SBATCH -t 2:00:00
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail
D=${CENSV_ROOT}/analyses/sv_call_reliability
PY=${CENSV_ROOT}/.pixi/envs/default/bin/python
cd $D
echo "== A3 targets"; $PY code/a3_seqval/resolve_targets.py output/a3_seqval
echo "== A5";         $PY code/a5_robustness/robustness.py output/a5_robustness
echo "== A6";         $PY code/a6_count/reconcile.py output/a6_count
echo "== A4";         $PY code/a4_determinism/determinism.py output/a4_determinism
echo LIGHT_DONE
