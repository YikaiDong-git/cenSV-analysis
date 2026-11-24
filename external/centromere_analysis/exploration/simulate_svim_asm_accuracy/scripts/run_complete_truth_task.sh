#!/usr/bin/env bash
: "${CENSV_COLLAB_ROOT:?set CENSV_COLLAB_ROOT (see README)}"
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail
cd ${CENSV_COLLAB_ROOT}/exploration/simulate_svim_asm_accuracy
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR="${TMPDIR:-/tmp}/complete_truth_mpl_${SLURM_JOB_ID:-local}"
task_python=${CENSV_ROOT}/.pixi/envs/default/bin/python
case "$1" in
  prepare) exec "$task_python" scripts/complete_truth_benchmark.py prepare ;;
  run) exec "$task_python" scripts/complete_truth_benchmark.py run --task "${SLURM_ARRAY_TASK_ID:?}" ;;
  score)
    "$task_python" scripts/evaluate_complete_truth.py "${@:2}"
    exec "$task_python" scripts/build_simulation_report.py --allow-partial
    ;;
  continue) exec "$task_python" scripts/submit_complete_truth.py remaining ;;
  *) exit 2 ;;
esac
