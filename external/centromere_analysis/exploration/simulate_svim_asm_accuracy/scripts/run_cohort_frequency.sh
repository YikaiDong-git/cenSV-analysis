#!/usr/bin/env bash
: "${CENSV_COLLAB_ROOT:?set CENSV_COLLAB_ROOT (see README)}"
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail
cd ${CENSV_COLLAB_ROOT}/exploration/simulate_svim_asm_accuracy
case "$(hostname -s)" in node03|node04|node08|node12) ;; *) echo 'Node outside authorized allowlist' >&2; exit 2;; esac
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
task_scratch=/scratch/jwein22/sv_cohort_frequency
mkdir -p "$task_scratch"
TMPDIR=$(mktemp -d "$task_scratch/job_${SLURM_JOB_ID}_XXXXXX")
export TMPDIR
trap 'rm -rf -- "$TMPDIR"' EXIT
export PYTHONPYCACHEPREFIX="$TMPDIR/pycache" MPLCONFIGDIR="$TMPDIR/mpl"
task_python=${CENSV_ROOT}/.pixi/envs/default/bin/python
if [[ "$1" == pilot ]]; then
    "$task_python" scripts/cohort_frequency_simulation.py pilot
elif [[ "$1" == run ]]; then
    "$task_python" scripts/cohort_frequency_simulation.py run --task "${SLURM_ARRAY_TASK_ID:?}"
elif [[ "$1" == analyze ]]; then
    task_children=()
    for task_chrom in chr1 chr7 chr13 chr19; do
        "$task_python" scripts/cohort_frequency_metrics.py --manifest results/cohort_frequency_v1/manifest.json --chrom "$task_chrom" &
        task_children+=("$!")
    done
    for task_child in "${task_children[@]}"; do wait "$task_child"; done
    "$task_python" scripts/cohort_frequency_metrics.py --manifest results/cohort_frequency_v1/manifest.json --aggregate
else
    echo 'Unknown cohort action' >&2
    exit 2
fi
