#!/usr/bin/env bash
: "${CENSV_COLLAB_ROOT:?set CENSV_COLLAB_ROOT (see README)}"
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail
cd ${CENSV_COLLAB_ROOT}/exploration/simulate_svim_asm_accuracy
case "$(hostname -s)" in node03|node04|node08|node12) ;; *) exit 2;; esac
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
task_scratch=/scratch/jwein22/sv_cohort_joint
mkdir -p "$task_scratch"
TMPDIR=$(mktemp -d "$task_scratch/job_${SLURM_JOB_ID}_XXXXXX")
export TMPDIR
trap 'rm -rf -- "$TMPDIR"' EXIT
export PYTHONPYCACHEPREFIX="$TMPDIR/pycache" MPLCONFIGDIR="$TMPDIR/mpl"
task_python=${CENSV_ROOT}/.pixi/envs/default/bin/python
task_root=results/cohort_frequency_v1

joint_pipeline() {
    "$task_python" - <<'PY' || return
from cohort_frequency_metrics import ensure_complete_cohort
ensure_complete_cohort('results/cohort_frequency_v1/manifest.json')
PY
    "$task_python" scripts/canonicalize_cohort_raw_scores.py || return
    "$task_python" scripts/prepare_joint_cohort.py --root "$task_root" --workers 4 || return
    if [[ ! -f "$task_root/joint_analysis/analysis_complete.json" ]]; then
        "$task_python" - <<'PY' || return
import json
from pathlib import Path
from cohort_frequency_joint_metrics import analyze_joint
from cohort_frequency_joint_reference import rescore_joint_thresholds
root = Path('results/cohort_frequency_v1')
result = analyze_joint(root/'manifest.json', rescore_joint_thresholds)
(root/'joint_analysis/analysis_complete.json').write_text(json.dumps(dict(completed=True,
    rows=result['n_rows'], source_carrier_assertions=result['n_memberships']), indent=2)+'\n')
print('Joint frequency analysis complete', result['n_rows'], result['n_memberships'], flush=True)
PY
    fi
    "$task_python" scripts/validate_cohort_joint_semantics.py --manifest "$task_root/manifest.json" \
        --output "$task_root/joint_analysis/merged_semantics_validation.json" || return
    "$task_python" scripts/audit_cohort_joint_run.py --root "$task_root" || return
    "$task_python" scripts/stratify_cohort_frequency.py --manifest "$task_root/manifest.json" || return
    "$task_python" scripts/build_cohort_frequency_report.py --manifest "$task_root/manifest.json" \
        --output-dir "$task_root/report_frequency_joint_v1" || return
    if [[ ! -f results/multichrom_diagnostics_v2/report_frequency_v4/complete.json ]]; then
        "$task_python" scripts/build_multichrom_sv_report.py --run-dir results/multichrom_sv_v1 \
            --summary-dir results/multichrom_diagnostics_v2/summary \
            --output results/multichrom_diagnostics_v2/report_frequency_v4 \
            --cohort-report "$task_root/report_frequency_joint_v1" --finalize || return
    fi
    "$task_python" - <<'PY' || return
from pathlib import Path
import hashlib,json
r=Path('results/cohort_frequency_v1')
pdf=r/'report_frequency_joint_v1/latex_report/report.pdf'
(r/'joint_analysis/complete.json').write_text(json.dumps(dict(completed=True,
    pdf_sha256=hashlib.sha256(pdf.read_bytes()).hexdigest()),indent=2)+'\n')
PY
}

export PYTHONPATH="$PWD/scripts${PYTHONPATH:+:$PYTHONPATH}"
while true; do
    if joint_pipeline; then
        rm -f "$task_root/joint_stage_failed.flag"
        exit 0
    fi
    date -u > "$task_root/joint_stage_failed.flag"
    echo 'Joint stage failed; retaining this final allocation for a reviewed repair and retry flag.' >&2
    while [[ ! -f "$task_root/joint_stage_retry.flag" ]]; do sleep 30; done
    rm -f "$task_root/joint_stage_retry.flag"
done
