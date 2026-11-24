#!/bin/bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail

PROJ=${CENSV_ROOT}
UNIT=$PROJ/analyses/meqtl_replication
LIST=$UNIT/output/a1_lift/${1:?tag list filename required}

TAG=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
if [ -z "$TAG" ]; then
  echo "no tag at line ${SLURM_ARRAY_TASK_ID} of $LIST" >&2
  exit 1
fi

echo "=== $(date -Is)  task ${SLURM_ARRAY_TASK_ID}  tag ${TAG}  node $(hostname) ==="
cd "$PROJ"
~/.pixi/bin/pixi run python "$UNIT/code/01_lift/lift_targets_percpg.py" "$TAG"
echo "=== $(date -Is)  done ${TAG} ==="
