#!/bin/bash
#SBATCH --job-name=liftv2
#SBATCH --output=${CENSV_ROOT}/analyses/meqtl_replication/logs/v2_lift/array/%A_%a.log
#SBATCH --time=4:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -uo pipefail
export OMP_NUM_THREADS=1
P=${CENSV_ROOT}
U=$P/analyses/meqtl_replication
CHUNK=$U/output/v2_lift/chunks/chunk_$(printf '%03d' "$SLURM_ARRAY_TASK_ID").txt

echo "HOST $(hostname)  TASK $SLURM_ARRAY_TASK_ID  START $(date -u +%Y-%m-%dT%H:%M:%SZ)"
[ -f "$CHUNK" ] || { echo "FAIL chunk file missing: $CHUNK"; exit 1; }
N=$(wc -l < "$CHUNK")
echo "chunk holds $N haplotypes"

cd "$P" || exit 1
~/.pixi/bin/pixi run python "$U/code/07_lift_targets/lift_targets.py" $(cat "$CHUNK" | tr '\n' ' ')
rc=$?
echo "EXIT $rc  END $(date -u +%Y-%m-%dT%H:%M:%SZ)"
if [ $rc -eq 0 ]; then echo PASS; else echo "FAIL rc=$rc"; fi
exit $rc
