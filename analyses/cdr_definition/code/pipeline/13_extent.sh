#!/bin/bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"
#SBATCH --job-name=gao_ext
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --output=${CENSV_ROOT}/analyses/cdr_definition/logs/ext/gao_ext_%A_%a.log
set -euo pipefail
export OMP_NUM_THREADS=1

B=${CENSV_WORKSPACE}/tools/CDR-Finder
PY=$B/env/.pixi/envs/default/bin/python
CF=$B/upstream/workflow/scripts/cdr_finder.py
SRC=${CENSV_ROOT}/analyses/cdr_calls/output/calls
U=${CENSV_ROOT}/analyses/cdr_definition
OUT=$U/output/callsets/ext_dips

CHUNK=24
START=$(( (SLURM_ARRAY_TASK_ID - 1) * CHUNK + 1 ))
TAGS=$(ls "$SRC" | sed -n "${START},$((START + CHUNK - 1))p")
[ -z "$TAGS" ] && { echo "no tags for task ${SLURM_ARRAY_TASK_ID}"; exit 0; }

for TAG in $TAGS; do
  IN=$SRC/$TAG/binned_aa.bed
  [ -s "$IN" ] || { echo "[skip] $TAG no binned_aa.bed"; continue; }
  mkdir -p "$OUT/$TAG"
  $PY "$CF" -i "$IN" -o "$OUT/$TAG/cdrs.bed" \
    --thr_height_perc_valley 0.34 --thr_prom_perc_valley 0.30 \
    --bp_merge 1 --bp_edge 500000 --baseline_avg_methyl 0.4 \
    --edge_height_heuristic min --extend_edges_std -1 2>/dev/null
  echo "[ok] $TAG  canonical=$(wc -l < "$SRC/$TAG/cdrs.bed")  extended=$(wc -l < "$OUT/$TAG/cdrs.bed")"
done
echo "EXT_CHUNK_DONE ${SLURM_ARRAY_TASK_ID}"
