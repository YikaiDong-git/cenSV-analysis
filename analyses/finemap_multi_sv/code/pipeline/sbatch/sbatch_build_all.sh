#!/bin/bash
#SBATCH --job-name=msv_build_v2
#SBATCH --array=1-21%6
#SBATCH --mem=32G
#SBATCH -c 4
#SBATCH --time=2:00:00
#SBATCH --output=${CENSV_ROOT}/analyses/finemap_multi_sv/logs/full_restart/build/build_v2.%A_%a.out
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"

set -euo pipefail

BASE=${CENSV_ROOT}/analyses/finemap_multi_sv
PIXI=~/.pixi/bin/pixi
PIXI_MANIFEST=${CENSV_WORKSPACE}/_shared/tools/bayesian_finemap/v1_20260515/pixi.toml
REGION_MAP=$BASE/output/full_restart/manifests/region_map_v2_with_admit.tsv

cd "$BASE"
mkdir -p logs/full_restart/build output/full_restart/bundles

ROW=$(awk -v idx=$SLURM_ARRAY_TASK_ID -F'\t' 'NR>1 && $7>=1 {n++; if (n==idx) print $0}' "$REGION_MAP")
if [ -z "$ROW" ]; then
  echo "no row for task=$SLURM_ARRAY_TASK_ID (after filtering n_admit_v2>=1)"
  exit 0
fi

CHROM=$(echo "$ROW" | awk -F'\t' '{print $1}')
REGION_ID=$(echo "$ROW" | awk -F'\t' '{print $2}')
RSTART=$(echo "$ROW" | awk -F'\t' '{print $3}')
REND=$(echo "$ROW" | awk -F'\t' '{print $4}')
N_ADMIT=$(echo "$ROW" | awk -F'\t' '{print $7}')

CHROM_TAG="chr${CHROM}"
OUT_DIR="$BASE/output/full_restart/bundles/${CHROM_TAG}"
mkdir -p "$OUT_DIR"

BUNDLE_PATH="$OUT_DIR/${REGION_ID}_bundle.npz"
if [ -f "$BUNDLE_PATH" ]; then
  echo "[$(date)] bundle v2 exists, skipping: $BUNDLE_PATH"
  exit 0
fi

echo "[$(date)] task=$SLURM_ARRAY_TASK_ID V2 building $REGION_ID (chr$CHROM, ${RSTART}-${REND}, n_admit_v2=$N_ADMIT)"

$PIXI run --manifest-path "$PIXI_MANIFEST" python \
    "$BASE/code/pipeline/02_build_bundles/build_multisv_bundle.py" \
    --chrom "$CHROM_TAG" \
    --region-id "$REGION_ID" \
    --region-start $RSTART \
    --region-end $REND \
    --out-dir "$OUT_DIR"

echo "[$(date)] DONE V2 bundle $REGION_ID"
ls -la "$BUNDLE_PATH"
