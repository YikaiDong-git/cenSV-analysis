#!/bin/bash
#SBATCH --job-name=msv_fm_v2
#SBATCH --mem=48G
#SBATCH -c 4
#SBATCH --time=4:00:00
#SBATCH --output=${CENSV_ROOT}/analyses/finemap_multi_sv/logs/full_restart/finemap/fm_v2.%A_%a.out
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"

set -euo pipefail

BASE=${CENSV_ROOT}/analyses/finemap_multi_sv
MANIFEST=$BASE/output/full_restart/manifests/finemap_manifest_v2.tsv

cd "$BASE"
mkdir -p logs/full_restart/finemap output/full_restart/finemap

ROW=$(awk -v idx=$SLURM_ARRAY_TASK_ID -F'\t' 'NR==idx+1' "$MANIFEST")
CHROM=$(echo "$ROW" | awk -F'\t' '{print $2}')
REGION_ID=$(echo "$ROW" | awk -F'\t' '{print $3}')
PHENO=$(echo "$ROW" | awk -F'\t' '{print $4}')
STUDY=$(echo "$ROW" | awk -F'\t' '{print $5}')

CHROM_TAG="chr${CHROM}"
BUNDLE="$BASE/output/full_restart/bundles/${CHROM_TAG}/${REGION_ID}_bundle.npz"
OUT_DIR="$BASE/output/full_restart/finemap/${CHROM_TAG}"

if [ ! -f "$BUNDLE" ]; then
  echo "ERROR: V2 bundle missing for $REGION_ID: $BUNDLE"
  exit 1
fi

bash "$BASE/code/pipeline/03_finemap/prepare_and_finemap_one.sh" \
    "$BUNDLE" "$PHENO" "$STUDY" "$OUT_DIR" "$REGION_ID"
