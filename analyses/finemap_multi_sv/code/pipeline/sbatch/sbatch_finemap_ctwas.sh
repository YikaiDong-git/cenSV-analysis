#!/bin/bash
#SBATCH --job-name=msv_ctwas
#SBATCH --mem=48G
#SBATCH -c 4
#SBATCH --time=4:00:00
#SBATCH --output=${CENSV_ROOT}/analyses/finemap_multi_sv/logs/full_restart/finemap_ctwas/fm_ctwas.%A_%a.out
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"

set -euo pipefail

BASE=${CENSV_ROOT}/analyses/finemap_multi_sv
PIXI=~/.pixi/bin/pixi
PIXI_MANIFEST=${CENSV_WORKSPACE}/_shared/tools/bayesian_finemap/v1_20260515/pixi.toml
MANIFEST=$BASE/output/full_restart/manifests/finemap_manifest_v2.tsv

cd "$BASE"
mkdir -p logs/full_restart/finemap_ctwas output/full_restart/finemap_ctwas

ROW=$(awk -v idx=$SLURM_ARRAY_TASK_ID -F'\t' 'NR==idx+1' "$MANIFEST")
CHROM=$(echo "$ROW" | awk -F'\t' '{print $2}')
REGION_ID=$(echo "$ROW" | awk -F'\t' '{print $3}')
PHENO=$(echo "$ROW" | awk -F'\t' '{print $4}')
STUDY=$(echo "$ROW" | awk -F'\t' '{print $5}')

CHROM_TAG="chr${CHROM}"
BUNDLE="$BASE/output/full_restart/bundles/${CHROM_TAG}/${REGION_ID}_bundle.npz"
OUT_DIR="$BASE/output/full_restart/finemap_ctwas/${CHROM_TAG}"
mkdir -p "$OUT_DIR"

TMP_BASE=${SLURM_TMPDIR:-/tmp}
WORK="$TMP_BASE/msv_ctwas_${SLURM_JOB_ID:-$$}_${REGION_ID}_${PHENO}"
mkdir -p "$WORK"
trap "rm -rf $WORK" EXIT

DIAG_PATH="$OUT_DIR/${REGION_ID}__${PHENO}__diag.tsv"
if [ -f "$DIAG_PATH" ]; then
  echo "[$(date)] ctwas diag exists, skipping: $DIAG_PATH"
  exit 0
fi

if [ ! -f "$BUNDLE" ]; then
  echo "ERROR: V2 bundle missing for $REGION_ID: $BUNDLE"
  exit 1
fi

echo "[$(date)] cTWAS-faithful prepare jwbundle for $REGION_ID pheno=$PHENO study=$STUDY"
$PIXI run --manifest-path "$PIXI_MANIFEST" python \
    "$BASE/code/pipeline/07_ctwas_faithful/prepare_jwbundle_ctwas.py" \
    --bundle "$BUNDLE" \
    --pheno "$PHENO" \
    --study "$STUDY" \
    --out-dir "$WORK"

JW_PATH="$WORK/${REGION_ID}__${PHENO}__jw_ctwas.npz"
if [ ! -f "$JW_PATH" ]; then
  echo "ERROR: ctwas jwbundle not produced at $JW_PATH"
  exit 1
fi

echo "[$(date)] cTWAS-faithful finemap (L=10, max_iter=100, n=length(z)) $REGION_ID $PHENO"
$PIXI run --manifest-path "$PIXI_MANIFEST" Rscript \
    "$BASE/code/pipeline/07_ctwas_faithful/R/run_susie_rss_ctwas.R" \
    "$JW_PATH" "$OUT_DIR" "$PHENO" "$REGION_ID"

echo "[$(date)] DONE cTWAS-faithful $REGION_ID pheno=$PHENO"
ls -la "$DIAG_PATH"
