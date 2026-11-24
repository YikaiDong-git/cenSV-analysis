#!/bin/bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"

set -euo pipefail

BUNDLE=$1
PHENO=$2
STUDY=$3
OUT_DIR=$4
REGION_ID=$5

BASE=${CENSV_ROOT}/analyses/finemap_multi_sv
PIXI=~/.pixi/bin/pixi
PIXI_MANIFEST=${CENSV_WORKSPACE}/_shared/tools/bayesian_finemap/v1_20260515/pixi.toml

TMP_BASE=${SLURM_TMPDIR:-/tmp}
WORK="$TMP_BASE/msv_v2_${SLURM_JOB_ID:-$$}_${REGION_ID}_${PHENO}"
mkdir -p "$WORK"
trap "rm -rf $WORK" EXIT

mkdir -p "$OUT_DIR"

DIAG_PATH="$OUT_DIR/${REGION_ID}__${PHENO}__diag.tsv"
if [ -f "$DIAG_PATH" ]; then
  echo "[$(date)] diag exists, skipping: $DIAG_PATH"
  exit 0
fi

echo "[$(date)] V2 prepare jwbundle for $REGION_ID pheno=$PHENO study=$STUDY"
$PIXI run --manifest-path "$PIXI_MANIFEST" python \
    "$BASE/code/02_prepare_jwbundle_multisv/prepare_jwbundle_multisv.py" \
    --bundle "$BUNDLE" \
    --pheno "$PHENO" \
    --study "$STUDY" \
    --out-dir "$WORK"

JW_PATH="$WORK/${REGION_ID}__${PHENO}__jw.npz"
if [ ! -f "$JW_PATH" ]; then
  echo "ERROR: jwbundle not produced at $JW_PATH"
  exit 1
fi

echo "[$(date)] V2 finemap susie_rss (L=5, max_iter=500) $REGION_ID pheno=$PHENO"
export V2_CANONICAL_QC_R_DIR=${CENSV_ROOT}/analyses/sv_catalog/code/q8_validation/qc/R
$PIXI run --manifest-path "$PIXI_MANIFEST" Rscript \
    "$BASE/code/pipeline/03_finemap/R/run_susie_rss_multisv.R" \
    "$JW_PATH" "$OUT_DIR" "$PHENO" "$REGION_ID"

echo "[$(date)] DONE V2 $REGION_ID pheno=$PHENO"
ls -la "$DIAG_PATH"
