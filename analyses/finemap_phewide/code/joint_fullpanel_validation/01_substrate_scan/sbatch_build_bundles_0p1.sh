#!/bin/bash
#SBATCH --job-name=msv_build_0p1
#SBATCH --array=1-21%8
#SBATCH --mem=96G
#SBATCH -c 4
#SBATCH --time=3:00:00
#SBATCH --output=${CENSV_SCRATCH}/ssv_joint/logs/build0p1.%A_%a.out
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_SCRATCH:?set CENSV_SCRATCH (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"
set -euo pipefail
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMEXPR_NUM_THREADS=4

M=${CENSV_ROOT}/analyses/finemap_multi_sv
PIXI=~/.pixi/bin/pixi
PM=${CENSV_WORKSPACE}/_shared/tools/bayesian_finemap/v1_20260515/pixi.toml
REGION_MAP=$M/output/full_restart/manifests/region_map_v2_with_admit.tsv
ADMIT=$M/output/full_restart/admit/admit_decompose_cvr2_0p1.tsv
OUT_BASE=${CENSV_SCRATCH}/ssv_joint/bundles_0p1

ROW=$(awk -v idx=$SLURM_ARRAY_TASK_ID -F'\t' 'NR>1 && $7>=1 {n++; if (n==idx) print $0}' "$REGION_MAP")
if [ -z "$ROW" ]; then echo "no row for task=$SLURM_ARRAY_TASK_ID"; exit 0; fi
CHROM=$(echo "$ROW" | awk -F'\t' '{print $1}')
REGION_ID=$(echo "$ROW" | awk -F'\t' '{print $2}')
RSTART=$(echo "$ROW" | awk -F'\t' '{print $3}')
REND=$(echo "$ROW" | awk -F'\t' '{print $4}')
CHROM_TAG="chr${CHROM}"
OUT_DIR="$OUT_BASE/${CHROM_TAG}"
mkdir -p "$OUT_DIR"
BUNDLE_PATH="$OUT_DIR/${REGION_ID}_bundle.npz"
if [ -f "$BUNDLE_PATH" ]; then echo "[$(date)] bundle_0p1 exists, skip: $BUNDLE_PATH"; exit 0; fi

echo "[$(date)] task=$SLURM_ARRAY_TASK_ID build_0p1 $REGION_ID (chr$CHROM ${RSTART}-${REND}) node=$(hostname)"
$PIXI run --manifest-path "$PM" env OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 python \
    "$M/code/pipeline/02_build_bundles/build_multisv_bundle.py" \
    --chrom "$CHROM_TAG" --region-id "$REGION_ID" \
    --region-start $RSTART --region-end $REND \
    --admit-tsv "$ADMIT" --out-dir "$OUT_DIR"
echo "[$(date)] DONE build_0p1 $REGION_ID"
ls -la "$BUNDLE_PATH"
