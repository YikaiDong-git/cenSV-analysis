#!/bin/bash
#SBATCH --job-name=scan0p1
#SBATCH --mem=48G
#SBATCH -c 2
#SBATCH --time=8:00:00
#SBATCH --output=${CENSV_SCRATCH}/ssv_joint/logs/scan0p1.%A_%a.out
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_SCRATCH:?set CENSV_SCRATCH (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"
set -euo pipefail
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2

BASE=${CENSV_ROOT}/analyses/finemap_multi_sv
P10=$BASE/code/pipeline/10_block_multisv_464trait_scan
FAN=${CENSV_SCRATCH}/ssv_joint
BUILDALL=$P10/code/build_all_regions_drop_srss.py
DRIVER=$P10/code/ctwas_group_em_srss.R
TABIX_MANIFEST=$P10/inputs/tabix_manifest_464.tsv
REGION_LIST=$P10/inputs/region_list.tsv
TRAIT_MANIFEST=$P10/trait_manifest_464.tsv
BUNDLES_ROOT=$FAN/bundles_0p1
PIXI=~/.pixi/bin/pixi
PIXI_MANIFEST=${CENSV_WORKSPACE}/_shared/tools/bayesian_finemap/v1_20260515/pixi.toml

mkdir -p $FAN/logs $FAN/scan_0p1

ROW=$(awk -v i=$((SLURM_ARRAY_TASK_ID + 1)) 'NR==i' "$TRAIT_MANIFEST")
PHENO=$(echo "$ROW" | awk -F'\t' '{print $1}')
STUDY=$(echo "$ROW" | awk -F'\t' '{print $2}')
echo "[$(date)] scan0p1 task=$SLURM_ARRAY_TASK_ID pheno=$PHENO study=$STUDY node=$(hostname)"
if [ -z "$PHENO" ] || [ -z "$STUDY" ]; then echo "ERROR empty pheno/study"; exit 1; fi

OUT_DIR="$FAN/scan_0p1/$STUDY"
mkdir -p "$OUT_DIR"
TRACE="$OUT_DIR/${PHENO}__param_trace.tsv"
if [ -f "$TRACE" ]; then echo "[$(date)] trace exists for $PHENO, skipping"; exit 0; fi

TMP_BASE=${SLURM_TMPDIR:-/tmp}
WORK="$TMP_BASE/scan0p1_${SLURM_JOB_ID:-$$}_${PHENO}"
mkdir -p "$WORK"
trap "rm -rf $WORK" EXIT
SKIPLOG="$OUT_DIR/${PHENO}__skipped_regions.tsv"

echo "[$(date)] build-all 21 regions (cv2>=0.1) for $PHENO"
$PIXI run --manifest-path "$PIXI_MANIFEST" env OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 VECLIB_MAXIMUM_THREADS=2 python "$BUILDALL" \
  --region-list "$REGION_LIST" --bundles-root "$BUNDLES_ROOT" \
  --pheno "$PHENO" --study "$STUDY" --out-dir "$WORK" \
  --tabix-manifest "$TABIX_MANIFEST" --skip-log "$SKIPLOG"

for f in "$WORK"/*__jw_drop.npz; do mv "$f" "${f%__jw_drop.npz}__jw.npz"; done

echo "[$(date)] group-EM + estimate_s_rss (niter=30 L=10) for $PHENO"
$PIXI run --manifest-path "$PIXI_MANIFEST" env OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 VECLIB_MAXIMUM_THREADS=2 Rscript "$DRIVER" "$WORK" "$OUT_DIR" "$PHENO" 30 10

echo "[$(date)] DONE scan0p1 $PHENO"
ls -la "$TRACE" 2>/dev/null || echo "NO TRACE WRITTEN"
