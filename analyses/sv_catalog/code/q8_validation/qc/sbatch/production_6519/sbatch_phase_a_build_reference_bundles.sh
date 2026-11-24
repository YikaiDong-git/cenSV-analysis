#!/bin/bash
#SBATCH --job-name=p6519_pA
#SBATCH --output=${CENSV_ROOT}/analyses/sv_catalog/logs/qc/production_6519/phaseA_%A_%a.log
#SBATCH --error=${CENSV_ROOT}/analyses/sv_catalog/logs/qc/production_6519/phaseA_%A_%a.err
#SBATCH --time=3:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=2
#SBATCH --array=1-326%100
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"
#SBATCH --exclude=node01,node02,node04,node05,node08,node09,node10,node11,node12,node13

set -euo pipefail
sleep $((RANDOM % 15))

PROJ=${CENSV_ROOT}/analyses/sv_catalog
PROJROOT=${CENSV_ROOT}
CODE_DIR=$PROJ/code/ref_compare_proxy_vs_dense
PIXI=${CENSV_WORKSPACE}/_shared/tools/bayesian_finemap/v1_20260515
PY=$PIXI/.pixi/envs/default/bin/python

MANIFEST=$PROJ/output/ref_compare_proxy_vs_dense/sweep_6519_locked_v3/manifest/MASTER_MANIFEST_6519.tsv
BUNDLE_ROOT=$PROJ/output/ref_compare_proxy_vs_dense/jwbundle_v3_locked

K=20
BASE=$(( (SLURM_ARRAY_TASK_ID - 1) * K ))

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

export PATH=${CENSV_ROOT}/.pixi/envs/default/bin:$PATH

cd $PIXI

n_built=0
n_skip=0
HEARTBEAT_FILE=$PROJ/logs/qc/production_6519/STATUS_worker_heartbeat
mkdir -p $(dirname "$HEARTBEAT_FILE")
echo "=========================================="
echo "[$(date '+%F %T')] phaseA jid=${SLURM_ARRAY_JOB_ID} task=${SLURM_ARRAY_TASK_ID} node=$(hostname)"
echo "[phaseA] K=$K  BASE=$BASE  range=$((BASE+1)).. $((BASE+K))"
echo "=========================================="

for i in $(seq 0 $((K - 1))); do
  PAIR_ID=$((BASE + i + 1))
  LINE=$(awk -F'\t' -v r=$((PAIR_ID + 1)) 'NR==r' "$MANIFEST")
  if [ -z "$LINE" ]; then
    echo "[phaseA] pair_id=$PAIR_ID past EOF; stopping inner loop"
    break
  fi

  SV_LOCUS=$(echo "$LINE"   | cut -f2)
  SV_ID=$(echo "$LINE"      | cut -f3)
  SV_CHROM=$(echo "$LINE"   | cut -f4)
  SVPOS=$(echo "$LINE"      | cut -f5)
  STRATEGY=$(echo "$LINE"   | cut -f6)
  ARM=$(echo "$LINE"        | cut -f7)
  WIN_START=$(echo "$LINE"  | cut -f8)
  WIN_END=$(echo "$LINE"    | cut -f9)
  BUNDLE_REL=$(echo "$LINE" | cut -f10)

  OUT_DIR=$BUNDLE_ROOT/$BUNDLE_REL

  if [ -f "$OUT_DIR/reference_bundle.npz.sha256" ] && \
     [ -s "$OUT_DIR/reference_bundle.npz" ] && \
     [ -f "$OUT_DIR/manifest.json" ]; then
    n_skip=$((n_skip + 1))
    echo "[skip] pair_id=$PAIR_ID  $SV_LOCUS  $ARM  (bundle complete)"
    continue
  fi
  if [ -d "$OUT_DIR" ]; then
    rm -f "$OUT_DIR/reference_bundle.npz" \
          "$OUT_DIR/reference_bundle.npz.sha256" \
          "$OUT_DIR/manifest.json" \
          "$OUT_DIR/snp_meta.tsv.gz"
  fi

  mkdir -p "$OUT_DIR"
  echo "[build] pair_id=$PAIR_ID  $SV_LOCUS  $ARM  strategy=$STRATEGY  win=[$WIN_START,$WIN_END]"

  $PY $CODE_DIR/build_our_bundle_generic_resid.py \
      --sv-id "$SV_ID" \
      --sv-chrom "$SV_CHROM" \
      --svpos "$SVPOS" \
      --imputation-strategy "$STRATEGY" \
      --local-window-start "$WIN_START" \
      --local-window-end "$WIN_END" \
      --matched-ens "" \
      --qtl-type eqtl \
      --out-dir "$OUT_DIR" \
      --maf-min 0.01

  if [ ! -s "$OUT_DIR/reference_bundle.npz" ] || \
     [ ! -f "$OUT_DIR/reference_bundle.npz.sha256" ] || \
     [ ! -f "$OUT_DIR/manifest.json" ] || \
     [ ! -s "$OUT_DIR/snp_meta.tsv.gz" ]; then
    echo "[FATAL] pair_id=$PAIR_ID build incomplete; sidecar/npz/manifest/snp_meta missing"
    exit 3
  fi
  n_built=$((n_built + 1))
  touch "$HEARTBEAT_FILE"
done

echo "=========================================="
echo "[$(date '+%F %T')] phaseA task=$SLURM_ARRAY_TASK_ID DONE  built=$n_built  skipped=$n_skip"
echo "=========================================="
