#!/bin/bash
#SBATCH --job-name=msv_bcac
#SBATCH --mem=64G
#SBATCH -c 4
#SBATCH --time=2:00:00
#SBATCH --output=${CENSV_ROOT}/analyses/finemap_multi_sv/logs/full_restart/bcac/bcac.%j.out
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"

set -euo pipefail

BASE=${CENSV_ROOT}/analyses/finemap_multi_sv
PIXI=~/.pixi/bin/pixi
PIXI_MANIFEST=${CENSV_WORKSPACE}/_shared/tools/bayesian_finemap/v1_20260515/pixi.toml

cd "$BASE"
mkdir -p logs/full_restart/bcac output/full_restart/bcac_replication

PROCESSED=$BASE/inputs/bcac_2017/processed_b38/bcac_erneg_chr1_b38.tsv.gz
if [ ! -f "$PROCESSED" ]; then
  echo "[$(date)] Phase C step 1: extract + liftover BCAC chr1 b37→b38"
  $PIXI run --manifest-path "$PIXI_MANIFEST" python \
      "$BASE/code/pipeline/06_bcac_replication/extract_bcac_erneg_chr1_b38.py"
else
  echo "[$(date)] BCAC b38 TSV already exists: $PROCESSED"
fi

JW_PATH=$BASE/output/full_restart/bcac_replication/chr1_super1_118296444_148907004__bcac_erneg__jw.npz
if [ ! -f "$JW_PATH" ]; then
  echo "[$(date)] Phase C step 2: prepare BCAC jwbundle"
  $PIXI run --manifest-path "$PIXI_MANIFEST" python \
      "$BASE/code/pipeline/06_bcac_replication/prepare_jwbundle_bcac.py"
else
  echo "[$(date)] jwbundle already exists"
fi

OUT_DIAG=$BASE/output/full_restart/bcac_replication/chr1_super1_118296444_148907004__bcac_erneg__diag.tsv
if [ ! -f "$OUT_DIAG" ]; then
  echo "[$(date)] Phase C step 3: finemap susie_rss L=5 max_iter=500"
  export V2_CANONICAL_QC_R_DIR=${CENSV_ROOT}/analyses/sv_catalog/code/q8_validation/qc/R
  $PIXI run --manifest-path "$PIXI_MANIFEST" Rscript \
      "$BASE/code/pipeline/03_finemap/R/run_susie_rss_multisv.R" \
      "$JW_PATH" "$BASE/output/full_restart/bcac_replication" \
      "bcac_erneg" "chr1_super1_118296444_148907004"
else
  echo "[$(date)] diag already exists"
fi

echo "[$(date)] Phase C DONE. chr1_sv616 result:"
grep -E "chr1_sv616" "$OUT_DIAG" 2>/dev/null | head -5
