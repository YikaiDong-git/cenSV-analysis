#!/bin/bash
#SBATCH --job-name=cs_cond_chr3
#SBATCH --mem=32G
#SBATCH -c 4
#SBATCH --nodes=1
#SBATCH --exclude=node01,node02,node03,node04,node09,node10,node11,node12,node13,node14
#SBATCH --time=1:00:00
#SBATCH --output=${CENSV_ROOT}/analyses/finemap_multi_sv/validation/ssimp_impg/work/diag_chr3/cs_condition_chr3.%j.out
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"
set -euo pipefail
P=${CENSV_WORKSPACE}
BASE=$P/Project1_Centromere/analyses/finemap_multi_sv
V=$BASE/validation/ssimp_impg
PIXI=~/.pixi/bin/pixi
PIXI_MANIFEST=$P/_shared/tools/bayesian_finemap/v1_20260515/pixi.toml
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
export V2_CANONICAL_QC_R_DIR=$P/Project1_Centromere/analyses/sv_catalog/code/q8_validation/qc/R
if [ ! -f "$V/work/diag_chr3/diag_chr3_matched.npz" ]; then
  echo "[$(date)] npz missing -> rebuild"; BT=$P/_shared/envs/workspace/v1/.pixi/envs/biotools/bin; $BT/python $V/scripts/diag_chr3_build.py
fi
echo "[$(date)] credible-set conditioning"
$PIXI run --manifest-path "$PIXI_MANIFEST" Rscript $V/scripts/cs_condition_chr3.R
echo "[$(date)] DONE"
