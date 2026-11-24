#!/bin/bash
#SBATCH -J ssv_val -c 1 --mem=4G -t 0:30:00
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"
set -e
NEW=${CENSV_ROOT}/analyses/finemap_phewide
BUN=${CENSV_ROOT}/analyses/finemap_multi_sv/output/bundles/chr16/chr16_super1_31371622_48974015_bundle.npz
ZM=$NEW/output/zmat/chr16_super1_31371622_48974015__zmat.npz
PIXI=${CENSV_WORKSPACE}/_shared/tools/bayesian_finemap/v1_20260515/pixi.toml
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd $NEW
for arm in impute matched; do
  $HOME/.pixi/bin/pixi run --manifest-path $PIXI Rscript code/fused_run.R \
    "$BUN" "$ZM" hitbundle $arm 74 78 chr16_super1_31371622_48974015 output/validate_chr16_$arm.tsv 2>&1 \
    | grep -vE 'WARN|deprecated|project|platforms|replace|╭|│|·|╰|─' || true
done
echo "VALIDATE DONE impute=$(tail -n +2 output/validate_chr16_impute.tsv 2>/dev/null|wc -l) matched=$(tail -n +2 output/validate_chr16_matched.tsv 2>/dev/null|wc -l)"
