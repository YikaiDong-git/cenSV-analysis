#!/bin/bash
#SBATCH -J ssv_pz -c 1 --mem=8G -t 1:00:00
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"
set -e
NEW=${CENSV_ROOT}/analyses/finemap_phewide
PIXI=${CENSV_WORKSPACE}/_shared/tools/bayesian_finemap/v1_20260515/pixi.toml
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1  # avoid BLAS thrash (impg_fill 5s->0.02s)
cd $NEW
BUNDLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" manifests/region_bundles.txt)
mkdir -p output/zmat
echo "[pz $SLURM_ARRAY_TASK_ID] $(basename $BUNDLE) START $(date +%H:%M:%S)"
$HOME/.pixi/bin/pixi run --manifest-path $PIXI python code/precompute_z.py \
    "$BUNDLE" manifests/traits_all464.tsv output/zmat 2>&1 \
    | grep -vE 'Pixi task|project|platforms|replace|environment|╭|│|·|╰|─' || true
echo "[pz $SLURM_ARRAY_TASK_ID] DONE $(date +%H:%M:%S)"
