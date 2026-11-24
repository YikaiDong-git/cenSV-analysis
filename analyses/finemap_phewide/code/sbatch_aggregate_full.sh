#!/bin/bash
#SBATCH -J ssv_agg -c 2 --mem=32G -t 1:30:00
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"
set -e
NEW=${CENSV_ROOT}/analyses/finemap_phewide
PIXI=${CENSV_WORKSPACE}/_shared/tools/bayesian_finemap/v1_20260515/pixi.toml
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd $NEW
$HOME/.pixi/bin/pixi run --manifest-path $PIXI python code/aggregate_full.py \
    "output/full/*/*/*_diag.tsv" output/full_landscape.txt output/full_catalog.tsv 2>&1 \
    | grep -vE 'WARN|deprecated|project|platforms|replace|╭|│|·|╰|─' || true
touch $NEW/logs/FULL_AGG_DONE
echo "AGGREGATE_FULL DONE $(date +%H:%M:%S)"
