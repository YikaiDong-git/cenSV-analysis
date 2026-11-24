#!/bin/bash
#SBATCH --job-name=coloc_bway
#SBATCH --mem=16G
#SBATCH -c 2
#SBATCH --nodes=1
#SBATCH --exclude=node01,node02,node03,node04,node09,node10,node11,node12,node13,node14
#SBATCH --time=00:20:00
#SBATCH --output=${CENSV_ROOT}/analyses/finemap_multi_sv/validation/ssimp_impg/work/coloc_pros1/coloc_broadaway.%j.out
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"
set -euo pipefail
P=${CENSV_WORKSPACE}
BASE=$P/Project1_Centromere/analyses/finemap_multi_sv
V=$BASE/validation/ssimp_impg
BT=$P/_shared/envs/workspace/v1/.pixi/envs/biotools/bin
PIXI=~/.pixi/bin/pixi
RSTATS_MANIFEST=$P/_shared/envs/workspace/v1/pixi.toml
echo "[$(date)] prep broadaway"
$BT/python $V/scripts/coloc_broadaway_prep.py
echo "[$(date)] coloc broadaway"
$PIXI run --manifest-path "$RSTATS_MANIFEST" -e rstats Rscript $V/scripts/coloc_broadaway.R
echo "[$(date)] DONE_COLOC_BROADAWAY"
