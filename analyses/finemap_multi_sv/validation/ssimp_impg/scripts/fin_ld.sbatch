#!/bin/bash
#SBATCH --job-name=fin_ld_redo
#SBATCH --mem=16G
#SBATCH -c 2
#SBATCH --nodes=1
#SBATCH --exclude=node01,node02,node03,node04,node09,node10,node11,node12,node13,node14
#SBATCH --time=00:20:00
#SBATCH --output=${CENSV_ROOT}/analyses/finemap_multi_sv/validation/ssimp_impg/work/diag_chr3/fin_ld_redo.%j.out
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"
set -euo pipefail
P=${CENSV_WORKSPACE}
V=$P/Project1_Centromere/analyses/finemap_multi_sv/validation/ssimp_impg
BT=$P/_shared/envs/workspace/v1/.pixi/envs/biotools/bin
$BT/python $V/scripts/fin_ld.py
echo DONE
