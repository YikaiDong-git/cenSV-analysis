#!/bin/bash
#SBATCH --job-name=chr10_qtl1e6
#SBATCH --array=0-1
#SBATCH --nodes=1
#SBATCH --mem=8G
#SBATCH -c 2
#SBATCH --time=2:00:00
#SBATCH --output=${CENSV_ROOT}/analyses/cdr_qtl/logs/chr10_qtl1e6_%A_%a.log
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail
P=${CENSV_ROOT}
WU=$P/analyses/cdr_qtl
PHENOS=(cdr_pos_offset cdr_pos_iqr)
PH=${PHENOS[$SLURM_ARRAY_TASK_ID]}
cd "$P"
echo "[$(date +%T)] chr10 $PH @1e6"
~/.pixi/bin/pixi run python3 "$WU/code/scan/b6_cdr_pos_qtl_perchrom.py" \
  --chrom chr10 --phenotype "$PH" --n-perm 1000000 --out-dir results/cdr_qtl_corrected_1e6
echo "[$(date +%T)] DONE chr10 $PH @1e6"
