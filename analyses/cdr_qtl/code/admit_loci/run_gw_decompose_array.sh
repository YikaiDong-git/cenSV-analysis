#!/bin/bash
#SBATCH --job-name=gw_decomp
#SBATCH --output=gw_decomp_%A_%a.log
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --time=6:00:00
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail
WU=${CENSV_ROOT}/analyses/cdr_qtl
CHROM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$WU/code/admit_loci/gw_chrom_list.txt")
[ -z "$CHROM" ] && { echo "no chrom for index $SLURM_ARRAY_TASK_ID"; exit 1; }
exec bash "$WU/code/admit_loci/run_decompose.sh" "$CHROM"
