#!/bin/bash
#SBATCH --job-name=canon_scan
#SBATCH --output=canon_scan_%A_%a.log
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=6:00:00
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail
WU=${CENSV_ROOT}/analyses/cdr_qtl
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$WU/code/admit_loci/scan_tasks.txt")
CHROM=$(echo "$LINE" | awk '{print $1}')
PHENO=$(echo "$LINE" | awk '{print $2}')
[ -z "$CHROM" ] && { echo "no task for index $SLURM_ARRAY_TASK_ID"; exit 1; }
if [ "$PHENO" = "cdr_pos_offset_deepest" ]; then
  PF="$WU/output/per_hap_chrom_phenotypes_deepest.tsv"
else
  PF="$WU/output/per_hap_chrom_phenotypes.tsv"
fi
echo "[$(date)] canonical scan $CHROM / $PHENO (pheno-file=$(basename $PF))"
cd ${CENSV_ROOT}
~/.pixi/bin/pixi run python "$WU/code/scan/b6_cdr_pos_qtl_perchrom.py" \
  --chrom "$CHROM" --phenotype "$PHENO" \
  --sv-source "$WU/canonical_geno/perchrom" \
  --pheno-file "$PF" \
  --n-perm 100000 \
  --out-dir "$WU/results/canonical_scan_1e5"
echo "[$(date)] DONE $CHROM / $PHENO"
