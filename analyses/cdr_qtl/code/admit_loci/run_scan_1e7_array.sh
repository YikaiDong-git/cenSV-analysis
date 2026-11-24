#!/bin/bash
#SBATCH --job-name=canon_1e7
#SBATCH --output=canon_1e7_%A_%a.log
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --time=6:00:00
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail
WU=${CENSV_ROOT}/analyses/cdr_qtl
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$WU/results/canonical_admit/canonical_floor_cells.txt")
CHROM=$(echo "$LINE" | awk '{print $1}'); PHENO=$(echo "$LINE" | awk '{print $2}')
[ -z "$CHROM" ] && { echo "no cell for $SLURM_ARRAY_TASK_ID"; exit 1; }
if [ "$PHENO" = "cdr_pos_offset_deepest" ]; then PF="$WU/output/per_hap_chrom_phenotypes_deepest.tsv"; else PF="$WU/output/per_hap_chrom_phenotypes.tsv"; fi
echo "[$(date)] 1e7 canonical $CHROM / $PHENO"
cd ${CENSV_ROOT}
~/.pixi/bin/pixi run python "$WU/code/scan/b6_cdr_pos_qtl_perchrom.py" \
  --chrom "$CHROM" --phenotype "$PHENO" \
  --sv-source "$WU/canonical_geno/perchrom" --pheno-file "$PF" \
  --n-perm 10000000 \
  --out-dir "$WU/results/canonical_scan_1e7"
echo "[$(date)] DONE $CHROM / $PHENO"
