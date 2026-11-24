#!/bin/bash
#SBATCH --job-name=full_scan
#SBATCH --output=full_scan_%A_%a.log
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=6:00:00
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail
WU=${CENSV_ROOT}/analyses/cdr_qtl_full_panel
CANON=${CENSV_ROOT}/analyses/cdr_qtl
NPERM=${NPERM:-100000}
OUTDIR=${OUTDIR:-$WU/output/scan_1e5}
TASKFILE=${TASKFILE:-$WU/code/scan/scan_tasks.txt}

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TASKFILE")
CHROM=$(echo "$LINE" | awk '{print $1}')
PHENO=$(echo "$LINE" | awk '{print $2}')
[ -z "$CHROM" ] && { echo "no task for index $SLURM_ARRAY_TASK_ID"; exit 1; }

if [ "$PHENO" = "cdr_pos_offset_deepest" ]; then
  PF="$CANON/output/per_hap_chrom_phenotypes_deepest.tsv"
else
  PF="$CANON/output/per_hap_chrom_phenotypes.tsv"
fi
echo "[$(date)] full-substrate scan $CHROM / $PHENO  nperm=$NPERM  (pheno-file=$(basename $PF))"
cd ${CENSV_ROOT}
~/.pixi/bin/pixi run python "$CANON/code/scan/b6_cdr_pos_qtl_perchrom.py" \
  --chrom "$CHROM" --phenotype "$PHENO" \
  --sv-source "$WU/geno/perchrom" \
  --pheno-file "$PF" \
  --n-perm "$NPERM" \
  --extra-cols sv_locus,geno_group,status,cv_r2,strategy,arm \
  --out-dir "$OUTDIR"
echo "[$(date)] DONE $CHROM / $PHENO"
