#!/bin/bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
#SBATCH --job-name=gao_gw
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --output=${CENSV_ROOT}/analyses/cdr_definition/logs/gw/gao_gw_%A_%a.log
set -euo pipefail
export OMP_NUM_THREADS=1

PROJ=${CENSV_ROOT}
U=$PROJ/analyses/cdr_definition
WU=$PROJ/analyses/cdr_qtl_full_panel
CANON=$PROJ/analyses/cdr_qtl
NPERM=${NPERM:-100000}
TASKFILE=${TASKFILE:-$U/output/tasks/gw_chroms.txt}
OUTBASE=${OUTBASE:-$U/output/scan/scan_gw_1e5}
PHENOS=${PHENOS:-"cdr_pos_offset cdr_pos_offset_deepest n_cdrs cdrs_per_mb cdr_fraction median_cdr_width"}

CHROM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TASKFILE" | awk '{print $1}')
[ -z "$CHROM" ] && { echo "no chromosome at task index ${SLURM_ARRAY_TASK_ID}"; exit 1; }

cd "$PROJ"
for ARM in dip gao; do
  COMPAT=$U/output/compat/compat_${ARM}_splice
  for PHENO in $PHENOS; do
    if [ "$PHENO" = "cdr_pos_offset_deepest" ]; then
      PF=$COMPAT/per_hap_chrom_phenotypes_deepest.tsv
    else
      PF=$COMPAT/per_hap_chrom_phenotypes.tsv
    fi
    [ -s "$PF" ] || { echo "missing phenotype file: $PF"; exit 1; }
    OUTDIR=$OUTBASE/$ARM
    mkdir -p "$OUTDIR"
    echo "[$(date)] arm=$ARM chrom=$CHROM pheno=$PHENO nperm=$NPERM"
    ~/.pixi/bin/pixi run python "$CANON/code/scan/b6_cdr_pos_qtl_perchrom.py" \
      --chrom "$CHROM" --phenotype "$PHENO" \
      --sv-source "$WU/geno/perchrom" \
      --pheno-file "$PF" \
      --n-perm "$NPERM" \
      --seed 42 \
      --extra-cols sv_locus,geno_group,status,cv_r2,strategy,arm \
      --out-dir "$OUTDIR"
  done
done
echo "[$(date)] CHROM_DONE $CHROM"
