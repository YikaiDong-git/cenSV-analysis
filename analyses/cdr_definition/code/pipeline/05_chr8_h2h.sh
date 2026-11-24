#!/bin/bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
#SBATCH --job-name=gao_chr8_h2h
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --output=${CENSV_ROOT}/analyses/cdr_definition/logs/gao_chr8_h2h_%j.log
set -euo pipefail
export OMP_NUM_THREADS=1

PROJ=${CENSV_ROOT}
U=$PROJ/analyses/cdr_definition
WU=$PROJ/analyses/cdr_qtl_full_panel
CANON=$PROJ/analyses/cdr_qtl
NPERM=${NPERM:-100000}
CHROM=${CHROM:-chr8}

cd "$PROJ"
for ARM in dip gao; do
  COMPAT=$U/output/compat/compat_${ARM}_splice
  for PHENO in cdr_pos_offset cdr_pos_offset_deepest mean_5mC; do
    if [ "$PHENO" = "cdr_pos_offset_deepest" ]; then
      PF=$COMPAT/per_hap_chrom_phenotypes_deepest.tsv
    else
      PF=$COMPAT/per_hap_chrom_phenotypes.tsv
    fi
    OUTDIR=$U/output/scan/scan_h2h/$ARM
    mkdir -p "$OUTDIR"
    echo "[$(date)] === $ARM / $CHROM / $PHENO  (nperm=$NPERM) ==="
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
echo "[$(date)] H2H_DONE"
