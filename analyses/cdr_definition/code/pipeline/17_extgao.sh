#!/bin/bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
#SBATCH --job-name=gao_extgao
#SBATCH --mem=24G
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --output=${CENSV_ROOT}/analyses/cdr_definition/logs/gao_extgao_%j.log
set -euo pipefail
export OMP_NUM_THREADS=1

PROJ=${CENSV_ROOT}
U=$PROJ/analyses/cdr_definition
WU=$PROJ/analyses/cdr_qtl_full_panel
CANON=$PROJ/analyses/cdr_qtl
cd "$PROJ"

echo "[$(date)] === phenotypes on extent-corrected + 150-kb-merged calls ==="
~/.pixi/bin/pixi run python "$U/code/pipeline/03_recompute_phenotypes.py" \
  --dips-root "$U/output/callsets/ext_gao_cdrs" \
  --out-name per_hap_chrom_phenotypes_extgao.tsv

echo "[$(date)] === compat shim for the extgao arm ==="
~/.pixi/bin/pixi run python "$U/code/pipeline/17b_compat_extgao.py"

echo "[$(date)] === chr8 scan, extgao arm ==="
COMPAT=$U/output/compat/compat_extgao_splice
for PHENO in cdr_pos_offset cdr_pos_offset_deepest mean_5mC; do
  if [ "$PHENO" = "cdr_pos_offset_deepest" ]; then
    PF=$COMPAT/per_hap_chrom_phenotypes_deepest.tsv
  else
    PF=$COMPAT/per_hap_chrom_phenotypes.tsv
  fi
  ~/.pixi/bin/pixi run python "$CANON/code/scan/b6_cdr_pos_qtl_perchrom.py" \
    --chrom chr8 --phenotype "$PHENO" \
    --sv-source "$WU/geno/perchrom" \
    --pheno-file "$PF" \
    --n-perm 100000 --seed 42 \
    --extra-cols sv_locus,geno_group,status,cv_r2,strategy,arm \
    --out-dir "$U/output/scan/scan_h2h/extgao"
done
echo "[$(date)] EXTGAO_DONE"
