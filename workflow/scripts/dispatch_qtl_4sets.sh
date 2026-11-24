#!/bin/bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail
PROJ=${CENSV_ROOT}
cd $PROJ
SUFFIX="mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
COV_STRING="sexTrue_gpc5_peer"
EQTL_PEER=60
SQTL_PEER=15
CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
WINDOWS=(10000000 0)
CEN_BED=resources/grch38_centromeres.bed
EXPR=resources/MAGE/expression_bed.bed.gz
EQTL_COV=resources/MAGE/eqtl_covariates.tab.gz
SPLICE=resources/MAGE/splicing_bed.bed.gz
PHENO_GRP=resources/MAGE/splicing_phenotype_groups.txt
SQTL_COV=resources/MAGE/sqtl_covariates.tab.gz

mkdir -p logs/qtl
N=0

declare -A R2=( [cv01_maf005]=0.1 [cv01_maf010]=0.1 [cv03_maf005]=0.3 [cv03_maf010]=0.3 )
declare -A MAF=( [cv01_maf005]=0.05 [cv01_maf010]=0.10 [cv03_maf005]=0.05 [cv03_maf010]=0.10 )
SETS=(cv01_maf005 cv01_maf010 cv03_maf005 cv03_maf010)
ENSEMBLES=(best8 best6)

for SET in "${SETS[@]}"; do
  MIN_R2=${R2[$SET]}; MIN_MAF=${MAF[$SET]}
  for EN_STRAT in "${ENSEMBLES[@]}"; do
    IMP_BASE=results/imputation/decompose/cosine_nj/${EN_STRAT}_${SUFFIX}
    EQTL_BASE=results/eqtl_overlap/decompose/cosine_nj_${SET}/${EN_STRAT}_${SUFFIX}
    SQTL_BASE=results/sqtl_overlap/decompose/cosine_nj_${SET}/${EN_STRAT}_${SUFFIX}
    for ch in "${CHROMS[@]}"; do
      IMP_DIR=${IMP_BASE}/chr${ch}
      OTHER=$(chs=(); for c in "${CHROMS[@]}"; do [[ "$c" != "$ch" ]] && chs+=("chr$c"); done; IFS=","; echo "${chs[*]}")
      for W in "${WINDOWS[@]}"; do
        OUT=${EQTL_BASE}/chr${ch}/${COV_STRING}${EQTL_PEER}/window_${W}
        mkdir -p $OUT
        sbatch --mem=16G -c 1 --time=2:00:00 --job-name=ceq_${SET}_${EN_STRAT}_chr${ch}_w${W} \
          --output=logs/qtl/ceq_${SET}_${EN_STRAT}_chr${ch}_w${W}_%A.log \
          --wrap="~/.pixi/bin/pixi run python3 workflow/scripts/qtl_overlap.py --imputation-dir $IMP_DIR --centromere-bed $CEN_BED --output $OUT --chromosome chr${ch} --qtl-type eqtl --expression $EXPR --eqtl-covariates $EQTL_COV --regmode within-chrom --regmode-para $W --min-r2 $MIN_R2 --min-maf $MIN_MAF" >/dev/null
        N=$((N+1))
        OUT=${SQTL_BASE}/chr${ch}/${COV_STRING}${SQTL_PEER}/window_${W}
        mkdir -p $OUT
        sbatch --mem=16G -c 1 --time=2:00:00 --job-name=csq_${SET}_${EN_STRAT}_chr${ch}_w${W} \
          --output=logs/qtl/csq_${SET}_${EN_STRAT}_chr${ch}_w${W}_%A.log \
          --wrap="~/.pixi/bin/pixi run python3 workflow/scripts/qtl_overlap.py --imputation-dir $IMP_DIR --centromere-bed $CEN_BED --output $OUT --chromosome chr${ch} --qtl-type sqtl --splicing-bed $SPLICE --phenotype-groups $PHENO_GRP --sqtl-covariates $SQTL_COV --regmode within-chrom --regmode-para $W --min-r2 $MIN_R2 --min-maf $MIN_MAF" >/dev/null
        N=$((N+1))
      done
      OUT=${EQTL_BASE}/chr${ch}/${COV_STRING}${EQTL_PEER}/trans
      mkdir -p $OUT
      sbatch --mem=32G -c 1 --time=6:00:00 --job-name=teq_${SET}_${EN_STRAT}_chr${ch} \
        --output=logs/qtl/teq_${SET}_${EN_STRAT}_chr${ch}_%A.log \
        --wrap="~/.pixi/bin/pixi run python3 workflow/scripts/qtl_overlap.py --imputation-dir $IMP_DIR --centromere-bed $CEN_BED --output $OUT --chromosome chr${ch} --qtl-type eqtl --expression $EXPR --eqtl-covariates $EQTL_COV --regmode cross-chrom --regmode-para $OTHER --min-r2 $MIN_R2 --min-maf $MIN_MAF" >/dev/null
      N=$((N+1))
      OUT=${SQTL_BASE}/chr${ch}/${COV_STRING}${SQTL_PEER}/trans
      mkdir -p $OUT
      sbatch --mem=32G -c 1 --time=6:00:00 --job-name=tsq_${SET}_${EN_STRAT}_chr${ch} \
        --output=logs/qtl/tsq_${SET}_${EN_STRAT}_chr${ch}_%A.log \
        --wrap="~/.pixi/bin/pixi run python3 workflow/scripts/qtl_overlap.py --imputation-dir $IMP_DIR --centromere-bed $CEN_BED --output $OUT --chromosome chr${ch} --qtl-type sqtl --splicing-bed $SPLICE --phenotype-groups $PHENO_GRP --sqtl-covariates $SQTL_COV --regmode cross-chrom --regmode-para $OTHER --min-r2 $MIN_R2 --min-maf $MIN_MAF" >/dev/null
      N=$((N+1))
    done
  done
done

echo "dispatched $N QTL jobs (expected 4 sets x 2 ens x 138 = 1104)"
