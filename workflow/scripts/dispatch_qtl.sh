#!/bin/bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail
PROJ=${CENSV_ROOT}
cd $PROJ
SUFFIX="mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
EN_STRAT=best8
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
IMP_BASE=results/imputation/decompose/cosine_nj/${EN_STRAT}_${SUFFIX}
EQTL_BASE=results/eqtl_overlap/decompose/cosine_nj/${EN_STRAT}_${SUFFIX}
SQTL_BASE=results/sqtl_overlap/decompose/cosine_nj/${EN_STRAT}_${SUFFIX}

mkdir -p logs/qtl
N=0

for ch in "${CHROMS[@]}"; do
  IMP_DIR=${IMP_BASE}/chr${ch}
  OTHER=$(IFS=,; chs=(); for c in "${CHROMS[@]}"; do [[ "$c" != "$ch" ]] && chs+=("chr$c"); done; echo "${chs[*]}")
  OTHER=$(echo "$OTHER" | tr ' ' ',')

  for W in "${WINDOWS[@]}"; do
    OUT=${EQTL_BASE}/chr${ch}/${COV_STRING}${EQTL_PEER}/window_${W}
    mkdir -p $OUT
    sbatch --mem=16G -c 1 --time=2:00:00 \
      --job-name=ceqtl_chr${ch}_w${W} \
      --output=logs/qtl/ceqtl_chr${ch}_w${W}_%A.log \
      --wrap="~/.pixi/bin/pixi run python3 workflow/scripts/qtl_overlap.py \
        --imputation-dir $IMP_DIR --centromere-bed $CEN_BED --output $OUT \
        --chromosome chr${ch} --qtl-type eqtl \
        --expression $EXPR --eqtl-covariates $EQTL_COV \
        --regmode within-chrom --regmode-para $W"
    N=$((N+1))
  done

  for W in "${WINDOWS[@]}"; do
    OUT=${SQTL_BASE}/chr${ch}/${COV_STRING}${SQTL_PEER}/window_${W}
    mkdir -p $OUT
    sbatch --mem=16G -c 1 --time=2:00:00 \
      --job-name=csqtl_chr${ch}_w${W} \
      --output=logs/qtl/csqtl_chr${ch}_w${W}_%A.log \
      --wrap="~/.pixi/bin/pixi run python3 workflow/scripts/qtl_overlap.py \
        --imputation-dir $IMP_DIR --centromere-bed $CEN_BED --output $OUT \
        --chromosome chr${ch} --qtl-type sqtl \
        --splicing-bed $SPLICE --phenotype-groups $PHENO_GRP --sqtl-covariates $SQTL_COV \
        --regmode within-chrom --regmode-para $W"
    N=$((N+1))
  done

  OUT=${EQTL_BASE}/chr${ch}/${COV_STRING}${EQTL_PEER}/trans
  mkdir -p $OUT
  sbatch --mem=32G -c 1 --time=6:00:00 \
    --job-name=teqtl_chr${ch} \
    --output=logs/qtl/teqtl_chr${ch}_%A.log \
    --wrap="~/.pixi/bin/pixi run python3 workflow/scripts/qtl_overlap.py \
      --imputation-dir $IMP_DIR --centromere-bed $CEN_BED --output $OUT \
      --chromosome chr${ch} --qtl-type eqtl \
      --expression $EXPR --eqtl-covariates $EQTL_COV \
      --regmode cross-chrom --regmode-para $OTHER"
  N=$((N+1))

  OUT=${SQTL_BASE}/chr${ch}/${COV_STRING}${SQTL_PEER}/trans
  mkdir -p $OUT
  sbatch --mem=32G -c 1 --time=6:00:00 \
    --job-name=tsqtl_chr${ch} \
    --output=logs/qtl/tsqtl_chr${ch}_%A.log \
    --wrap="~/.pixi/bin/pixi run python3 workflow/scripts/qtl_overlap.py \
      --imputation-dir $IMP_DIR --centromere-bed $CEN_BED --output $OUT \
      --chromosome chr${ch} --qtl-type sqtl \
      --splicing-bed $SPLICE --phenotype-groups $PHENO_GRP --sqtl-covariates $SQTL_COV \
      --regmode cross-chrom --regmode-para $OTHER"
  N=$((N+1))
done

echo "dispatched $N QTL jobs (expected 138)"
