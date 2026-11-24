#!/bin/bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail
PROJ=${CENSV_ROOT}
cd $PROJ
SUFFIX="mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
VCF="variant_calling/svim_asm/cohort_svim_asm_${SUFFIX}.vcf.gz"
OUT_BASE="sv_extraction/decompose/${SUFFIX}"
PIXIRUN="~/.pixi/bin/pixi run"
CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
mkdir -p logs/decompose
N=0
mkdir -p ${OUT_BASE}
for ch in "${CHROMS[@]}"; do
  
  sbatch --mem=16G -c 4 --time=4:00:00 \
    --job-name=dec_chr${ch} \
    --output=logs/decompose/dec_chr${ch}_%A.log \
    --wrap="~/.pixi/bin/pixi run python3 workflow/scripts/haplogroup_clustering.py \
      --vcf ${VCF} \
      --output ${OUT_BASE} \
      --chrom chr${ch} \
      --method cosine_nj \
      --min-samples 2 \
      --k-range 2,30 \
      --decompose"
  N=$((N+1))
done
echo "Submitted $N decomposition jobs"
