#!/bin/bash
: "${CENSV_COLLAB_ROOT:?set CENSV_COLLAB_ROOT (see README)}"
: "${CENSV_LAB:?set CENSV_LAB (see README)}"
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${HOME:?set HOME (see README)}"
set -euo pipefail

DISC=${CENSV_ROOT}/analyses/ld_maps
SCRIPT="$DISC/code/plot_sv_ld_styles.py"
LOGS="$DISC/logs"
PIXI_PY="${HOME}/.pixi/bin/pixi run --manifest-path ${CENSV_ROOT}/pixi.toml python"

CHROM=${1:-chr1}
WINDOW_BP=${2:-100000}
case $WINDOW_BP in
    100000)   WLBL=100kb ;;
    1000000)  WLBL=1mb ;;
    10000000) WLBL=10mb ;;
    *)        WLBL="${WINDOW_BP}bp" ;;
esac
OUT="$DISC/output/styles_${CHROM}_${WLBL}"
HOR_START=$(awk -v c="$CHROM" 'BEGIN{best=0} $1==c { len=$3-$2; if (len > best){best=len; s=$2; e=$3 } } END{ print s }' ${CENSV_COLLAB_ROOT}/resources/chm13v2.0_censat_v2.1.active.bed)
SNP_START=$((HOR_START - WINDOW_BP))
SNP_END=$((HOR_START - 1))
SNP_REGION="${CHROM}:${SNP_START}-${SNP_END}"
CHR_DIR=${CENSV_ROOT}/sv_extraction/decompose/mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200/${CHROM}
SNP_VCF=${CENSV_LAB}/shared/1KG/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_${CHROM}.recalibrated_variants.vcf.gz
ACTIVE_BED=${CENSV_COLLAB_ROOT}/resources/chm13v2.0_censat_v2.1.active.bed

mkdir -p "$OUT" "$LOGS"
echo "$CHROM snp_region=$SNP_REGION window_bp=$WINDOW_BP -> $OUT"

for v in v01 v02 v03 v04 v05 v06 v07 v08 v09 v10; do
    job_name="style_${v}_${CHROM}_${WLBL}"
    sbatch \
        --job-name="$job_name" \
        --output="${LOGS}/${job_name}_%j.log" \
        --mem=12G -c 2 --time=00:30:00 \
        --wrap="$PIXI_PY $SCRIPT --style $v \
            --chr-dir $CHR_DIR \
            --output-dir $OUT \
            --snp-vcf $SNP_VCF \
            --snp-region $SNP_REGION \
            --active-bed $ACTIVE_BED \
            --upstream-extension-bp $WINDOW_BP \
            --cross-max-snp-plot 160 \
            --cross-max-sv-plot 80 \
            --snp-min-maf 0.05 \
            --top-snp-fraction 0.20 \
            --allow-sv-only-fallback"
done
echo "10 jobs submitted"
