#!/bin/bash
: "${CENSV_COLLAB_ROOT:?set CENSV_COLLAB_ROOT (see README)}"
: "${CENSV_LAB:?set CENSV_LAB (see README)}"
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${HOME:?set HOME (see README)}"

set -euo pipefail

WINDOW_LABEL=${1:-100kb}
WINDOW_BP=${2:-100000}

DISC_DIR=${CENSV_ROOT}/analyses/ld_maps
SCRIPT="$DISC_DIR/code/plot_sv_ld_pubfig.py"
LOGS_DIR="$DISC_DIR/logs"
OUT_BASE="$DISC_DIR/output"

SV_BASE=${CENSV_ROOT}/sv_extraction/decompose/mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200
KG_BASE=${CENSV_LAB}/shared/1KG
ACTIVE_BED=${CENSV_COLLAB_ROOT}/resources/chm13v2.0_censat_v2.1.active.bed
PIXI_PY="${HOME}/.pixi/bin/pixi run --manifest-path ${CENSV_ROOT}/pixi.toml python"

mkdir -p "$LOGS_DIR" "$OUT_BASE"

echo "submitting plot_sv_ld_pubfig — WINDOW_LABEL=$WINDOW_LABEL WINDOW_BP=$WINDOW_BP"
printf "%-6s %-12s %-12s %-30s\n" chrom hor_start hor_end snp_region

CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX)

for chrom in "${CHROMS[@]}"; do
    read -r hor_start hor_end < <(
        awk -v c="$chrom" 'BEGIN{best=0} $1==c { len=$3-$2; if (len > best){ best=len; s=$2; e=$3 } } END{ print s, e }' "$ACTIVE_BED"
    )
    if [ -z "$hor_start" ] || [ -z "$hor_end" ]; then
        echo "[skip] $chrom — no active HOR span in BED"
        continue
    fi
    snp_start=$(( hor_start - WINDOW_BP ))
    snp_end=$(( hor_start - 1 ))
    if [ "$snp_start" -lt 1 ]; then snp_start=1; fi
    snp_region="${chrom}:${snp_start}-${snp_end}"

    chr_dir="${SV_BASE}/${chrom}"
    snp_vcf="${KG_BASE}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_${chrom}.recalibrated_variants.vcf.gz"
    if [ ! -e "$snp_vcf" ]; then
        echo "[skip] $chrom — 1KG VCF missing: $snp_vcf"
        continue
    fi
    out_sub="${OUT_BASE}/${WINDOW_LABEL}/${chrom}"
    mkdir -p "$out_sub"

    printf "%-6s %-12s %-12s %-30s\n" "$chrom" "$hor_start" "$hor_end" "$snp_region"

    job_name="ldpub_${WINDOW_LABEL}_${chrom}"
    log_file="${LOGS_DIR}/${job_name}_%j.log"

    sbatch \
        --job-name="$job_name" \
        --output="$log_file" \
        --mem=12G -c 2 --time=02:00:00 \
        --wrap="$PIXI_PY $SCRIPT \
            --chr-dir $chr_dir \
            --output-dir $out_sub \
            --snp-vcf $snp_vcf \
            --snp-region $snp_region \
            --active-bed $ACTIVE_BED \
            --upstream-extension-bp $WINDOW_BP \
            --cross-max-snp-plot 160 \
            --cross-max-sv-plot 80 \
            --snp-min-maf 0.01 \
            --allow-sv-only-fallback"
done

echo
echo "submission done; logs in $LOGS_DIR"
echo "outputs land in $OUT_BASE/$WINDOW_LABEL/<chrom>/"
