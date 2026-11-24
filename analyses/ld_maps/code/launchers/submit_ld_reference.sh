#!/bin/bash
: "${CENSV_COLLAB_ROOT:?set CENSV_COLLAB_ROOT (see README)}"
: "${CENSV_LAB:?set CENSV_LAB (see README)}"
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${HOME:?set HOME (see README)}"
set -euo pipefail
P=${CENSV_ROOT}
DISC=$P/analyses/ld_maps
FRAMES=$DISC/output/_censv_frames
OUT=$DISC/output/_ld_reference_g38
LOGS=$DISC/logs/ldrefg38
SV_BASE=$P/sv_extraction/decompose/mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200
KG=${CENSV_LAB}/shared/1KG
ACTIVE_BED=${CENSV_COLLAB_ROOT}/resources/chm13v2.0_censat_v2.1.active.bed
SCRIPT=$DISC/code/plot_sv_ld_pubfig.py
PIXI_PY="${HOME}/.pixi/bin/pixi run --manifest-path $P/pixi.toml python"

mkdir -p "$OUT" "$LOGS"
n=0
for chrom in $(cat "$FRAMES/chroms.txt"); do
    sj="$FRAMES/$chrom/summary.json"
    [ -f "$sj" ] || { echo "[skip] $chrom - no summary.json"; continue; }
    vcf="$KG/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_${chrom}.recalibrated_variants.vcf.gz"
    [ -e "$vcf" ] || { echo "[skip] $chrom - no VCF"; continue; }
    for side in L R; do
        region=$(python3 -c "import json,sys; f=json.load(open('$sj'))['flanks'].get('$side') or {}; print(f.get('region') or '')")
        [ -n "$region" ] || { echo "[skip] $chrom $side - no window"; continue; }
        sub="$OUT/$chrom/$side"
        mkdir -p "$sub"
        sbatch --parsable \
            --job-name="ldref_${chrom}_${side}" \
            --output="$LOGS/ldref_${chrom}_${side}_%j.log" \
            --nodes=1 \
            --exclude=node01,node02,node04,node05,node06,node10,node11,node12,node13,node14 \
            --mem=12G -c 2 --time=02:00:00 \
            --wrap="$PIXI_PY $SCRIPT \
                --chr-dir $SV_BASE/$chrom \
                --output-dir $sub \
                --snp-vcf $vcf \
                --snp-region $region \
                --active-bed $ACTIVE_BED \
                --upstream-extension-bp 100000 \
                --cross-max-snp-plot 160 \
                --cross-max-sv-plot 80 \
                --snp-min-maf 0.01 \
                --allow-sv-only-fallback" >/dev/null
        n=$((n + 1))
    done
done
echo "submitted $n reference panels; outputs in $OUT/<chrom>/<side>/"
