#!/bin/bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${HOME:?set HOME (see README)}"
set -euo pipefail
P=${CENSV_ROOT}
DISC=$P/analyses/ld_maps
FRAMES=$DISC/output/_censv_frames
OUT=$DISC/output/_ld_zoom_tiers
LOGS=$DISC/logs/ldzoom
SCRIPT=$DISC/code/measure_ld_zoom_tiers.py
PIXI_PY="${HOME}/.pixi/bin/pixi run --manifest-path $P/pixi.toml python"

mkdir -p "$OUT" "$LOGS"
n=0
for chrom in $(cat "$FRAMES/chroms.txt"); do
    [ -f "$FRAMES/$chrom/summary.json" ] || { echo "[skip] $chrom - no summary.json"; continue; }
    sbatch --parsable \
        --job-name="ldzoom_${chrom}" \
        --output="$LOGS/ldzoom_${chrom}_%j.log" \
        --nodes=1 \
        --exclude=node01,node02,node04,node05,node08,node09,node10,node13 \
        --mem=8G -c 1 --time=03:00:00 \
        --wrap="$PIXI_PY $SCRIPT --chrom $chrom --out $OUT" >/dev/null
    n=$((n + 1))
done
echo "submitted $n tier measurements; outputs in $OUT/<chrom>_tiers.json"
