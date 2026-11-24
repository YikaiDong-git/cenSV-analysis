#!/usr/bin/env bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -uo pipefail
P=${CENSV_ROOT}/subproject_methylation
D=${CENSV_ROOT}/analyses/cdr_calls
cd $P
n=0; got=0
while read SAMPLE HAP; do
  TAG=${SAMPLE}_${HAP}
  DEST=raw/censat/${TAG}_hprc_r2_v1.active.centromeres.bed
  [ -s "$DEST" ] && continue
  URL=$(awk -F'\t' -v s=$SAMPLE -v h=$HAP '$1==s && $2==h{print $12; exit}' manifest/methylation_manifest.tsv)
  [ -n "$URL" ] || { echo "no URL $TAG"; continue; }
  n=$((n+1))
  ~/.pixi/bin/pixi run aws s3 cp --no-sign-request --quiet "$URL" "$DEST" && got=$((got+1))
done < $D/code/all_haps.txt
echo "attempted=$n downloaded=$got"
ls raw/censat/*.active.centromeres.bed | wc -l
