#!/usr/bin/env bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"
set -uo pipefail
P=${CENSV_ROOT}/subproject_methylation
D=${CENSV_ROOT}/analyses/cdr_calls
B=${CENSV_WORKSPACE}/tools/CDR-Finder
BT=$B/env/.pixi/envs/default/bin/bedtools
cd $P
SAMPLES="HG00733 HG02818 NA19036 NA19240"
for SAMPLE in $SAMPLES; do for HAP in hap1 hap2; do
  TAG=${SAMPLE}_${HAP}; AH=raw/censat/${TAG}_hprc_r2_v1.as_hor.bed
  if [ ! -s "$AH" ]; then
    URL=$(awk -F'\t' -v s=$SAMPLE -v h=$HAP '$1==s&&$2==h{print $13;exit}' manifest/methylation_manifest.tsv)
    [ -n "$URL" ] && ~/.pixi/bin/pixi run aws s3 cp --no-sign-request --quiet "$URL" "$AH"
  fi
  [ -s "$AH" ] || { echo "no as_hor $TAG"; continue; }
  cut -f1-3 "$AH" | sort -k1,1 -k2,2n | $BT merge -d 1000 -i - | awk '($3-$2)>=100000' > $D/inputs/active_array/${TAG}.bed
  echo "$TAG as_hor active-array: $(wc -l < $D/inputs/active_array/${TAG}.bed) regions"
done; done
for SAMPLE in $SAMPLES; do for HAP in hap1 hap2; do
  bash $D/code/run_cdr_one.sh $SAMPLE $HAP || echo "[recover-fail] ${SAMPLE}_${HAP}"
done; done
echo RECOVER_DONE
