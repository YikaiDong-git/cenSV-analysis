#!/usr/bin/env bash
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"
set -euo pipefail
ROOT=${CENSV_WORKSPACE}
SH=$ROOT/_shared/datasets/centromere_globe_logsdon2025/v_20260611
B=$ROOT/tools/CDR-Finder
D=$ROOT/Project1_Centromere/analyses/cdr_calls
BT=$B/env/.pixi/envs/default/bin/bedtools
GLOBE=$D/inputs/globe_aa
mkdir -p "$GLOBE" "$D/inputs/active_array"

if [ ! -s "$GLOBE/.done" ]; then
  unzip -o -j "$SH/cen_annot_HPRC.zip" '*AS-HOR_stv_row.ort.bed' -d "$GLOBE" >/dev/null
  touch "$GLOBE/.done"
fi

for TAG in "$@"; do
  SAMPLE="${TAG%_*}"; HAP="${TAG##*_}"
  case "$HAP" in hap1) H=1;; hap2) H=2;; *) H="${HAP#hap}";; esac
  cat "$GLOBE"/*.bed \
    | awk -F'\t' -v S="$SAMPLE" -v H="$H" '{n=split($1,a,"_"); key=a[n]; pat="^"S"#"H"#"; if(key ~ pat){split(key,b,":"); print b[1]"\t"$2"\t"$3}}' \
    | sort -k1,1 -k2,2n | $BT merge -i - > "$D/inputs/active_array/${TAG}.bed"
  echo "$TAG: $(wc -l < "$D/inputs/active_array/${TAG}.bed") active-array regions, $(awk '{s+=$3-$2}END{printf "%.1f",s/1e6}' "$D/inputs/active_array/${TAG}.bed") Mb"
done
