#!/usr/bin/env bash
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"
set -euo pipefail
SAMPLE=$1; HAP=$2; TAG="${SAMPLE}_${HAP}"
ROOT=${CENSV_WORKSPACE}
P=$ROOT/Project1_Centromere/subproject_methylation
B=$ROOT/tools/CDR-Finder
D=$ROOT/Project1_Centromere/analyses/cdr_calls
PY=$B/env/.pixi/envs/default/bin/python
BT=$B/env/.pixi/envs/default/bin/bedtools
CW=$P/CDR-Finder/workflow/scripts/calculate_windows.py
CF=$B/upstream/workflow/scripts/cdr_finder.py

BM=$P/raw/bedmethyl_filt/${TAG}.active.bed
OUT=$D/output/calls/${TAG}; mkdir -p "$OUT"
MTMP=/tmp/${TAG}_m_$$.bed

[ -s "$BM" ] || { echo "[FAIL] $TAG no bedmethyl"; exit 1; }
AA="$OUT/aa_union.bed"
cat $D/inputs/active_array/${SAMPLE}_hap1.bed $D/inputs/active_array/${SAMPLE}_hap2.bed 2>/dev/null | sort -k1,1 -k2,2n > "$AA"
[ -s "$AA" ] || { echo "[FAIL] $TAG no globe active-array (either hap)"; exit 1; }

awk -F'\t' '$4=="m"' "$BM" > "$MTMP"
$PY "$CW" --target_bed "$AA" --methylation_tsv "$MTMP" --window_size 5000 -p 1 > "$OUT/binned_freq.bed" 2> "$OUT/calc_windows.log"
rm -f "$MTMP"
cut -f1-4 "$OUT/binned_freq.bed" | sort -k1,1 -k2,2n | $BT intersect -a - -b "$AA" -f 1 -wa -u > "$OUT/binned_aa.bed"
[ -s "$OUT/binned_aa.bed" ] || { echo "[FAIL] $TAG empty binned_aa (no methylation in active array / contig mismatch)"; exit 1; }
$PY "$CF" -i "$OUT/binned_aa.bed" --thr_height_perc_valley 0.34 --thr_prom_perc_valley 0.30 \
  --bp_merge 1 --bp_edge 500000 --baseline_avg_methyl 0.4 --edge_height_heuristic min \
  > "$OUT/cdrs.bed" 2> "$OUT/cdr_finder.log"
NC=$(wc -l < "$OUT/cdrs.bed"); NCHR=$(cut -f1 "$OUT/cdrs.bed" | sort -u | wc -l)
NWA=$(wc -l < "$OUT/binned_aa.bed")
MM=$(awk '{s+=$4;n++} END{if(n)printf "%.2f",s/n; else print "NA"}' "$OUT/binned_freq.bed")
printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$SAMPLE" "$HAP" "$NC" "$NCHR" "$MM" "$NWA" > "$OUT/summary.tsv"
echo "[OK] $TAG  CDRs=$NC contigs=$NCHR mean_methyl=$MM windows_in_array=$NWA"
