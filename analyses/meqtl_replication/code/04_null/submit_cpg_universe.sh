#!/bin/bash
: "${CENSV_COLLAB_ROOT:?set CENSV_COLLAB_ROOT (see README)}"
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -uo pipefail
U=${CENSV_ROOT}/analyses/meqtl_replication
R=${CENSV_COLLAB_ROOT}/exploration/meQTL/reports/chm13_observed_sv_meqtl_mac8_files
OUTROOT=$U/output/null_design/cpg_universe
mkdir -p "$OUTROOT" "$U/logs/b2_cpg_universe" "$U/code/04_null"

awk -F'\t' 'NR>1 && $15!="" && $15!="NA" {print $15}' \
    "$R/all_discovery_conditioning_consensus.tsv" | sort -u > "$OUTROOT/_tasklist.txt"
N=$(wc -l < "$OUTROOT/_tasklist.txt")
echo "tasklist: $N chunks"
head -3 "$OUTROOT/_tasklist.txt"

MISSING=0
while read -r rel; do
  [ -f "${CENSV_COLLAB_ROOT}/exploration/meQTL/$rel" ] || {
    echo "MISSING $rel"; MISSING=$((MISSING+1)); }
done < "$OUTROOT/_tasklist.txt"
echo "missing chunks: $MISSING"
[ "$MISSING" -gt 0 ] && { echo "ABORT: source chunks missing"; exit 1; }

TOTAL_GB=$(while read -r rel; do
  stat -c%s "${CENSV_COLLAB_ROOT}/exploration/meQTL/$rel"
done < "$OUTROOT/_tasklist.txt" | awk '{s+=$1} END{printf "%.1f", s/1073741824}')
echo "total input: ${TOTAL_GB} GB"

KEEP=$(sinfo -h -N -o "%n %t %C" \
       | awk '{split($3,a,"/"); if ($2=="mix" || $2=="idle") print $1, a[2]}' \
       | sort -k2 -rn | head -6 | awk '{print $1}' | sort -u)
ALL=$(sinfo -h -N -o "%n" | sort -u)
EXCL=$(comm -23 <(echo "$ALL") <(echo "$KEEP") | paste -sd, -)
echo "keeping nodes: $(echo $KEEP | tr '\n' ' ')"
echo "excluding    : $EXCL"

JID=$(sbatch --parsable --array=1-${N}%8 --exclude="$EXCL" "$U/code/04_null/build_cpg_universe.sh")
echo "SUBMITTED $JID  (array 1-${N}%8)"
echo "$JID" > "$OUTROOT/_jobid.txt"
