#!/bin/bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -uo pipefail
P=${CENSV_ROOT}
U=$P/analyses/meqtl_replication
CH=$U/output/headline_svs/chunks
NCHUNK=${1:-16}

mkdir -p "$CH" "$U/logs/headline/array" "$U/output/headline_svs/per_hap"
rm -f "$CH"/chunk_*.txt

TAGS=$U/output/a1_lift/all_tags.txt
N=$(wc -l < "$TAGS")
echo "haplotypes: $N   chunks: $NCHUNK"

i=0
while read -r tag; do
  [ -z "$tag" ] && continue
  idx=$(( i % NCHUNK + 1 ))
  echo "$tag" >> "$CH/chunk_$(printf '%03d' $idx).txt"
  i=$((i+1))
done < "$TAGS"

echo "chunk sizes:"
wc -l "$CH"/chunk_*.txt | head -20
TOTAL=$(cat "$CH"/chunk_*.txt | wc -l)
echo "total tags across chunks: $TOTAL  (must equal $N)"
[ "$TOTAL" -ne "$N" ] && { echo "ABORT chunk split lost tags"; exit 1; }
DUP=$(cat "$CH"/chunk_*.txt | sort | uniq -d | wc -l)
[ "$DUP" -ne 0 ] && { echo "ABORT $DUP duplicated tags across chunks"; exit 1; }

echo "target set:"
wc -l "$U/output/targets/targets_v2_lift.tsv"

CUR=$(squeue -h -u "$USER" -r | wc -l)
echo "my current queue depth: $CUR  (cap 1000, shared across projects)"
[ "$CUR" -gt 900 ] && { echo "ABORT queue too deep"; exit 1; }

KEEP=$(sinfo -h -N -o "%n %t %C" \
       | awk '{split($3,a,"/"); if ($2=="mix" || $2=="idle") print $1, a[2]}' \
       | sort -k2 -rn | head -6 | awk '{print $1}' | sort -u)
ALL=$(sinfo -h -N -o "%n" | sort -u)
EXCL=$(comm -23 <(echo "$ALL") <(echo "$KEEP") | paste -sd, -)
echo "keeping: $(echo $KEEP | tr '\n' ' ')"

JID=$(sbatch --parsable --array=1-${NCHUNK}%8 --exclude="$EXCL" "$U/code/09_headline/lift_headline_array.sh")
echo "SUBMITTED $JID  (array 1-${NCHUNK}%8, ~$((N / NCHUNK)) haplotypes per task)"
echo "$JID" > "$U/output/headline_svs/_run/_lift_jobid.txt"
