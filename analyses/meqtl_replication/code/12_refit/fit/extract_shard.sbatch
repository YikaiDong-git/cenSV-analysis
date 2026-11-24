#!/bin/bash
#SBATCH --job-name=cpgext
#SBATCH --output=${CENSV_ROOT}/analyses/meqtl_replication/logs/refit/ext_%A_%a.log
#SBATCH --time=1:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
: "${CENSV_COLLAB_ROOT:?set CENSV_COLLAB_ROOT (see README)}"
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -uo pipefail
P=${CENSV_ROOT}
U=$P/analyses/meqtl_replication
OUT=$U/output/arms/refit_genotype_complete
MAT=${CENSV_COLLAB_ROOT}/exploration/meQTL/results/meqtl_cis/chm13_native_same_chrom/matrices

SHARD=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$OUT/r2_shards.txt")
[ -n "$SHARD" ] || { echo "FAIL no shard for task ${SLURM_ARRAY_TASK_ID}"; exit 1; }
BASE=${SHARD%.methylation.tsv.gz}
CHROM=${BASE%%_*}
IDS="$OUT/wanted/$CHROM/$BASE.ids"
DEST_DIR="$OUT/extracted/$CHROM"
mkdir -p "$DEST_DIR"
DEST="$DEST_DIR/$BASE.tsv"

echo "HOST $(hostname)  SHARD $SHARD  START $(date -u +%Y-%m-%dT%H:%M:%SZ)"
zcat "$MAT/$SHARD" | awk -F'\t' -v OFS='\t' -v idsfile="$IDS" '
BEGIN{ while((getline l < idsfile) > 0) want[l]=1 }
NR==1{ print; next }
($4 in want){ print }
' > "$DEST"
rc=$?

GOT=$(( $(wc -l < "$DEST") - 1 ))
WANT=$(wc -l < "$IDS")
echo "EXTRACTED $GOT / $WANT  rc=$rc"
[ $rc -eq 0 ] && echo PASS || echo "FAIL rc=$rc"
echo "END $(date -u +%Y-%m-%dT%H:%M:%SZ)"
exit $rc
