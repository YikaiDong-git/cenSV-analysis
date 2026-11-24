#!/bin/bash
#SBATCH --job-name=fg_region
#SBATCH --time=0:20:00
#SBATCH --mem=2G
#SBATCH -c 1
#SBATCH --output=${CENSV_ROOT}/analyses/disease_replication/finngen_replication/logs/region_%A_%a.log
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail
PB=${CENSV_ROOT}/centromere_twas/finngen/processed_b38aligned
OUT=${CENSV_ROOT}/analyses/disease_replication/finngen_replication/output/region_extracts
LIST=${CENSV_ROOT}/analyses/disease_replication/finngen_replication/manifests/all_endpoints.txt
mkdir -p $OUT
P=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $LIST)
[ -z "$P" ] && { echo "no endpoint at line $SLURM_ARRAY_TASK_ID"; exit 2; }
TSV=$PB/$P/$P.hg38_b38aligned.tsv
[ -f "$TSV" ] || { echo "missing $TSV"; exit 3; }
TMP=$OUT/$P.region.tsv.tmp.$$
head -1 "$TSV" > $TMP
awk -F'\t' -v RS5=43983397 -v RE5=52838788 -v RS8=41863936 -v RE8=49169911 '
  NR==1{next}
  { n=split($1,a,"_"); if(n<4) next; c=a[1]; p=a[2]+0;
    if((c=="chr5"&&p>=RS5&&p<=RE5)||(c=="chr8"&&p>=RS8&&p<=RE8)) print }
' "$TSV" >> $TMP
mv -f $TMP $OUT/$P.region.tsv
echo "[$P] region rows: $(( $(wc -l < $OUT/$P.region.tsv) - 1 ))"
