#!/bin/bash
#SBATCH --job-name=cpguniv
#SBATCH --output=${CENSV_ROOT}/analyses/meqtl_replication/logs/b2_cpg_universe/array/%A_%a.log
#SBATCH --time=1:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
: "${CENSV_COLLAB_ROOT:?set CENSV_COLLAB_ROOT (see README)}"
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -uo pipefail

export OMP_NUM_THREADS=1
U=${CENSV_ROOT}/analyses/meqtl_replication
MXROOT=${CENSV_COLLAB_ROOT}/exploration/meQTL
TASKLIST=$U/output/null_design/cpg_universe/_tasklist.txt
OUTROOT=$U/output/null_design/cpg_universe

REL=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TASKLIST")
[ -z "$REL" ] && { echo "FAIL empty task line ${SLURM_ARRAY_TASK_ID}"; exit 1; }

SRC=$MXROOT/$REL
BASE=$(basename "$REL" .methylation.tsv.gz)
CHROM=${BASE%%_*}
DEST_DIR=$OUTROOT/$CHROM
DEST=$DEST_DIR/${BASE}.cpg.tsv.gz
mkdir -p "$DEST_DIR"

echo "TASK ${SLURM_ARRAY_TASK_ID}  chunk=$BASE  chrom=$CHROM"
echo "SRC  $SRC"
[ -f "$SRC" ] || { echo "FAIL source missing"; exit 1; }

if [ -f "$DEST" ]; then
  echo "SKIP already present: $DEST"
  echo "PASS"
  exit 0
fi

TMP=$DEST.tmp.$SLURM_JOB_ID.$SLURM_ARRAY_TASK_ID

zcat "$SRC" | awk -F'\t' '
  NR==1 { next }
  {
    n=0; s=0; ss=0;
    for (i=6; i<=NF; i++) {
      v=$i;
      if (v=="" || v=="NA" || v=="nan" || v=="NaN") continue;
      n++; s+=v; ss+=v*v;
    }
    if (n>0) { m=s/n; var=ss/n-m*m; if (var<0) var=0; sd=sqrt(var); }
    else     { m="NA"; sd="NA"; }
    printf "%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\n", $1,$2,$3,$4,$5,n,m,sd;
  }
' | gzip -c > "$TMP"

rc=${PIPESTATUS[0]}
if [ "$rc" -ne 0 ]; then
  echo "FAIL zcat rc=$rc"
  rm -f "$TMP"
  exit 1
fi

NROW=$(zcat "$TMP" | wc -l)
if [ "$NROW" -lt 1 ]; then
  echo "FAIL zero rows emitted"
  rm -f "$TMP"
  exit 1
fi

mv "$TMP" "$DEST"
echo "rows=$NROW  dest=$DEST"
echo "PASS"
