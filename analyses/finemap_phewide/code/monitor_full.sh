#!/bin/bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
NEW=${CENSV_ROOT}/analyses/finemap_phewide
JID=3828631
LOG=$NEW/logs/full_progress.log
TS=$(date "+%Y-%m-%d %H:%M:%S")
Q=$(squeue -h -j $JID -r 2>/dev/null | wc -l)
RUN=$(squeue -h -j $JID -r -t R 2>/dev/null | wc -l)
D=$(ls $NEW/output/full/*/*/*_diag.tsv 2>/dev/null | wc -l)
echo "$TS queue=$Q running=$RUN diag_files=$D/776" >> $LOG
if [ "$Q" -eq 0 ] && [ "$D" -ge 700 ] && [ ! -f $NEW/logs/FULL_DONE ]; then
  cd $NEW
  AGG=$(sbatch --parsable code/sbatch_aggregate_full.sh 2>/dev/null)
  echo "$TS FULL_ARRAY_DONE diag=$D -- submitted aggregate sbatch=$AGG" >> $LOG
  touch $NEW/logs/FULL_DONE
fi
