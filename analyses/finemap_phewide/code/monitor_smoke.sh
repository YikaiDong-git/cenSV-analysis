#!/bin/bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
NEW=${CENSV_ROOT}/analyses/finemap_phewide
LOG=$NEW/logs/STATUS_smoke_progress.log
Q=$(squeue -h -u "$USER" -r -n ssv_smoke | wc -l)
NDIAG=$(ls $NEW/output/smoke/*__*__diag.tsv 2>/dev/null | wc -l)
ROWS=$(tail -q -n +2 $NEW/output/smoke/*__*__diag.tsv 2>/dev/null | wc -l)
echo "$(date '+%F %T') queue=$Q diag_files=$NDIAG/14 total_rows=$ROWS" >> $LOG
if [ "$Q" -eq 0 ] && [ "$NDIAG" -ge 14 ]; then touch $NEW/output/smoke/_SMOKE_DONE; fi
