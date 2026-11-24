#!/bin/bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
NEW=${CENSV_ROOT}/analyses/finemap_phewide
cd $NEW
LOG=logs/full_v2/progress.log
TS=$(date '+%F %T')
V2Q=$(squeue -h -u "$USER" -r --name=ssv_v2 2>/dev/null | wc -l)
V2R=$(squeue -h -u "$USER" -r --name=ssv_v2 -t R 2>/dev/null | wc -l)
D=$(ls output/full_v2/*/*/*/*_diag.tsv 2>/dev/null | wc -l)
echo "$TS v2queue=$V2Q running=$V2R diag=$D/3273" >> $LOG
if ! tmux has-session -t ssv_v2_submit 2>/dev/null && [ ! -f logs/full_v2/ROLLING_DONE ]; then
  RESUME=$(grep -oE 'offset=[0-9]+ n=[0-9]+' logs/full_v2/rolling_pernode.log 2>/dev/null | sed 's/offset=//;s/n=//' | awk '{print $1+$2}' | sort -n | tail -1)
  RESUME=${RESUME:-0}
  tmux new -d -s ssv_v2_submit "START_AT=$RESUME bash $NEW/code/rolling_pernode.sh >> $NEW/logs/full_v2/rolling_pernode.log 2>&1"
  echo "$TS SELF-HEAL: rolling tmux was dead -> restarted START_AT=$RESUME" >> $LOG
fi
if [ -f logs/full_v2/ROLLING_DONE ] && [ "$V2Q" -eq 0 ] && [ "$D" -ge 3000 ] && [ ! -f logs/full_v2/FULL_V2_DONE ]; then
  AGG=$(sbatch --parsable code/sbatch_aggregate.sh 2>/dev/null)
  echo "$TS V2_COMPLETE diag=$D -- submitted aggregate=$AGG" >> $LOG
  touch logs/full_v2/FULL_V2_DONE
fi
