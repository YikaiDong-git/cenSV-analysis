#!/bin/bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -u
NEW=${CENSV_ROOT}/analyses/finemap_phewide
cd $NEW
TASKS=$(awk -F'\t' 'NR>1 && $4=="hitbundle"' manifests/full_tasks_v2.tsv | wc -l)
ROUND=450                 # array 1..450  (<1001 cap); logs/full_v2/rNN/ holds <=450 files
THROTTLE=48               # concurrent cap. Lowered from 96: with deterministic stagger this keeps slurmd alive.
CAP=1001; SAFETY=80
LOWATER=100               # top up next round only when my v2 queue < this  => ~1 round active => ~96 running
DONE=${START_AT:-0}
ALL=$(sinfo -h -N -o '%n' | sort -u); EXCL=$(echo "$ALL" | grep -vE '^node0[1247]$' | paste -sd,)   # pin node01/02/04/07 (healthiest at launch; avoid sick node03 load80/node05 load59)
RN=$(( DONE / ROUND ))
echo "$(date '+%F %T') ROLLING START done=$DONE tasks=$TASKS excl=$EXCL"
while [ $DONE -lt $TASKS ]; do
  V2Q=$(squeue -h -u "$USER" -r --name=ssv_v2 2>/dev/null | wc -l)
  ALLQ=$(squeue -h -u "$USER" -r 2>/dev/null | wc -l)
  N=$((TASKS-DONE)); [ $N -gt $ROUND ] && N=$ROUND
  if [ "$V2Q" -lt "$LOWATER" ] && [ $((ALLQ + N + SAFETY)) -le $CAP ]; then
    RN=$((RN+1)); LOGD=logs/full_v2/r$(printf '%02d' $RN); mkdir -p "$LOGD"
    JID=$(sbatch --parsable --exclude=$EXCL --export=ALL,OFFSET=$DONE,STAGGER_MOD=$THROTTLE \
          --array=1-${N}%${THROTTLE} -o "$LOGD/o_%a.out" -e "$LOGD/o_%a.out" code/sbatch_full.sh)
    echo "$(date '+%F %T') SUBMIT round=$RN offset=$DONE n=$N jid=$JID (v2q=$V2Q allq=$ALLQ)"
    DONE=$((DONE+N))
  else
    echo "$(date '+%F %T') WAIT done=$DONE v2q=$V2Q allq=$ALLQ"
  fi
  sleep 120
done
touch logs/full_v2/ROLLING_DONE
echo "$(date '+%F %T') ROLLING COMPLETE all $TASKS submitted"
