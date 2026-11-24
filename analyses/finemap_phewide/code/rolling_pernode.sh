#!/bin/bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -u
NEW=${CENSV_ROOT}/analyses/finemap_phewide
cd $NEW
TASKS=$(awk -F'\t' 'NR>1 && $4=="hitbundle"' manifests/full_tasks_v2.tsv | wc -l)
PERNODE=${PERNODE:-12}          # <=12 tasks/node (proven safe); the slurmd-startup ceiling, not a core limit
NNODES=${NNODES:-6}             # weekend: up to 6 nodes (396 concurrent cap); daytime use 4
ROUND=${ROUND:-144}             # chunk/node/round (<1001 array cap; 6*144=864 total queue < cap)
CAP=1001; SAFETY=80
DONE=${START_AT:-0}

pick_nodes() {
  sinfo -h -N -o '%n %t %C %O' | awk -v p=$PERNODE '
    {split($3,a,"/"); idle=a[2]; load=$4; if (idle >= p+2) print load, $1}' | sort -n | head -$NNODES | awk '{print $2}'
}
NODES=( $(pick_nodes) )
[ ${#NODES[@]} -eq 0 ] && { echo "$(date '+%F %T') NO HEALTHY NODES -- abort"; exit 1; }

declare -A NJID
while read -r nd jf; do
  [ -n "$nd" ] && [ -n "$jf" ] && [ "$nd" != "None" ] && NJID[$nd]=$jf
done < <(squeue -h -u "$USER" --name=ssv_v2 -t R -o '%N %F' 2>/dev/null | sort -u)
echo "$(date '+%F %T') PERNODE START done=$DONE tasks=$TASKS pernode=$PERNODE nodes=${NODES[*]} seeded=[${!NJID[*]}]"

RN=$(( DONE / ROUND ))
while [ $DONE -lt $TASKS ]; do
  ALLQ=$(squeue -h -u "$USER" -r 2>/dev/null | wc -l)
  for node in "${NODES[@]}"; do
    [ $DONE -ge $TASKS ] && break
    jid=${NJID[$node]:-}
    nq=0; [ -n "$jid" ] && nq=$(squeue -h -j "$jid" -r 2>/dev/null | wc -l)
    if [ -z "$jid" ] || [ "$nq" -eq 0 ]; then          # node free (no round, or round fully drained)
      N=$((TASKS-DONE)); [ $N -gt $ROUND ] && N=$ROUND
      if [ $((ALLQ + N + SAFETY)) -le $CAP ]; then
        RN=$((RN+1)); LOGD=logs/full_v2/p$(printf '%03d' $RN); mkdir -p "$LOGD"
        NEWJID=$(sbatch --parsable --nodelist=$node --nodes=1 --export=ALL,OFFSET=$DONE,STAGGER_MOD=$PERNODE \
                 --array=1-${N}%${PERNODE} -o "$LOGD/o_%a.out" -e "$LOGD/o_%a.out" code/sbatch_full.sh)
        echo "$(date '+%F %T') SUBMIT node=$node round=$RN offset=$DONE n=$N jid=$NEWJID (allq=$ALLQ)"
        NJID[$node]=$NEWJID; DONE=$((DONE+N)); ALLQ=$((ALLQ+N))
      fi
    fi
  done
  sleep 90
done
touch logs/full_v2/ROLLING_DONE
echo "$(date '+%F %T') PERNODE COMPLETE all $TASKS submitted"
