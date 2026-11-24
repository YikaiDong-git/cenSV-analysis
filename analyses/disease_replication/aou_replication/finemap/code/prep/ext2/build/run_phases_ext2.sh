#!/usr/bin/env bash
set -u
FM=/home/dataproc/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication/aou_replication/finemap
LOGS=$FM/logs_ext2
export PLINK2=/home/dataproc/bin/plink2
cd "$FM/code_ext2" || exit 1

wait_for_phase () {   # $1 = log file, $2 = human label, $3 = optional PID to wait on
  local log=$1 label=$2 pid=${3:-}
  if [ -n "$pid" ]; then
    while kill -0 "$pid" 2>/dev/null; do sleep 60; done
  fi
  if ! grep -q "GROUP BATCH DONE" "$log" 2>/dev/null; then
    echo "GATE FAIL [$label]: no 'GROUP BATCH DONE' in $log — stopping the chain."; return 1
  fi
  local bad
  bad=$(grep -cE "fail=[1-9]" "$log" 2>/dev/null || true)
  if [ "${bad:-0}" -ne 0 ]; then
    echo "GATE FAIL [$label]: $bad group(s) reported failures — stopping the chain."; return 1
  fi
  echo "GATE PASS [$label]: $(grep -c ' susie=' "$log") groups, 0 failures."
  return 0
}

run_phase () {        # $1 = strata, $2 = workers, $3 = log basename, $4 = label
  echo "=== launching $4 (strata=$1, workers=$2) at $(date -u +%H:%M:%SZ) ==="
  FINEMAP_GROUP_WORKERS=$2 python3 run_finemap_group_batch.py --strata "$1" > "$LOGS/$3" 2>&1
  local rc=$?
  echo "$4 exited rc=$rc at $(date -u +%H:%M:%SZ)"
  return $rc
}

echo "driver: waiting on Phase A (already running)"
wait_for_phase "$LOGS/phaseA_amr_afr.log" "Phase A (AMR+AFR)" "${PHASE_A_PID:-}" || exit 1

run_phase EUR 3 phaseB_eur.log "Phase B (EUR)" \
  && wait_for_phase "$LOGS/phaseB_eur.log" "Phase B (EUR)" || exit 1

run_phase all-anc 3 phaseC_allanc.log "Phase C (all-anc)" \
  && wait_for_phase "$LOGS/phaseC_allanc.log" "Phase C (all-anc)" || exit 1

echo "ALL_PHASES_DONE $(date -u +%H:%M:%SZ)"
echo "next: Gate 0 completeness check, then code_ext2/collate_finemap.py"
