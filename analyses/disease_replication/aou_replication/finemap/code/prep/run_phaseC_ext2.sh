#!/usr/bin/env bash
set -u
FM=/home/dataproc/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication/aou_replication/finemap
LOGS=$FM/logs_ext2
export PLINK2=/home/dataproc/bin/plink2
cd "$FM/code_ext2" || exit 1

wave () {   # $1 label, $2 workers, $3 log, $4... extra args
  local label=$1 w=$2 log=$3; shift 3
  echo "=== $label (workers=$w) start $(date -u +%H:%M:%SZ) ==="
  FINEMAP_GROUP_WORKERS=$w python3 run_finemap_group_batch.py --strata all-anc "$@" > "$LOGS/$log" 2>&1
  local rc=$?
  echo "=== $label rc=$rc $(date -u +%H:%M:%SZ) ==="
  if [ $rc -ne 0 ] || ! grep -q "GROUP BATCH DONE" "$LOGS/$log"; then
    echo "GATE FAIL [$label] -- stopping the chain (partial tree is never built upon)"; return 1
  fi
  local bad; bad=$(grep -cE "fail=[1-9]" "$LOGS/$log" || true)
  [ "${bad:-0}" -ne 0 ] && { echo "GATE FAIL [$label]: $bad group(s) failed"; return 1; }
  echo "GATE PASS [$label]"
  return 0
}

wave "wave1 chr19"          1 phaseC_w1_chr19.log  --chroms chr19        || exit 1
wave "wave2 chr1"           1 phaseC_w2_chr1.log   --chroms chr1         || exit 1
wave "wave3 remaining 17"   3 phaseC_w3_rest.log   --exclude-chroms chr19,chr1 || exit 1
echo "PHASE_C_ALL_WAVES_DONE $(date -u +%H:%M:%SZ)"
