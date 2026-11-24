#!/bin/bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"

set -u

PROJ=${CENSV_ROOT}/analyses/sv_catalog
LOG=$PROJ/logs/qc/production_6519/STATUS_progress.log
HEARTBEAT=$PROJ/logs/qc/production_6519/STATUS_worker_heartbeat
BUNDLE_ROOT=$PROJ/output/ref_compare_proxy_vs_dense/jwbundle_v3_locked
SWEEP_ROOT=$PROJ/output/ref_compare_proxy_vs_dense/sweep_6519_locked_v3/diag

TARGET_BUNDLE=6519
TARGET_PHENO=6519

mkdir -p "$(dirname "$LOG")"
if [ ! -f "$LOG" ]; then
  printf "timestamp\tqueue_total\tqueue_running\tqueue_pending\tnbundle\tpct_bundle\tndiag_FG\tpct_FG\tndiag_NEALE\tpct_NEALE\tnote\n" > "$LOG"
fi

now=$(date '+%F %T')

qtot=$(squeue -h -u "$USER" -r 2>/dev/null | wc -l)
qrun=$(squeue -h -u "$USER" -t R -r 2>/dev/null | wc -l)
qpd=$((qtot - qrun))

nb=$(ls -d $BUNDLE_ROOT/chr*/*__best*/reference_bundle.npz 2>/dev/null | wc -l)
pb=$(awk -v n=$nb -v t=$TARGET_BUNDLE 'BEGIN{ printf "%.1f", (t>0)?100*n/t:0 }')

nd_fg=$(ls $SWEEP_ROOT/I9_HYPTENSESS/chr*/*__best*__diag.tsv 2>/dev/null | wc -l)
pd_fg=$(awk -v n=$nd_fg -v t=$TARGET_PHENO 'BEGIN{ printf "%.1f", (t>0)?100*n/t:0 }')

nd_neale=$(ls $SWEEP_ROOT/30210_irnt/chr*/*__best*__diag.tsv 2>/dev/null | wc -l)
pd_neale=$(awk -v n=$nd_neale -v t=$TARGET_PHENO 'BEGIN{ printf "%.1f", (t>0)?100*n/t:0 }')

note=""
if [ -f "$HEARTBEAT" ]; then
  hb_age=$(( $(date +%s) - $(stat -c %Y "$HEARTBEAT") ))
  if [ "$qtot" -gt 0 ] && [ "$hb_age" -gt 1800 ]; then
    note="ALERT_STALE_${hb_age}s"
    touch "$PROJ/logs/qc/production_6519/ALERT_STALE_$(date +%Y%m%dT%H%M%S)"
  fi
fi

last_line=$(tail -n 1 "$LOG" 2>/dev/null | cut -f2-)
new_line=$(printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$now" "$qtot" "$qrun" "$qpd" "$nb" "$pb" "$nd_fg" "$pd_fg" "$nd_neale" "$pd_neale" "$note")
new_tail=$(echo "$new_line" | cut -f2-)
if [ "$last_line" != "$new_tail" ]; then
  echo "$new_line" >> "$LOG"
fi

echo "$new_line"
