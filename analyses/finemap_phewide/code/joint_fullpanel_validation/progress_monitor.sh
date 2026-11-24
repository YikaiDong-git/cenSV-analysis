#!/bin/bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_SCRATCH:?set CENSV_SCRATCH (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"
FAN=${CENSV_SCRATCH}/ssv_joint
CJ=${CENSV_ROOT}/analyses/finemap_phewide/code/joint_fullpanel_validation
PIXI=~/.pixi/bin/pixi
PM=${CENSV_WORKSPACE}/_shared/tools/bayesian_finemap/v1_20260515/pixi.toml
LOG=$FAN/STATUS_progress.log
TS=$(date "+%Y-%m-%d %H:%M:%S")

QALL=$(squeue -h -u ydon268 -r | wc -l)
A_DIAG=$(ls $FAN/fanout_out/0p3/*_fanout_diag.tsv 2>/dev/null | wc -l)
A_CHUNKS=$(wc -l < $FAN/chunks_0p3.txt 2>/dev/null || echo 0)
A_Q=$(squeue -h -u ydon268 -r -n jfan 2>/dev/null | wc -l)
B_DIAG=$(ls $FAN/fanout_out/0p1/*_fanout_diag.tsv 2>/dev/null | wc -l)
B_CHUNKS=$(wc -l < $FAN/chunks_0p1.txt 2>/dev/null || echo 0)
B_Q=$(squeue -h -u ydon268 -r -n jfan0p1 2>/dev/null | wc -l)

echo "[$TS] Qall=$QALL | A_fanout(0p3) diag=$A_DIAG/$A_CHUNKS q=$A_Q | B_fanout(0p1) diag=$B_DIAG/$B_CHUNKS q=$B_Q" >> $LOG

if [ "$A_Q" = "0" ] && [ "$A_DIAG" -ge 1 ] && [ "$A_CHUNKS" -gt 0 ] && [ "$A_DIAG" -ge $((A_CHUNKS*97/100)) ] && [ ! -f $FAN/fanout_out/agg_0p3/DONE ]; then
  $PIXI run --manifest-path $PM python $CJ/aggregate_joint_fanout.py 0p3 > $FAN/fanout_out/agg_0p3/verdict_0p3.txt 2>&1
  touch $FAN/fanout_out/agg_0p3/DONE
  echo "[$TS] AGG_0p3_DONE" >> $LOG
fi
if [ "$A_Q" = "0" ] && [ "$B_Q" = "0" ] && [ "$A_DIAG" -ge $((A_CHUNKS*97/100)) ] && [ -f $FAN/fanout_out/agg_0p3/DONE ] && [ ! -f $FAN/fanout_out/agg_0p1_FULL/DONE ]; then
  mkdir -p $FAN/fanout_out/agg_0p1_FULL
  $PIXI run --manifest-path $PM python $CJ/finalize_0p1.py > $FAN/fanout_out/agg_0p1_FULL/verdict_0p1_FULL.txt 2>&1
  touch $FAN/fanout_out/agg_0p1_FULL/DONE
  echo "[$TS] AGG_0p1_FULL_DONE  ==> ALL_FANOUTS_AGGREGATED" >> $LOG
fi
