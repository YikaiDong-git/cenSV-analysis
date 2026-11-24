#!/bin/bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"
NEW=${CENSV_ROOT}/analyses/finemap_phewide
PIXI=${CENSV_WORKSPACE}/_shared/tools/bayesian_finemap/v1_20260515/pixi.toml
PIXIRUN="$HOME/.pixi/bin/pixi run --manifest-path $PIXI"
JID=3826876
LOG=$NEW/logs/keyint_progress.log
TS=$(date "+%Y-%m-%d %H:%M:%S")
Q=$(squeue -h -j $JID -r 2>/dev/null | wc -l)
RUN=$(squeue -h -j $JID -r -t R 2>/dev/null | wc -l)
D=$(ls $NEW/output/keyinterest/*__*__diag.tsv 2>/dev/null | wc -l)
ROWS=$(cat $NEW/output/keyinterest/*__*__diag.tsv 2>/dev/null | grep -vc '^region')
echo "$TS queue=$Q running=$RUN diag_files=$D/24 rows=$ROWS" >> $LOG
if [ "$Q" -eq 0 ] && [ "$D" -ge 24 ] && [ ! -f $NEW/logs/KEYINT_DONE ]; then
  echo "$TS KEYINT_SMOKE_DONE diag=$D rows=$ROWS -- running aggregate" >> $LOG
  cd $NEW && $PIXIRUN python code/aggregate_keyinterest.py "output/keyinterest/*__*__diag.tsv" \
      manifests/keyinterest_svs.tsv manifests/traits_all464.tsv \
      output/keyinterest_wide.tsv output/keyinterest_digest.txt output/keyinterest_candidates.tsv \
      >> $LOG 2>&1
  echo "$TS AGGREGATE_DONE" >> $LOG
  touch $NEW/logs/KEYINT_DONE
fi
