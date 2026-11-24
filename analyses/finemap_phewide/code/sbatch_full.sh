#!/bin/bash
#SBATCH -J ssv_v2 -c 1 --mem=6G -t 24:00:00
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -e
NEW=${CENSV_ROOT}/analyses/finemap_phewide
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd $NEW
MOD=${STAGGER_MOD:-12}
SLOT=$(( (SLURM_ARRAY_TASK_ID - 1) % MOD ))
sleep $(( SLOT * 300 / MOD + RANDOM % 5 + 120 * ${SLURM_RESTART_COUNT:-0} ))   # ~5min wave regardless of throttle; requeue backoff
source $NEW/code/pixi_env_hook.sh                      # activate env ONCE (pre-generated; no per-task pixi solve)

OFFSET=${OFFSET:-0}                                    # rolling-submit supplies OFFSET; array indices stay <=MaxArraySize
LINE=$((OFFSET + SLURM_ARRAY_TASK_ID))                 # -> covers all 3273 manifest lines across rounds
ROW=$(sed -n "$((LINE + 1))p" ${MANIFEST:-manifests/full_tasks_v2.tsv})
IFS=$'\t' read -r TID BUN ZM MODE ARM CS CE REGION MEM TIME <<< "$ROW"
[ -z "$REGION" ] && { echo "[v2] EMPTY ROW at LINE=$LINE (offset=$OFFSET id=$SLURM_ARRAY_TASK_ID) -- skip"; exit 0; }

SHARD=$(printf 's%02d' $((CS / 500)))                  # extra shard level: <=1000 diag files/folder even for chr19 (1089 blocks)
OUT=$NEW/output/full_v2/$REGION/$ARM/$SHARD; mkdir -p $OUT
DIAG=$OUT/${MODE}_blk$(printf '%06d' $CS)_${CE}_diag.tsv
echo "[v2 task $TID] $REGION $MODE $ARM [$CS,$CE) START $(date +%H:%M:%S) on $(hostname)"
Rscript code/fused_run.R $BUN $ZM $MODE $ARM $CS $CE $REGION $DIAG 2>&1 \
    | grep -vE 'WARN|deprecated|project|platforms|replace|╭|│|·|╰|─' || true
echo "[v2 task $TID] DONE rows=$(tail -n +2 $DIAG 2>/dev/null | wc -l) $(date +%H:%M:%S)"
