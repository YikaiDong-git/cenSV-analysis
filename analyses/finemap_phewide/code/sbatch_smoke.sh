#!/bin/bash
#SBATCH -J ssv_smoke -c 4 --mem=20G -t 4:00:00
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_SCRATCH:?set CENSV_SCRATCH (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"
set -e
NEW=${CENSV_ROOT}/analyses/finemap_phewide
PIXI=${CENSV_WORKSPACE}/_shared/tools/bayesian_finemap/v1_20260515/pixi.toml
PY="$HOME/.pixi/bin/pixi run --manifest-path $PIXI python"
RS="$HOME/.pixi/bin/pixi run --manifest-path $PIXI Rscript"
cd $NEW
ROW=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" manifests/smoke_tasks.tsv)
IFS=$'\t' read -r TID RB CHROM START END COLS ARM <<< "$ROW"
REGION=$(basename "$RB" | sed 's/_bundle.npz//')
WORK=${CENSV_SCRATCH}/ssv_smoke/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p $WORK/jw output/smoke
MAN=$WORK/manifest.tsv; : > $MAN
SMOKE_LOCI=$(awk -F'\t' -v rb="$RB" '$4==rb{print $1}' manifests/smoke_svs.tsv | paste -sd,)
echo "[task $TID] region=$REGION arm=$ARM cols=$COLS smoke=$SMOKE_LOCI START=$(date +%H:%M:%S)"
$PY code/build_cells.py --bundle "$RB" --chrom $CHROM --start $START --end $END --arm $ARM \
    --traits manifests/traits_all464.tsv --smoke-cols "$COLS" --work-dir $WORK/jw --manifest-out $MAN 2>&1 \
    | grep -vE 'WARN|deprecated|project|platforms|replace|╭|│|·|╰|─' || true
echo "[task $TID] built $(wc -l < $MAN) jw bundles; fitting $(date +%H:%M:%S)"
$RS code/run_susie_batch.R $MAN output/smoke/${REGION}__${ARM}__diag.tsv "$SMOKE_LOCI" 2>&1 \
    | grep -vE 'WARN|deprecated|project|platforms|replace|╭|│|·|╰|─' || true
rm -rf $WORK
echo "[task $TID] DONE region=$REGION arm=$ARM $(date +%H:%M:%S) rows=$(tail -n +2 output/smoke/${REGION}__${ARM}__diag.tsv 2>/dev/null | wc -l)"
