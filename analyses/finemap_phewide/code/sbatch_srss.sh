#!/bin/bash
#SBATCH -J ssv_srss -c 2 --mem=12G -t 0:30:00
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_SCRATCH:?set CENSV_SCRATCH (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"
set -e
NEW=${CENSV_ROOT}/analyses/finemap_phewide
PIXI=${CENSV_WORKSPACE}/_shared/tools/bayesian_finemap/v1_20260515/pixi.toml
PY="$HOME/.pixi/bin/pixi run --manifest-path $PIXI python"
RS="$HOME/.pixi/bin/pixi run --manifest-path $PIXI Rscript"
cd $NEW
WORK=${CENSV_SCRATCH}/ssv_srss/${SLURM_JOB_ID}
mkdir -p $WORK output/keyinterest
$PY code/srss_build.py output/keyinterest_candidates.tsv $WORK $WORK/srss_manifest.tsv 2>&1 \
    | grep -vE 'WARN|deprecated|project|platforms|replace|╭|│|·|╰|─' || true
$RS code/srss_run.R $WORK/srss_manifest.tsv output/keyinterest_srss.tsv 2>&1 \
    | grep -vE 'WARN|deprecated|project|platforms|replace|╭|│|·|╰|─' || true
rm -rf $WORK
echo "=== keyinterest_srss.tsv ==="; cat output/keyinterest_srss.tsv
