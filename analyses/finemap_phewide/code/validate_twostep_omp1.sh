#!/bin/bash
#SBATCH -J ssv_ts1 -c 1 --mem=6G -t 0:30:00
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_SCRATCH:?set CENSV_SCRATCH (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"
set -e
NEW=${CENSV_ROOT}/analyses/finemap_phewide
BUN=${CENSV_ROOT}/analyses/finemap_multi_sv/output/bundles/chr16/chr16_super1_31371622_48974015_bundle.npz
PIXI=${CENSV_WORKSPACE}/_shared/tools/bayesian_finemap/v1_20260515/pixi.toml
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd $NEW
WORK=${CENSV_SCRATCH}/val_ts1/$SLURM_JOB_ID; mkdir -p $WORK
for arm in impute matched; do
  MAN=$WORK/man_$arm.tsv; : > $MAN
  $HOME/.pixi/bin/pixi run --manifest-path $PIXI python code/build_cells.py --bundle "$BUN" \
     --chrom 16 --start 31371622 --end 48974015 --arm $arm --traits manifests/traits_all464.tsv \
     --smoke-cols "74,77" --modes hitbundle --work-dir $WORK/jw_$arm --manifest-out $MAN 2>&1 \
     | grep -vE 'WARN|deprecated|project|platforms|replace|╭|│|·|╰|─' || true
  $HOME/.pixi/bin/pixi run --manifest-path $PIXI Rscript code/run_susie_batch.R \
     $MAN output/twostep_omp1_$arm.tsv "" 2>&1 \
     | grep -vE 'WARN|deprecated|project|platforms|replace|╭|│|·|╰|─' || true
done
rm -rf $WORK
echo "TWOSTEP_OMP1 DONE impute=$(tail -n +2 output/twostep_omp1_impute.tsv 2>/dev/null|wc -l) matched=$(tail -n +2 output/twostep_omp1_matched.tsv 2>/dev/null|wc -l)"
