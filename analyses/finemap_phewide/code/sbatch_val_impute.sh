#!/bin/bash
#SBATCH -J ssv_valimp -c 4 --mem=16G -t 0:45:00
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"
set -e
NEW=${CENSV_ROOT}/analyses/finemap_phewide
PIXI=${CENSV_WORKSPACE}/_shared/tools/bayesian_finemap/v1_20260515/pixi.toml
BASE=${CENSV_ROOT}/analyses/finemap_multi_sv
BUN=$BASE/output/bundles/chr4/chr4_super1_47409879_53011266_bundle.npz
cd $NEW
rm -f output/val_impute/manifest.tsv output/val_impute/diag_val.tsv
~/.pixi/bin/pixi run --manifest-path $PIXI python code/build_cells.py --bundle $BUN --chrom 4 --start 47409879 --end 53011266 --arm impute --traits manifests/traits_val.tsv --smoke-cols 89,123 --work-dir output/val_impute/jw --manifest-out output/val_impute/manifest.tsv 2>&1 | grep -vE 'WARN|deprecated|project|platforms|replace|╭|│|·|╰|─'
echo "=== manifest ($(wc -l < output/val_impute/manifest.tsv) rows) ==="; cut -f2-7 output/val_impute/manifest.tsv
~/.pixi/bin/pixi run --manifest-path $PIXI Rscript code/run_susie_batch.R output/val_impute/manifest.tsv output/val_impute/diag_val.tsv "" 2>&1 | grep -vE 'WARN|deprecated|project|platforms|replace|╭|│|·|╰|─'
echo "=== DIAG header + anchor SV rows ==="; head -1 output/val_impute/diag_val.tsv; grep -E 'chr4_sv349|chr4_sv572' output/val_impute/diag_val.tsv
echo "=== VAL_DONE ==="
