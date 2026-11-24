#!/bin/bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
#SBATCH --job-name=gao_pheno
#SBATCH --mem=24G
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --output=${CENSV_ROOT}/analyses/cdr_definition/logs/gao_pheno_%j.log
set -euo pipefail
export OMP_NUM_THREADS=1
PROJ=${CENSV_ROOT}
U=$PROJ/analyses/cdr_definition
SHIPPED=$PROJ/analyses/cdr_phenotypes/output/phenotypes/per_hap_chrom_phenotypes_frames.tsv
cd "$PROJ"

echo "[$(date)] === 1. reproduction gate: same code, original dips ==="
~/.pixi/bin/pixi run python "$U/code/pipeline/03_recompute_phenotypes.py" \
    --out-name per_hap_chrom_phenotypes_reproduced.tsv

if diff -q "$SHIPPED" "$U/output/phenotypes/per_hap_chrom_phenotypes_reproduced.tsv" > /dev/null; then
  echo "REPRO_GATE_PASS  the fork reproduces the shipped table byte-for-byte"
else
  echo "REPRO_GATE_FAIL  the fork changed something other than the call set"
  diff "$SHIPPED" "$U/output/phenotypes/per_hap_chrom_phenotypes_reproduced.tsv" | head -20
  exit 1
fi

echo "[$(date)] === 2. Gao CDR unit ==="
~/.pixi/bin/pixi run python "$U/code/pipeline/03_recompute_phenotypes.py" \
    --dips-root "$U/output/callsets/gao_cdrs" \
    --out-name per_hap_chrom_phenotypes_gao.tsv

echo "[$(date)] PHENO_DONE"
