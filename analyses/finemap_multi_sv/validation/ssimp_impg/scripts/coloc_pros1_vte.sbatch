#!/bin/bash
#SBATCH --job-name=coloc_pros1
#SBATCH --mem=32G
#SBATCH -c 4
#SBATCH --nodes=1
#SBATCH --exclude=node01,node02,node03,node04,node09,node10,node11,node12,node13,node14
#SBATCH --time=1:00:00
#SBATCH --output=${CENSV_ROOT}/analyses/finemap_multi_sv/validation/ssimp_impg/work/coloc_pros1/coloc_pros1_vte.%j.out
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"
set -euo pipefail
P=${CENSV_WORKSPACE}
BASE=$P/Project1_Centromere/analyses/finemap_multi_sv
V=$BASE/validation/ssimp_impg
BT=$P/_shared/envs/workspace/v1/.pixi/envs/biotools/bin
EQGZ=$P/_shared/datasets/eQTL/eqtlgen_phase1_cis/v_20191211/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz
PIXI=~/.pixi/bin/pixi
RSTATS_MANIFEST=$P/_shared/envs/workspace/v1/pixi.toml
WD=$V/work/coloc_pros1; mkdir -p $WD
if [ ! -s $WD/eqtlgen_PROS1_hg19.tsv ]; then
  echo "[$(date)] extract eQTLGen PROS1"
  zcat $EQGZ | awk -F'\t' 'NR==1 || $9=="PROS1"' > $WD/eqtlgen_PROS1_hg19.tsv
fi
echo "[extract] eqtlgen_PROS1_hg19.tsv rows=$(wc -l < $WD/eqtlgen_PROS1_hg19.tsv)"
echo "[$(date)] prep"
$BT/python $V/scripts/coloc_pros1_vte_prep.py
echo "[$(date)] coloc"
$PIXI run --manifest-path "$RSTATS_MANIFEST" -e rstats Rscript $V/scripts/coloc_pros1_vte.R
echo "[$(date)] DONE_COLOC"
