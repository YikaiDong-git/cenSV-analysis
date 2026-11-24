#!/bin/bash
#SBATCH --job-name=fg_joint
#SBATCH --time=0:30:00
#SBATCH --mem=16G
#SBATCH -c 4
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"
: "${HOME:?set HOME (see README)}"
set -euo pipefail
REPL=${CENSV_ROOT}/analyses/disease_replication/finngen_replication
MS=${CENSV_ROOT}/analyses/finemap_multi_sv
PIXI=${HOME}/.pixi/bin/pixi
MAN=${CENSV_WORKSPACE}/_shared/tools/bayesian_finemap/v1_20260515/pixi.toml
PY=${CENSV_WORKSPACE}/_shared/tools/bayesian_finemap/v1_20260515/.pixi/envs/default/bin/python
PREP=$REPL/scripts/prepare_jwbundle_finngen.py
SUSIE=$MS/code/03_run_susie_rss_multisv/R/run_susie_rss_multisv.R
JW=$REPL/output/joint/jwbundles; OUT=$REPL/output/joint/susie; RX=$REPL/output/region_extracts
mkdir -p $JW $OUT
echo "=== R lib check ==="
$PIXI run --manifest-path $MAN Rscript -e 'for(p in c("reticulate","susieR","digest")) suppressMessages(library(p,character.only=TRUE)); cat("R_LIBS_OK\n")'
run_one(){
  PHENO=$1; BUNDLE=$2
  echo "=== JOINT $PHENO ==="
  $PY $PREP $BUNDLE $RX/$PHENO.region.tsv $PHENO $JW
  STEM=$(basename $BUNDLE | sed 's/_bundle.npz//')
  $PIXI run --manifest-path $MAN Rscript $SUSIE $JW/${STEM}__${PHENO}__jw.npz $OUT $PHENO
}
B5=$MS/output/bundles/chr5/chr5_super1_43983397_52838788_bundle.npz
B8=$MS/output/bundles/chr8/chr8_super1_41863936_49169911_bundle.npz
run_one H8_VERTIGO      $B5
run_one H8_HL_SEN_NAS   $B8
run_one H8_HL_MIX_NAS   $B8
echo "ALL_JOINT_DONE"
