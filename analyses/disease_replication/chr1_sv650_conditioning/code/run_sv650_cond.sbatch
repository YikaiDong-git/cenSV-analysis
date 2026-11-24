#!/bin/bash
#SBATCH --job-name=sv650cond
#SBATCH --mem=16G
#SBATCH -c 4
#SBATCH --nodes=1
#SBATCH --time=0:30:00
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"
set -uo pipefail
P=${CENSV_WORKSPACE}
BASE=$P/Project1_Centromere/analyses/finemap_multi_sv
V=$BASE/validation/ssimp_impg
BT=$P/_shared/envs/workspace/v1/.pixi/envs/biotools/bin
PIXI=$HOME/.pixi/bin/pixi
PIXI_MANIFEST=$P/_shared/tools/bayesian_finemap/v1_20260515/pixi.toml
GD=$P/Project1_Centromere/analyses/sv_catalog/output/ref_compare_proxy_vs_dense/gwas_bgzip_tabix
export V2_CANONICAL_QC_R_DIR=$P/Project1_Centromere/analyses/sv_catalog/code/q8_validation/qc/R
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4

REGION=chr1_super1_118296444_148907004
CHROM=1; START=118296444; END=148907004
BUNDLE=$BASE/output/bundles/chr1/${REGION}_bundle.npz
SV=chr1_sv650_INS_ori_best8_fksize500k
STUDY=finngen
WORK=$P/Project1_Centromere/analyses/disease_replication/chr1_sv650_conditioning
mkdir -p $WORK/multi $WORK/single

PHENOS="C3_BREAST_EXALLC C3_BREAST_ERNEG_EXALLC C3_BREAST_ERPLUS_EXALLC"

for PHENO in $PHENOS; do
  BGZ=$GD/$STUDY/$PHENO.tsv.bgz
  [ -f "$BGZ" ] || { echo "[skip] no gwas $BGZ"; continue; }
  echo "===== $PHENO ====="
  $BT/python $V/scripts/build_matched_jwbundle.py --bundle $BUNDLE --gwas-bgz $BGZ \
      --chrom $CHROM --start $START --end $END --pheno $PHENO --region-id $REGION --out-dir $WORK/multi
  JWm=$WORK/multi/${REGION}__${PHENO}__jw.npz
  [ -f "$JWm" ] && $PIXI run --manifest-path $PIXI_MANIFEST Rscript \
      $BASE/code/03_run_susie_rss_multisv/R/run_susie_rss_multisv.R $JWm $WORK/multi $PHENO $REGION
  $BT/python $WORK/code/build_single_sv_jwbundle.py --bundle $BUNDLE --gwas-bgz $BGZ \
      --chrom $CHROM --start $START --end $END --pheno $PHENO --region-id $REGION \
      --out-dir $WORK/single --only-sv-locus $SV
  JWs=$WORK/single/${REGION}__${PHENO}__jw.npz
  [ -f "$JWs" ] && $PIXI run --manifest-path $PIXI_MANIFEST Rscript \
      $BASE/code/03_run_susie_rss_multisv/R/run_susie_rss_multisv.R $JWs $WORK/single $PHENO $REGION
done

echo "=== SUMMARY chr1_sv650: multi-SV vs single-SV (sv_pip) ==="
printf "%-28s %-14s %-8s %-26s %-8s\n" PHENO multiSV_pip G_multi singleSV_pip_topSNP p_snps
for PHENO in $PHENOS; do
  m=$WORK/multi/${REGION}__${PHENO}__diag.tsv
  s=$WORK/single/${REGION}__${PHENO}__diag.tsv
  sp=$WORK/single/${REGION}__${PHENO}__pip_snp.tsv
  mpip=$( [ -f "$m" ] && awk -F'\t' 'NR>1 && $3 ~ /_sv650_/ {print $9}' "$m" || echo NA )
  gm=$( [ -f "$m" ] && awk -F'\t' 'NR==2{print $11}' "$m" || echo NA )
  spip=$( [ -f "$s" ] && awk -F'\t' 'NR>1 && $3 ~ /_sv650_/ {print $9}' "$s" || echo NA )
  psn=$( [ -f "$s" ] && awk -F'\t' 'NR==2{print $10}' "$s" || echo NA )
  topsnp=$( [ -f "$sp" ] && tail -n +2 "$sp" | sort -t$'\t' -k4 -gr | head -1 | awk -F'\t' '{printf "%.3f@%s", $4, $3}' || echo NA )
  printf "%-28s %-14s %-8s SVpip=%-6s top=%-12s %-8s\n" "$PHENO" "$mpip" "$gm" "$spip" "$topsnp" "$psn"
done
echo DONE_SV650COND
