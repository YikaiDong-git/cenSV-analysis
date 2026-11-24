#!/bin/bash
#SBATCH --job-name=p6519_pBresc
#SBATCH --output=${CENSV_ROOT}/analyses/sv_catalog/logs/qc/production_6519/phaseB_rescue_%A_%a.log
#SBATCH --error=${CENSV_ROOT}/analyses/sv_catalog/logs/qc/production_6519/phaseB_rescue_%A_%a.err
#SBATCH --time=6:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"
#SBATCH --exclude=node03,node04,node06,node07,node08,node09,node10,node11,node12,node13

set -euo pipefail
sleep $((RANDOM % 4))

PROJ=${CENSV_ROOT}/analyses/sv_catalog
CODE_DIR_REF=$PROJ/code/ref_compare_proxy_vs_dense
V2=$PROJ/code/q8_validation/qc
PIXI=${CENSV_WORKSPACE}/_shared/tools/bayesian_finemap/v1_20260515
PY=$PIXI/.pixi/envs/default/bin/python
RSCRIPT=$PIXI/.pixi/envs/default/bin/Rscript

: "${PHENO:?PHENO env var required}"
: "${GWAS_PATH:?GWAS_PATH env var required}"
: "${STUDY:?STUDY env var required (neale|finngen)}"

MANIFEST=${MANIFEST_OVERRIDE:-$PROJ/output/ref_compare_proxy_vs_dense/sweep_6519_locked_v3/manifest/MISSING_RESCUE_${PHENO}.tsv}
BUNDLE_ROOT=$PROJ/output/ref_compare_proxy_vs_dense/jwbundle_v3_locked
SWEEP_ROOT=$PROJ/output/ref_compare_proxy_vs_dense/sweep_6519_locked_v3/rescue

if [ ! -f "$GWAS_PATH" ]; then
  echo "[FATAL] GWAS_PATH missing: $GWAS_PATH"
  exit 4
fi
if [ ! -f "$MANIFEST" ]; then
  echo "[FATAL] MANIFEST missing: $MANIFEST"
  exit 5
fi

K=${K_OVERRIDE:-2}
BASE=$(( (SLURM_ARRAY_TASK_ID - 1) * K ))

P_MAX_FOR_RAUG=50000

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export V2_CANONICAL_QC_R_DIR=$V2/R

export PATH=${CENSV_ROOT}/.pixi/envs/default/bin:$PATH

TMP_BASE=/tmp/${USER}_p6519BR_${PHENO}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p $TMP_BASE
trap 'rm -rf "$TMP_BASE"' EXIT

cd $PIXI

n_done=0
n_skip=0
n_gap=0
n_built=0
n_too_wide=0
HEARTBEAT_FILE=$PROJ/logs/qc/production_6519/STATUS_rescue_heartbeat
mkdir -p $(dirname "$HEARTBEAT_FILE")

echo "=========================================="
echo "[$(date '+%F %T')] phaseB_rescue jid=${SLURM_ARRAY_JOB_ID} task=${SLURM_ARRAY_TASK_ID} node=$(hostname)"
echo "[rescue] PHENO=$PHENO  STUDY=$STUDY  GWAS=$GWAS_PATH"
echo "[rescue] K=$K  BASE=$BASE  range=$((BASE+1)).. $((BASE+K))  cap=$P_MAX_FOR_RAUG"
echo "[rescue] MANIFEST=$MANIFEST"
echo "=========================================="

for i in $(seq 0 $((K - 1))); do
  PAIR_ID_LOCAL=$((BASE + i + 1))
  LINE=$(awk -F'\t' -v r=$((PAIR_ID_LOCAL + 1)) 'NR==r' "$MANIFEST")
  if [ -z "$LINE" ]; then
    echo "[rescue] task-local pair_id=$PAIR_ID_LOCAL past EOF; stopping inner loop"
    break
  fi

  PAIR_ID=$(echo "$LINE"    | cut -f1)
  SV_LOCUS=$(echo "$LINE"   | cut -f2)
  SV_CHROM=$(echo "$LINE"   | cut -f4)
  ARM=$(echo "$LINE"        | cut -f7)
  BUNDLE_REL=$(echo "$LINE" | cut -f10)
  RESCUE_GROUP=$(echo "$LINE" | cut -f14)

  REF_BUNDLE=$BUNDLE_ROOT/$BUNDLE_REL
  OUT_DIR=$SWEEP_ROOT/$PHENO/$SV_CHROM
  OUT_DIAG=$OUT_DIR/${SV_LOCUS}__${ARM}__rescue_diag.tsv

  if [ -s "$OUT_DIAG" ]; then
    n_skip=$((n_skip + 1))
    continue
  fi

  if [ ! -s "$REF_BUNDLE/reference_bundle.npz" ]; then
    echo "[phase_a_missing] master_pair_id=$PAIR_ID  $SV_LOCUS  $ARM  (no Phase A bundle)"
    mkdir -p "$OUT_DIR"
    TMP_STUB=$OUT_DIAG.tmp.$$
    printf "phenotype\tsv_locus\tarm\tchrom\tstatus\trescue_group\tsv_pip\ttop_snp_pip\tsum_pip\n" >  "$TMP_STUB"
    printf "%s\t%s\t%s\t%s\tPHASE_A_MISSING\t%s\t\t\t\n" "$PHENO" "$SV_LOCUS" "$ARM" "$SV_CHROM" "$RESCUE_GROUP" >> "$TMP_STUB"
    mv "$TMP_STUB" "$OUT_DIAG"
    touch "$HEARTBEAT_FILE"
    continue
  fi

  mkdir -p "$OUT_DIR"

  SNP_META=$REF_BUNDLE/snp_meta.tsv.gz
  if [ ! -s "$SNP_META" ]; then
    echo "[FATAL] master_pair_id=$PAIR_ID snp_meta.tsv.gz missing at $SNP_META"
    exit 9
  fi
  N_PANEL_SNP=$(zcat "$SNP_META" | wc -l)
  N_PANEL_SNP=$((N_PANEL_SNP - 1))
  if [ "$N_PANEL_SNP" -gt "$P_MAX_FOR_RAUG" ]; then
    echo "[wide_64G] master_pair_id=$PAIR_ID  $SV_LOCUS  $ARM  p=$N_PANEL_SNP > cap=$P_MAX_FOR_RAUG; writing 64G stub"
    TMP_STUB=$OUT_DIAG.tmp.$$
    printf "phenotype\tsv_locus\tarm\tchrom\tstatus\trescue_group\tn_panel_snp\tsv_pip\ttop_snp_pip\tsum_pip\n" >  "$TMP_STUB"
    printf "%s\t%s\t%s\t%s\tR_AUG_TOO_LARGE_AT_64G\t%s\t%d\t\t\t\n" "$PHENO" "$SV_LOCUS" "$ARM" "$SV_CHROM" "$RESCUE_GROUP" "$N_PANEL_SNP" >> "$TMP_STUB"
    mv "$TMP_STUB" "$OUT_DIAG"
    n_too_wide=$((n_too_wide + 1))
    touch "$HEARTBEAT_FILE"
    continue
  fi

  PHENO_TMP=$TMP_BASE/${SV_LOCUS}__${ARM}
  mkdir -p "$PHENO_TMP"

  echo "[prep ] master_pair_id=$PAIR_ID  $SV_LOCUS  $ARM  p=$N_PANEL_SNP  grp=$RESCUE_GROUP"
  $PY $CODE_DIR_REF/prepare_jwbundle_pair_generic_v3.py \
      "$REF_BUNDLE" \
      "$GWAS_PATH" \
      "$PHENO_TMP" \
      --r-diagonal-mode wRw \
      --r-diagonal-multiplier 1.0

  META=$PHENO_TMP/meta.tsv
  if [ ! -s "$META" ]; then
    echo "[FATAL] master_pair_id=$PAIR_ID prepare wrote no meta.tsv at $META"
    exit 6
  fi

  STATUS=$(awk -F'\t' '$1=="status" {print $2; exit}' "$META")
  if [ "$STATUS" = "ORIENT_GAP" ]; then
    TMP_STUB=$OUT_DIAG.tmp.$$
    printf "phenotype\tsv_locus\tarm\tchrom\tstatus\trescue_group\tsv_pip\ttop_snp_pip\tsum_pip\n" >  "$TMP_STUB"
    printf "%s\t%s\t%s\t%s\tORIENT_GAP\t%s\t\t\t\n" "$PHENO" "$SV_LOCUS" "$ARM" "$SV_CHROM" "$RESCUE_GROUP" >> "$TMP_STUB"
    mv "$TMP_STUB" "$OUT_DIAG"
    n_gap=$((n_gap + 1))
    rm -rf "$PHENO_TMP"
    touch "$HEARTBEAT_FILE"
    continue
  fi
  if [ "$STATUS" != "OK" ]; then
    echo "[FATAL] master_pair_id=$PAIR_ID unexpected meta.tsv status=$STATUS"
    exit 7
  fi

  echo "[susie] master_pair_id=$PAIR_ID  $SV_LOCUS  $ARM  p=$N_PANEL_SNP"
  $RSCRIPT $V2/R/run_susie_rss.R \
      --bundle-dir "$PHENO_TMP" \
      --pheno      "$PHENO" \
      --sv-locus   "$SV_LOCUS" \
      --arm        "$ARM" \
      --chrom      "$SV_CHROM" \
      --out-dir    "$SWEEP_ROOT"

  WRITTEN=$SWEEP_ROOT/$PHENO/$SV_CHROM/${SV_LOCUS}__${ARM}__diag.tsv
  if [ ! -s "$WRITTEN" ]; then
    echo "[FATAL] master_pair_id=$PAIR_ID SuSiE produced no diag.tsv at $WRITTEN"
    exit 8
  fi
  mv "$WRITTEN" "$OUT_DIAG"

  rm -rf "$PHENO_TMP"
  n_done=$((n_done + 1))
  touch "$HEARTBEAT_FILE"
done

echo "=========================================="
echo "[$(date '+%F %T')] phaseB_rescue task=$SLURM_ARRAY_TASK_ID DONE  done=$n_done  orient_gap=$n_gap  too_wide=$n_too_wide  skipped=$n_skip"
echo "=========================================="
