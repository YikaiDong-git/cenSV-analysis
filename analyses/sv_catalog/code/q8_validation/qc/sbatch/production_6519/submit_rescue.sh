#!/bin/bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"

set -euo pipefail

PROJ=${CENSV_ROOT}/analyses/sv_catalog
SBATCH_DIR=$PROJ/code/q8_validation/qc/sbatch/production_6519
PHASE_B_RESCUE=$SBATCH_DIR/sbatch_phase_b_rescue.sh
MANIFEST_ROOT=$PROJ/output/ref_compare_proxy_vs_dense/sweep_6519_locked_v3/manifest

K=2
THROTTLE=50

GWAS_FG=${CENSV_ROOT}/centromere_twas/finngen/processed_b38aligned/I9_HYPTENSESS/I9_HYPTENSESS.hg38_b38aligned.tsv
GWAS_NEALE=${CENSV_WORKSPACE}/_shared/datasets/GWAS/neale_ukbb_round2/v20180801/processed_b38aligned/30210_irnt/30210_irnt.hg38_b38aligned.tsv

DRY=0
for arg in "$@"; do
  case "$arg" in
    --dry-run) DRY=1 ;;
    *) echo "[FATAL] unknown arg: $arg" >&2; exit 2 ;;
  esac
done

require_numeric_jid() {
  local label=$1
  local jid=$2
  if ! [[ "$jid" =~ ^[0-9]+$ ]]; then
    echo "[FATAL] $label submit returned non-numeric JID: '$jid'" >&2
    exit 9
  fi
}

existing_queue=$(squeue -h -u "$USER" -r 2>/dev/null | wc -l)

count_tasks() {
  local mani=$1
  local rows=$(($(wc -l < "$mani") - 1))   # subtract header
  echo $(( (rows + K - 1) / K ))           # ceil(rows/K)
}

MANI_FG=$MANIFEST_ROOT/MISSING_RESCUE_I9_HYPTENSESS.tsv
MANI_NEALE=$MANIFEST_ROOT/MISSING_RESCUE_30210_irnt.tsv
for m in "$MANI_FG" "$MANI_NEALE"; do
  if [ ! -f "$m" ]; then
    echo "[FATAL] rescue manifest missing: $m" >&2
    echo "        Run build_rescue_manifest.py first." >&2
    exit 3
  fi
done

N_FG=$(count_tasks "$MANI_FG")
N_NEALE=$(count_tasks "$MANI_NEALE")
projected=$((existing_queue + N_FG + N_NEALE))
echo "[preflight] current squeue=$existing_queue  FG_tasks=$N_FG  Neale_tasks=$N_NEALE  projected=$projected  cap=1000"
if [ "$projected" -gt 1000 ]; then
  echo "[FATAL] projected queue $projected > 1000 cap; aborting" >&2
  exit 4
fi

submit_pheno() {
  local PHENO=$1
  local STUDY=$2
  local GWAS=$3
  local N_TASKS=$4
  echo "[submit] rescue — $STUDY/$PHENO  GWAS=$GWAS  N_tasks=$N_TASKS%${THROTTLE}" >&2
  if [ ! -f "$GWAS" ]; then
    echo "[FATAL] GWAS missing: $GWAS" >&2
    exit 5
  fi
  local EXPORT="ALL,PHENO=$PHENO,STUDY=$STUDY,GWAS_PATH=$GWAS"
  if [ "$DRY" = "1" ]; then
    echo "[dry] sbatch --parsable --array=1-${N_TASKS}%${THROTTLE} --export=$EXPORT $PHASE_B_RESCUE" >&2
    echo "DRY_RESCUE_$PHENO"
  else
    local jid
    jid=$(sbatch --parsable --array=1-${N_TASKS}%${THROTTLE} --export="$EXPORT" "$PHASE_B_RESCUE")
    jid=${jid%%;*}
    require_numeric_jid "rescue/$PHENO" "$jid"
    echo "$jid"
  fi
}

JID_FG=$(submit_pheno I9_HYPTENSESS finngen "$GWAS_FG" "$N_FG")
echo "[submit] JID_RESCUE_FG=$JID_FG"

JID_NEALE=$(submit_pheno 30210_irnt neale "$GWAS_NEALE" "$N_NEALE")
echo "[submit] JID_RESCUE_NEALE=$JID_NEALE"

echo
echo "=========================================="
echo "RESCUE SUBMISSION SUMMARY"
echo "  Rescue — I9_HYPTENSESS (FG)    : $JID_FG  (array=1-${N_FG}%${THROTTLE})"
echo "  Rescue — 30210_irnt (Neale)    : $JID_NEALE  (array=1-${N_NEALE}%${THROTTLE})"
echo "  K=$K cells/task   --mem=64G   --time=6:00:00   P_MAX_FOR_RAUG=50000"
echo "  Output tree: output/ref_compare_proxy_vs_dense/sweep_6519_locked_v3/rescue/<pheno>/"
echo "=========================================="
echo
echo "Track with:"
echo "  squeue -j $JID_FG,$JID_NEALE"
