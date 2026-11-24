#!/bin/bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"

set -euo pipefail

PROJ=${CENSV_ROOT}/analyses/sv_catalog
SBATCH_DIR=$PROJ/code/q8_validation/qc/sbatch/production_6519
MONITOR=$SBATCH_DIR/monitor_progress.sh

PHASE_A=$SBATCH_DIR/sbatch_phase_a_build_reference_bundles.sh
PHASE_B=$SBATCH_DIR/sbatch_phase_b_prepare_susie.sh

GWAS_FG=${CENSV_ROOT}/centromere_twas/finngen/processed_b38aligned/I9_HYPTENSESS/I9_HYPTENSESS.hg38_b38aligned.tsv
GWAS_NEALE=${CENSV_WORKSPACE}/_shared/datasets/GWAS/neale_ukbb_round2/v20180801/processed_b38aligned/30210_irnt/30210_irnt.hg38_b38aligned.tsv

DRY=0
SKIP_A=0
for arg in "$@"; do
  case "$arg" in
    --dry-run) DRY=1 ;;
    --skip-phase-a) SKIP_A=1 ;;
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
projected=$((existing_queue + 326 + 652))
echo "[preflight] current squeue=$existing_queue  projected (Phase A run + Phase B held)=$projected  cap=1000"
if [ "$projected" -gt 1000 ] && [ "$SKIP_A" = "0" ]; then
  echo "[WARN] projected queue $projected > 1000 cap. Will rely on Phase A draining before Phase B PD lands."
fi

if [ "$SKIP_A" = "0" ]; then
  echo "[submit] Phase A — build 6,519 reference bundles (652 tasks %100)" >&2
  if [ "$DRY" = "1" ]; then
    echo "[dry] sbatch --parsable $PHASE_A" >&2
    JID_A=DRY_A
  else
    JID_A=$(sbatch --parsable "$PHASE_A")
    JID_A=${JID_A%%;*}   # strip any cluster suffix
    require_numeric_jid "Phase A" "$JID_A"
  fi
  echo "[submit] JID_A=$JID_A"
  DEP_FLAG="--dependency=afterany:$JID_A"
else
  echo "[submit] Phase A SKIPPED — Phase B will run with no SLURM dependency" >&2
  JID_A=NA
  DEP_FLAG=""
fi

submit_pheno() {
  local PHENO=$1
  local STUDY=$2
  local GWAS=$3
  echo "[submit] Phase B — $STUDY/$PHENO  GWAS=$GWAS" >&2
  if [ ! -f "$GWAS" ]; then
    echo "[FATAL] GWAS missing: $GWAS" >&2
    exit 3
  fi
  local EXPORT="ALL,PHENO=$PHENO,STUDY=$STUDY,GWAS_PATH=$GWAS"
  if [ "$DRY" = "1" ]; then
    echo "[dry] sbatch --parsable $DEP_FLAG --export=$EXPORT $PHASE_B" >&2
    echo "DRY_B_$PHENO"
  else
    local jid
    jid=$(sbatch --parsable $DEP_FLAG --export="$EXPORT" "$PHASE_B")
    jid=${jid%%;*}
    require_numeric_jid "Phase B/$PHENO" "$jid"
    echo "$jid"
  fi
}

JID_B_FG=$(submit_pheno I9_HYPTENSESS finngen "$GWAS_FG")
echo "[submit] JID_B_FG=$JID_B_FG"

JID_B_NEALE=$(submit_pheno 30210_irnt neale "$GWAS_NEALE")
echo "[submit] JID_B_NEALE=$JID_B_NEALE"

if [ "$DRY" = "0" ]; then
  echo "[cron] installing hourly :17 monitor entry (idempotent)" >&2
  CRON_LINE="17 * * * * $MONITOR >/dev/null 2>&1"
  (crontab -l 2>/dev/null | grep -v "$MONITOR" ; echo "$CRON_LINE") | crontab -
  echo "[cron] installed: $CRON_LINE" >&2
fi

echo
echo "=========================================="
echo "SUBMISSION SUMMARY"
echo "  Phase A (build bundles)        : $JID_A"
echo "  Phase B — I9_HYPTENSESS (FG)   : $JID_B_FG"
echo "  Phase B — 30210_irnt (Neale)   : $JID_B_NEALE"
echo "=========================================="
echo
echo "Track with:"
echo "  squeue -j $JID_A,$JID_B_FG,$JID_B_NEALE"
echo "  tail -f $PROJ/logs/qc/production_6519/STATUS_progress.log"
echo
echo "Cleanup (after completion):"
echo "  crontab -l | grep -v 'monitor_progress.sh' | crontab -"
