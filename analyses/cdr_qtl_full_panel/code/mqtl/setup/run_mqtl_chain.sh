#!/usr/bin/env bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail
P1=${CENSV_ROOT}
WU=$P1/analyses/cdr_qtl_full_panel
cd "$P1"
R() { echo "### $*"; ~/.pixi/bin/pixi run python -u "$@"; }

R "$WU/code/mqtl/lift/cis_mqtl_aggregate.py" per_hap_bam
for W in 10000 50000 100000; do R "$WU/code/mqtl/scan/cis_mqtl_scan.py" "$W" _bam; done
for W in 10000 50000 100000; do R "$WU/code/mqtl/scan/cis_mqtl_escalate.py" "$W" _bam 1000000 0.005; done
R "$WU/code/mqtl/scan/cis_mqtl_calibration.py" _bam
echo "CHAIN COMPLETE"
