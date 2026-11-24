#!/bin/bash
#SBATCH --job-name=nullmatch
#SBATCH --output=${CENSV_ROOT}/analyses/meqtl_replication/logs/b2_cpg_universe/nullmatch_%j.log
#SBATCH --time=2:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -uo pipefail
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2

P=${CENSV_ROOT}
U=$P/analyses/meqtl_replication

echo "HOST $(hostname)  START $(date -u +%Y-%m-%dT%H:%M:%SZ)"

NDONE=$(ls $U/output/null_design/cpg_universe/*/*.cpg.tsv.gz 2>/dev/null | wc -l)
NWANT=$(wc -l < $U/output/null_design/cpg_universe/_tasklist.txt)
echo "universe chunks: $NDONE / $NWANT"
if [ "$NDONE" -lt "$NWANT" ]; then
  echo "FAIL universe incomplete -- refusing to build a null from a partial universe"
  exit 1
fi

cd "$P" || exit 1
~/.pixi/bin/pixi run python "$U/code/04_null/build_matched_null.py"
rc=$?
echo "EXIT $rc  END $(date -u +%Y-%m-%dT%H:%M:%SZ)"
if [ "$rc" -eq 0 ]; then echo "PASS"; else echo "FAIL rc=$rc"; fi
exit $rc
