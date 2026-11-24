#!/bin/bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -u
c="$1"
P1=${CENSV_ROOT}
PY=$P1/.pixi/envs/default/bin/python
SC=$P1/analyses/cdr_qtl/code
for ph in mean_5mC n_cdrs cdrs_per_mb cdr_fraction median_cdr_width cdr_pos_offset cdr_pos_var cdr_pos_iqr; do
  echo "=== $c $ph ==="
  "$PY" "$SC/b6_cdr_pos_qtl_perchrom.py" \
    --chrom "$c" --phenotype "$ph" --n-perm 100000 --cohort all || echo "FAILED $c $ph"
done
echo "ALL_PHENOS_DONE $c"
