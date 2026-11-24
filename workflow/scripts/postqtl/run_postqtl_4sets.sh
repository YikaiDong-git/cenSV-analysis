#!/bin/bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail
cd ${CENSV_ROOT}
PIXIRUN=~/.pixi/bin/pixi
SETS=(cv01_maf005 cv01_maf010 cv03_maf005 cv03_maf010)
for SET in "${SETS[@]}"; do
  echo ">>> processing SET=$SET"
  ln -sfnT cosine_nj_${SET} results/eqtl_overlap/decompose/cosine_nj
  ln -sfnT cosine_nj_${SET} results/sqtl_overlap/decompose/cosine_nj
  $PIXIRUN run python3 workflow/scripts/postqtl/trans_bh.py       2>&1 | tail -5
  $PIXIRUN run python3 workflow/scripts/postqtl/list_all_genes.py 2>&1 | tail -5
  $PIXIRUN run python3 workflow/scripts/postqtl/merge_gene_annotations.py 2>&1 | tail -5
  $PIXIRUN run python3 workflow/scripts/postqtl/make_summary_table.py 2>&1 | tail -3
  DEST=results/analysis_qtl/final_tables_${SET}
  rm -rf $DEST
  mv results/analysis_qtl/final_tables $DEST
  mv results/postqtl_v2b results/postqtl_v2b_${SET}
  echo ">>> DONE $SET -> $DEST"
done
rm -f results/eqtl_overlap/decompose/cosine_nj
rm -f results/sqtl_overlap/decompose/cosine_nj
echo ALL_DONE
