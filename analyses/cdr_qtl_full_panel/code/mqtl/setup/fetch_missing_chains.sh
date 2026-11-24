#!/usr/bin/env bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -uo pipefail
P1=${CENSV_ROOT}
SUB=$P1/subproject_methylation
WU=$P1/analyses/cdr_qtl_full_panel
MIS=$WU/output/mqtl/missing_tags.txt
DEST=$SUB/raw/chains_chm13
PRE_ROOT="s3://human-pangenomics/submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION"

mkdir -p "$DEST"
ok=0; fail=0
cut -d_ -f1 "$MIS" | sort -u | while read -r S; do
  [ -z "$S" ] && continue
  LIST=$(~/.pixi/bin/pixi run aws s3 ls "$PRE_ROOT/$S/assemblies/freeze_2/annotation/chains/minigraph-cactus/" \
         --no-sign-request 2>/dev/null | awk '{print $NF}' | grep -E '_vs_CHM13\.chain\.gz$' | sort)
  if [ -z "$LIST" ]; then echo "NO_CHAINS_IN_S3 $S"; continue; fi
  for KEY in $LIST; do
    OUT="$DEST/${KEY%%_hprc_r2*}_vs_CHM13.chain.gz"     # e.g. HG01081_mat_vs_CHM13.chain.gz
    if [ -s "$OUT" ]; then echo "have $(basename "$OUT")"; continue; fi
    if ~/.pixi/bin/pixi run aws s3 cp --no-sign-request --quiet \
         "$PRE_ROOT/$S/assemblies/freeze_2/annotation/chains/minigraph-cactus/$KEY" "$OUT" 2>/dev/null; then
      if gzip -t "$OUT" 2>/dev/null; then echo "OK $(basename "$OUT")"; else echo "CORRUPT $(basename "$OUT")"; rm -f "$OUT"; fi
    else
      echo "DLFAIL $KEY"
    fi
  done
done
echo "--- chains on disk now: $(ls "$DEST" | wc -l) ---"
