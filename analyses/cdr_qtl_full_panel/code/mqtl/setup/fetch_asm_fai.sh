#!/usr/bin/env bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -uo pipefail
P1=${CENSV_ROOT}
SUB=$P1/subproject_methylation
DEST=$SUB/raw/asm_fai
MAN=$SUB/manifest/methylation_manifest.tsv
PRE_ROOT="s3://human-pangenomics/submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION"
mkdir -p "$DEST"
cd "$SUB"

ok=0; miss=0
while IFS=$'\t' read -r S HL HR _ _ _ VER _; do
  [ "$S" = "sample" ] && continue
  [ -z "${HR:-}" ] && HR="$HL"
  [ -z "${VER:-}" ] && VER="v1.0.1"
  OUT="$DEST/${S}_${HR}.fa.gz.fai"
  [ -s "$OUT" ] && { ok=$((ok+1)); continue; }
  GOT=0
  for V in "$VER" v1.0.1 v1.1 v1; do
    KEY="${S}_${HR}_hprc_r2_${V}.fa.gz.fai"
    if ~/.pixi/bin/pixi run aws s3 cp --no-sign-request --quiet \
         "$PRE_ROOT/$S/assemblies/freeze_2/$KEY" "$OUT" 2>/dev/null; then
      GOT=1; ok=$((ok+1)); break
    fi
  done
  if [ "$GOT" -eq 0 ]; then
    KEY=$(~/.pixi/bin/pixi run aws s3 ls "$PRE_ROOT/$S/assemblies/freeze_2/" --no-sign-request 2>/dev/null \
          | awk '{print $NF}' | grep -E "^${S}_${HR}_.*\.fa\.gz\.fai$" | tail -1)
    if [ -n "$KEY" ] && ~/.pixi/bin/pixi run aws s3 cp --no-sign-request --quiet \
         "$PRE_ROOT/$S/assemblies/freeze_2/$KEY" "$OUT" 2>/dev/null; then
      ok=$((ok+1))
    else
      echo "MISS ${S}_${HR}"; miss=$((miss+1))
    fi
  fi
done < "$MAN"
echo "fai on disk: $(ls "$DEST" | wc -l)   ok=$ok miss=$miss"
