#!/usr/bin/env bash
: "${CENSV_ROOT:?set CENSV_ROOT (see README)}"
set -euo pipefail

PROJECT_ROOT=${CENSV_ROOT}
DEST="${PROJECT_ROOT}/resources/hg002_q100"
URL="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz"
URL_FAI="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz.fai"
URL_GZI="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz.gzi"
URL_MD5="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz.md5"

mkdir -p "${DEST}"

ssh dtn01 "cd ${DEST} && \
  nohup bash -c '
    set -e
    curl -L -o hg002v1.1.fasta.gz.md5 ${URL_MD5}
    curl -L -o hg002v1.1.fasta.gz.fai ${URL_FAI}
    curl -L -o hg002v1.1.fasta.gz.gzi ${URL_GZI}
    curl -L -o hg002v1.1.fasta.gz ${URL}
    md5sum -c hg002v1.1.fasta.gz.md5
    echo \"[done] \$(date -Iseconds)\"
  ' > download.log 2>&1 &
  echo \"[dtn01 launched pid \$!]\"
"
echo "Download dispatched on dtn01. Check ${DEST}/download.log for progress."
