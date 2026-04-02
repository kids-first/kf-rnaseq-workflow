#!/usr/bin/env bash

set -exo

GENCODE_VER=$1
CTAT_RESOURCE_VER=$2
CTAT_FUSION_VER=$3
HLA_VER=$4
DFAM_VER=$5
PFAM_VER=$6

usage() {
  echo "Usage: \$(basename "$0") <GENCODE_VER> <CTAT_RESOURCE_VER> <CTAT_FUSION_VER> <HLA_VER> <DFAM_VER> <PFAM_VER>" >&2
  echo "Example: \$(basename "$0") 39 GRCh38_gencode_v44_CTAT_lib_Oct292023 0.3.0 3.43.0-alpha 3.1 36.0" >&2
}

# exit 2 if any required arg is missing/empty
if [ -z "$GENCODE_VER" ] || [ -z "$CTAT_RESOURCE_VER" ] || [ -z "$CTAT_FUSION_VER" ] || [ -z "$HLA_VER" ] || [ -z "$DFAM_VER" ] || [ -z "$PFAM_VER" ]; then
  usage
  exit 2
fi

fetch_and_echo() {
  local url="$1"
  echo "Downloading $url"
  wget -q "$url"
  echo "$url" >> gencode_reference_manifest.txt;
}

fetch_and_echo "https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/\${CTAT_RESOURCE_VER}.source.tar.gz"

fetch_and_echo "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_\${GENCODE_VER}/gencode.v\${GENCODE_VER}.primary_assembly.annotation.gtf.gz"
fetch_and_echo "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_\${GENCODE_VER}/gencode.v\${GENCODE_VER}.transcripts.fa.gz"
fetch_and_echo "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_\${GENCODE_VER}/GRCh38.primary_assembly.genome.fa.gz"

curl -s "https://api.github.com/repos/FusionAnnotator/CTAT_HumanFusionLib/releases/tags/v\${CTAT_FUSION_VER}" \
  | jq -r '.assets[].browser_download_url' \
  | while read -r url; do
      fetch_and_echo "$url";
    done

fetch_and_echo "https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/AnnotFilterRule.pm"

fetch_and_echo "https://github.com/ANHIG/IMGTHLA/raw/refs/tags/v\${HLA_VER}/hla.dat"

exts=("" ".h3f" ".h3i" ".h3m" ".h3p")
for ext in "\${exts[@]}"; do
  fetch_and_echo "http://dfam.org/releases/Dfam_\${DFAM_VER}/infrastructure/dfamscan/homo_sapiens_dfam.hmm\${ext}";
done

fetch_and_echo "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam\${PFAM_VER}/Pfam-A.hmm.gz"
