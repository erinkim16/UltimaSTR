#!/bin/bash
# Copy per-sample STR outputs out of a miniwdl run dir into a persistent results folder.
# Usage: collect_str_outputs.sh <run_dir> <dest_dir>
set -euo pipefail
RUN=${1:-_LAST}
DEST=${2:?give a destination, e.g. /project/rrg-shreejoy/<you>/Vincent_ASD/STR_Results}
mkdir -p "$DEST"
find -L "$RUN/out" -type f \( -name '*_genotypes.vcf.gz' -o -name '*_genotypes.vcf.gz.tbi' \
     -o -name '*_genotypes.bed' -o -name '*_summary.csv' -o -name '*_detailed.csv' \) \
     -exec cp -L -n {} "$DEST"/ \;
echo "$(ls "$DEST"/*_genotypes.vcf.gz | wc -l) samples in $DEST"
