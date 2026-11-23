#!/bin/bash
# slice_str.sh: cut a small BAM of the reads around one STR, for viewing in JBrowse 2.
# Usage: bash slice_str.sh <sample> <locus> <cram> [pad=500]
#   e.g. bash slice_str.sh ULTWGS_0002 HTT /path/to/ULTWGS_0002.cram
# Writes slices/<sample>_<locus>.bam (+ .bai) and slices/<locus>.bed (the repeat itself).
# These BAMs contain participant reads: keep them in the lab (gitignored, not on slides unanonymised).
set -euo pipefail
SAMPLE=$1; LOCUS=$2; CRAM=$3; PAD=${4:-500}
REF=$(ls references/*.fasta | head -1)
CAT=$(ls catalogs/*.female.json | head -1)
mkdir -p slices

# repeat coordinates from the catalog (first ReferenceRegion of the locus)
read CHR START END < <(python3 - "$CAT" "$LOCUS" <<'EOF'
import json, sys
for l in json.load(open(sys.argv[1])):
    if l["LocusId"] == sys.argv[2]:
        r = l["ReferenceRegion"]; r = r[0] if isinstance(r, list) else r
        c, se = r.split(":"); s, e = se.split("-"); print(c, s, e); break
else:
    sys.exit("locus not found in catalog")
EOF
)
echo "$LOCUS repeat: $CHR:$START-$END (padding $PAD)"
printf "%s\t%s\t%s\t%s\n" "$CHR" "$START" "$END" "$LOCUS" > "slices/$LOCUS.bed"

OUT="slices/${SAMPLE}_${LOCUS}.bam"
samtools view -b -T "$REF" -o "$OUT" "$CRAM" "$CHR:$((START-PAD))-$((END+PAD))"
samtools index "$OUT"
echo "wrote $OUT ($(samtools view -c "$OUT") reads); view at $CHR:$((START-PAD))-$((END+PAD))"
