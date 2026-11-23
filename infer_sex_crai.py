"""Infer sex from CRAM index (.crai) files only, like indexcov: no CRAM reading needed.

A .crai is gzipped text, one line per slice: ref_id, start, span, container_offset,
slice_offset, slice_size. Bytes of sequence per chromosome track read depth, so
  x_ratio = (bytes/bp on chrX) / median(bytes/bp on chr1-22)   ~1.0 female, ~0.5 male
  y_ratio = same for chrY                                        ~0.15 female (X-homologous reads), ~0.5 male
ref_id is the index into the CRAM header's @SQ list, assumed to follow the reference .fai order.

usage: python infer_sex_crai.py --fai ref.fasta.fai --out sex.tsv  sample_id=path.crai ...
   or: python infer_sex_crai.py --fai ref.fasta.fai --samples samples.tsv   (col1 id, col3 crai)
"""
import argparse, collections, gzip, statistics, sys

ap = argparse.ArgumentParser()
ap.add_argument("--fai", required=True)
ap.add_argument("--samples", help="TSV: sample_id, cram, crai[, ...] (samplesheet layout)")
ap.add_argument("--out", default="-")
ap.add_argument("pairs", nargs="*", help="sample_id=/path/to/file.crai")
args = ap.parse_args()

fai = [l.split("\t")[:2] for l in open(args.fai)]
names, lengths = [n for n, _ in fai], [int(L) for _, L in fai]
auto = {f"chr{i}" for i in range(1, 23)}

jobs = [p.split("=", 1) for p in args.pairs]
if args.samples:
    jobs += [(r[0], r[2]) for r in (l.rstrip("\n").split("\t") for l in open(args.samples)) if r[0]]
if not jobs:
    sys.exit("no samples given")

out = sys.stdout if args.out == "-" else open(args.out, "w")
out.write("sample_id\tx_ratio\ty_ratio\tinferred_sex\tcrai\n")
for sid, crai in jobs:
    size = collections.Counter()
    with gzip.open(crai, "rt") as f:
        for line in f:
            ref, _, _, _, _, sz = line.split("\t")
            ref = int(ref)
            if 0 <= ref < len(names):
                size[names[ref]] += int(sz)
    dens = {n: size[n] / L for n, L in zip(names, lengths) if n in auto | {"chrX", "chrY"}}
    med = statistics.median(dens[n] for n in auto)
    x, y = dens["chrX"] / med, dens["chrY"] / med
    # cut-offs calibrated on Vincent_ASD (41 people, 100% agreement with indexcov): F x 0.95-0.99, y 0.12-0.26; M x 0.46-0.52, y 0.40-0.63
    sex = "F" if x > 0.8 and y < 0.3 else "M" if 0.35 < x < 0.65 and y > 0.3 else "unknown"
    out.write(f"{sid}\t{x:.3f}\t{y:.3f}\t{sex}\t{crai}\n")
    out.flush()
