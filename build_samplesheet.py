"""Samplesheet for multi_sample_str_genotyper.wdl: sample_id  cram  crai  sex(M|F), no header.

One row per sample ID: the merged CRAM if one exists, else the single raw run.
IDs with several runs but no merged CRAM, unknown sex, or unreadable files are held back
(listed on stderr) instead of being guessed.

TAY:
  python build_samplesheet.py \
    --merged_dir /scratch/ehogan/wgs_processing/SNV/results/TAY/merged-samples \
    --runs_list /scratch/ehogan/wgs_processing/SNV/results/TAY/variant_calling_completed.txt \
    --sex_table runs/TAY/sex.tsv --out runs/TAY/samples.tsv
--sex_table: TSV with a header holding sample_id and a sex column (--sex_col, default
inferred_sex), values M/F/male/female. infer_sex_crai.py output works as-is.
"""
import argparse, collections, csv, glob, os, re, sys

ap = argparse.ArgumentParser()
ap.add_argument("--merged_dir")
ap.add_argument("--runs_list", action="append", default=[], help="file of raw run CRAM paths")
ap.add_argument("--id_regex", default=r"TAY01_CMH_\d+_\d+_[A-Z]+_\d+")
ap.add_argument("--sex_table", required=True)
ap.add_argument("--sex_col", default="inferred_sex")
ap.add_argument("--exclude", action="append", default=[])
ap.add_argument("--only_merged", action="store_true", help="skip single-run samples")
ap.add_argument("--out", default="samples.tsv")
args = ap.parse_args()
idre = re.compile(args.id_regex)

def sid_of(path):
    m = idre.search(os.path.basename(path))
    return m.group(0) if m else None

def index_of(cram):
    return next((c for c in (cram + ".crai", cram[:-5] + ".crai") if os.path.isfile(c)), None)

merged = {}
if args.merged_dir:
    for c in glob.glob(os.path.join(args.merged_dir, "*.cram")):
        merged[sid_of(c)] = c
runs = collections.defaultdict(set)
for lst in args.runs_list:
    for line in open(lst):
        p = line.strip()
        if p.endswith(".cram") and "_unmatched" not in p:
            runs[sid_of(p)].add(p)

SEX = {"m": "M", "male": "M", "f": "F", "female": "F"}
sex = {r["sample_id"]: SEX.get(r[args.sex_col].strip().lower())
       for r in csv.DictReader(open(args.sex_table), delimiter="\t")}

rows, held = [], []
for sid in sorted((set(merged) | set(runs)) - {None}):
    if sid in args.exclude:
        held.append(f"{sid}: excluded on request"); continue
    if sid in merged:
        cram, src = merged[sid], "merged"
    elif args.only_merged:
        continue
    elif len(runs[sid]) == 1:
        cram, src = next(iter(runs[sid])), "single"
    else:
        held.append(f"{sid}: {len(runs[sid])} runs but no merged CRAM"); continue
    crai = index_of(cram)
    if crai is None or not (os.access(cram, os.R_OK) and os.access(crai, os.R_OK)):
        held.append(f"{sid}: {cram} or its index missing/unreadable"); continue
    if sex.get(sid) is None:
        held.append(f"{sid}: no M/F sex in {args.sex_table}"); continue
    rows.append((sid, cram, crai, sex[sid], src))

with open(args.out, "w") as f:
    f.writelines("\t".join(r[:4]) + "\n" for r in rows)
c = collections.Counter((r[3], r[4]) for r in rows)
print(f"wrote {len(rows)} samples to {args.out}: "
      + ", ".join(f"{n} {s}/{src}" for (s, src), n in sorted(c.items())) + f"; held back {len(held)}")
for h in held:
    print("HELD BACK:", h, file=sys.stderr)
