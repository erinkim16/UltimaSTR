"""Vincent_ASD samplesheet from Ellie's sample_sex.tsv: sample_id  cram  crai  sex(M|F)."""
import argparse, csv, os, sys

ap = argparse.ArgumentParser()
ap.add_argument("--sex_table", default="/scratch/ehogan/wgs_processing/SNV/results/Vincent_ASD_v2/sample_sex.tsv")
ap.add_argument("--exclude", action="append", default=[], help="e.g. ULTWGS_0042")
ap.add_argument("--out", default="samples.tsv")
args = ap.parse_args()

rows, held = [], []
for r in csv.DictReader(open(args.sex_table), delimiter="\t"):
    sid = "_".join(r["sample_id"].split("_")[:2])          # ULTWGS_0001_01_SG01 -> ULTWGS_0001
    x, y = float(r["chrX_auto_ratio"]), float(r["chrY_auto_ratio"])
    if sid in args.exclude:
        held.append(f"{sid}: excluded on request"); continue
    # label must agree with the X/Y coverage ratios
    if r["inferred_sex"] == "female" and x > 0.8 and y < 0.2:
        sex = "F"
    elif r["inferred_sex"] == "male" and 0.35 < x < 0.65 and y > 0.25:
        sex = "M"
    else:
        held.append(f"{sid}: label {r['inferred_sex']} disagrees with X={x} Y={y}"); continue
    cram = r["cram"]
    crai = next((c for c in (cram + ".crai", cram[:-5] + ".crai") if os.path.isfile(c)), None)
    if not os.path.isfile(cram) or crai is None:
        held.append(f"{sid}: missing {cram} or its index"); continue
    rows.append((sid, cram, crai, sex, r["cram_source"]))

with open(args.out, "w") as f:
    f.writelines("\t".join(r[:4]) + "\n" for r in rows)
n_m = sum(r[3] == "M" for r in rows)
print(f"wrote {len(rows)} samples to {args.out} ({n_m} M, {len(rows) - n_m} F; "
      f"{sum(r[4] == 'merged' for r in rows)} merged)")
for h in held:
    print("HELD BACK:", h, file=sys.stderr)