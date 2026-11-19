"""Cohort filter for the table written by summarize_str.py.

Adds three columns and writes a short, sorted leads list.
  cohort_flag_rate  share of people flagged (PATHOGENIC_RANGE or INTERMEDIATE)
                    at this site. If most of the cohort is flagged, the
                    threshold or the counting doesn't fit the site
                    (PRE-MIR7-2, RAI1, PABPN1 sub-variant, RFC1, ...), so its
                    flags are demoted, not deleted.
  cohort_top5       yes if this person's largest allele is in the top 5% of
                    all alleles seen at this site in the cohort.
  lead              yes if flagged, at a site with flag rate <= --max-rate,
                    or, at the --cohort-sites (default: John's CA10, DIP2B_CAC,
                    PTCHD1, which have no published thresholds), if the call is
                    in the cohort top 5%. Pass --cohort-sites "" to turn it off.
  lead_type         threshold (flag vs. published thresholds) or
                    cohort_outlier (no thresholds; unusual vs. this cohort).
Leads file is sorted: top-5% calls first, then by spanning-read support.
Note: lower bounds (>=N) are treated as N, so a flanking-only call can
still be a lead; check its spanning_reads before believing it.
"""
import argparse, csv, collections, re

ap = argparse.ArgumentParser()
ap.add_argument("--table", default="vincent_str_calls.tsv")
ap.add_argument("--max-rate", type=float, default=0.5)
ap.add_argument("--cohort-sites", default="CA10,DIP2B_CAC,PTCHD1",
                help="comma-separated sites judged against the cohort (top 5%%) instead of thresholds")
args = ap.parse_args()
stem = args.table[:-4] if args.table.endswith(".tsv") else args.table
out, leads = stem + ".filtered.tsv", stem + ".leads.tsv"

cohort_sites = {x.strip() for x in args.cohort_sites.split(",") if x.strip()}
rows = list(csv.DictReader(open(args.table), delimiter="\t"))
alleles = lambda r: [int(x) for x in re.findall(r"\d+", r["genotype"])]
top = lambda r: max(alleles(r), default=0)
flagged = lambda r: "PATHOGENIC_RANGE" in r["flag"] or "INTERMEDIATE" in r["flag"]

n, k, pool = collections.Counter(), collections.Counter(), collections.defaultdict(list)
for r in rows:
    n[r["locus"]] += 1
    k[r["locus"]] += flagged(r)
    pool[r["locus"]] += alleles(r)
cut = {l: sorted(v)[int(0.95 * (len(v) - 1))] for l, v in pool.items() if v}

for r in rows:
    rate = k[r["locus"]] / n[r["locus"]]
    r["cohort_flag_rate"] = f"{rate:.2f}"
    r["cohort_top5"] = "yes" if alleles(r) and top(r) > cut[r["locus"]] else ""
    if flagged(r) and rate <= args.max_rate:
        r["lead"], r["lead_type"] = "yes", "threshold"
    elif r["locus"] in cohort_sites and r["cohort_top5"]:
        r["lead"], r["lead_type"] = "yes", "cohort_outlier"
    else:
        r["lead"], r["lead_type"] = "", ""

cols = list(rows[0].keys())
support = lambda r: max([int(x) for x in re.findall(r"\d+", r["spanning_reads"])], default=0)
L = sorted((r for r in rows if r["lead"]), key=lambda r: (not r["cohort_top5"], -support(r)))
for path, keep in ((out, rows), (leads, L)):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(keep)

noisy = sorted((l for l in n if k[l] / n[l] > args.max_rate), key=lambda l: -k[l] / n[l])
nt = sum(1 for r in L if r["lead_type"] == "threshold")
print(f"{sum(map(flagged, rows))} flagged calls -> {nt} threshold leads "
      f"({sum(1 for r in L if r['cohort_top5'] and r['lead_type'] == 'threshold')} also in the cohort top 5%)")
print(f"+ {len(L) - nt} cohort-outlier leads at the cohort sites:",
      ", ".join(sorted({r['locus'] for r in L if r['lead_type'] == 'cohort_outlier'})))
print(f"demoted sites (flag rate > {args.max_rate}):", ", ".join(f"{l} {k[l]}/{n[l]}" for l in noisy))
print("wrote", out, "and", leads)
