"""Per-site measurability from a summarize_str.py cohort table.
fully_spanned_pct: % of people whose every allele at the site is backed by spanning reads.
median_min_spanning: median over people of the spanning-read count of their weaker allele.
"""
import argparse, csv, collections, re, statistics as st
ap = argparse.ArgumentParser(); ap.add_argument("--table", default="vincent_str_calls.tsv")
ap.add_argument("--out", default="site_quality.tsv"); a = ap.parse_args()
by = collections.defaultdict(list)
for r in csv.DictReader(open(a.table), delimiter="\t"): by[r["locus"]].append(r)
res = []
for l, rs in by.items():
    full = sum(all(s == "SPANNING" for s in r["support"].split("/")) for r in rs)
    weak = [min([int(x) for x in re.findall(r"\d+", r["spanning_reads"])] or [0]) for r in rs]
    res.append((l, rs[0]["chrom"], len(rs), round(100 * full / len(rs)), st.median(weak)))
res.sort(key=lambda x: x[3])
with open(a.out, "w") as f:
    f.write("locus\tchrom\tn_called\tfully_spanned_pct\tmedian_min_spanning\n")
    for x in res: f.write("\t".join(map(str, x)) + "\n")
print(f"{sum(x[3] >= 80 for x in res)}/{len(res)} sites fully spanned in >=80% of people; "
      f"{sum(x[3] < 50 for x in res)} under 50%. Wrote {a.out}")
