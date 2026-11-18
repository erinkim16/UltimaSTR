"""What's normal in this cohort: per-site allele distribution, exact (spanning) alleles only.
Writes cohort_distributions.tsv: one row per site with the spread of exactly measured alleles
next to Broad's thresholds, and how many alleles sit above normal_max.
"""
import argparse, csv, collections, re, statistics as st
ap = argparse.ArgumentParser(); ap.add_argument("--table", default="vincent_str_calls.tsv")
ap.add_argument("--out", default="cohort_distributions.tsv"); a = ap.parse_args()
exact, lower, info = collections.defaultdict(list), collections.defaultdict(list), {}
for r in csv.DictReader(open(a.table), delimiter="\t"):
    info[r["locus"]] = (r["normal_max"], r["pathogenic_min"])
    gts = r["genotype"].strip("()").split(",") if r["genotype"] else []
    for g, s in zip(gts, r["support"].split("/")):
        n = int(re.sub(r"\D", "", g) or 0)
        (exact if s == "SPANNING" and not g.startswith(">=") else lower)[r["locus"]].append(n)
def q(v, p): v = sorted(v); return v[int(p * (len(v) - 1))]
with open(a.out, "w") as f:
    f.write("locus\tn_exact\tn_lower_bound\tmin\tp5\tmedian\tp95\tmax\tmode\tnormal_max\tpathogenic_min\texact_above_normal\n")
    for l in sorted(info):
        v = exact[l]; nm, pm = info[l]
        if not v: f.write(f"{l}\t0\t{len(lower[l])}\t\t\t\t\t\t\t{nm}\t{pm}\t\n"); continue
        above = sum(x > int(nm) for x in v) if nm else ""
        f.write(f"{l}\t{len(v)}\t{len(lower[l])}\t{min(v)}\t{q(v,.05)}\t{st.median(v):g}\t{q(v,.95)}\t{max(v)}\t"
                f"{collections.Counter(v).most_common(1)[0][0]}\t{nm}\t{pm}\t{above}\n")
print("wrote", a.out)
