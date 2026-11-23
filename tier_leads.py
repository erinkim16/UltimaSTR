"""Step 3: sort leads by strength of evidence.
Reads <stem>.filtered.tsv (filter_leads.py) and site_quality.tsv (site_quality.py),
looks at the LARGEST allele of each lead (the one that triggered the flag) and assigns:
  tier 1 strong   exact (spanning) count, >=3 spanning reads, >=3 units above normal_max
  tier 2 moderate exact but 1-2 reads, or >=3 above normal with fewer than 3 reads;
                  or a lower bound >= pathogenic_min at a site that is usually spanned
                  (>=80% fully spanned), so missing spanning reads is itself unusual
  tier 3 weak     only 1-2 units above normal_max (common small variation),
                  or a lower bound at a site that is often not spanned
Cohort-outlier leads (John's sites, which have no published thresholds) are
judged without thresholds and are never tier 1, since no disease size is known:
  tier 2  exact with >=3 spanning reads, or a lower bound at a site spanned in >=80%
  tier 3  the rest (exact with 1-2 reads, or a lower bound at a hard-to-read site)
Writes <stem>.tiered.tsv (leads only, tier column + reason), sorted by tier.
"""
import argparse, csv, re
ap = argparse.ArgumentParser()
ap.add_argument("--table", default="vincent_str_calls.filtered.tsv")
ap.add_argument("--quality", default="site_quality.tsv")
a = ap.parse_args()
q = {r["locus"]: float(r["fully_spanned_pct"]) for r in csv.DictReader(open(a.quality), delimiter="\t")}
rows = [r for r in csv.DictReader(open(a.table), delimiter="\t") if r["lead"] == "yes"]
nums = lambda s: [int(x) for x in re.findall(r"\d+", s)]
out = []
for r in rows:
    gts = r["genotype"].strip("()").split(",")
    sup = r["support"].split("/"); sp = nums(r["spanning_reads"])
    i = max(range(len(gts)), key=lambda k: nums(gts[k])[0] if nums(gts[k]) else 0)  # largest allele
    val = nums(gts[i])[0]; exact = sup[i] == "SPANNING" and not gts[i].startswith(">=")
    reads = sp[i] if i < len(sp) else 0
    nm = int(r["normal_max"]) if r["normal_max"] else None
    pm = int(r["pathogenic_min"]) if r["pathogenic_min"] else None
    excess = val - nm if nm is not None else 0
    spanned = q.get(r["locus"], 0)
    if r.get("lead_type") == "cohort_outlier":
        if exact and reads >= 3:
            t, why = 2, f"no thresholds; exact {val} with {reads} spanning reads, in cohort top 5%"
        elif not exact and spanned >= 80:
            t, why = 2, f"no thresholds; >= {val} with no spanning reads at a site spanned in {spanned:.0f}% of people, in cohort top 5%"
        else:
            t, why = 3, f"no thresholds; {'exact' if exact else '>='} {val} in cohort top 5%, but " + (
                f"only {reads} spanning read(s)" if exact else f"site spanned in only {spanned:.0f}% of people")
    elif exact and reads >= 3 and excess >= 3:
        t, why = 1, f"exact {val}, {reads} spanning reads, {excess} above normal"
    elif exact and excess >= 3:
        t, why = 2, f"exact {val} but only {reads} spanning read(s)"
    elif not exact and pm is not None and val >= pm and spanned >= 80:
        t, why = 2, f">= {val} (pathogenic range) with no spanning reads at a site spanned in {spanned:.0f}% of people"
    elif exact:
        t, why = 3, f"only {excess} above normal max {nm} (small, common variation)"
    else:
        t, why = 3, f"lower bound >= {val}; site spanned in only {spanned:.0f}% of people" if spanned < 80 else f"lower bound >= {val}, below pathogenic range"
    r["tier"], r["tier_reason"] = t, why
    out.append(r)
out.sort(key=lambda r: (r["tier"], r["cohort_flag_rate"], r["locus"]))
stem = a.table.replace(".filtered.tsv", "")
cols = ["tier", "sample", "locus", "genotype", "support", "spanning_reads", "normal_max", "pathogenic_min",
        "flag", "cohort_flag_rate", "cohort_top5", "lead_type", "tier_reason"]
with open(stem + ".tiered.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore", lineterminator="\n")
    w.writeheader(); w.writerows(out)
for t in (1, 2, 3): print(f"tier {t}: {sum(r['tier'] == t for r in out)}")
print("wrote", stem + ".tiered.tsv")
