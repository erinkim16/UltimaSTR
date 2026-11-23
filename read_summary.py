#!/usr/bin/env python3
"""Step 4: look at the individual reads behind a lead.

The tool's <sample>_detailed.csv has one row per (read, candidate repeat count): it lines each
read up against several versions of the site (repeat_count_in_aux = 10, 15, 16, ... copies) and
records how well each fits (alignment_score; score_ratio = score / perfect score).

For each read this script picks the best-fitting count:
  spanning read  -> the single best count (exact)
  flanking read  -> every count long enough scores the same (the read can't see the far end),
                    so the start of that flat tail is the lower bound ">=N"
Reads whose best score_ratio < 0.85 are the ones the tool throws away; they're shown as 'dropped'.

usage (on Trillium, from the repo):
  python3 read_summary.py --results results/Vincent_ASD                # default lead list
  python3 read_summary.py --results results/Vincent_ASD --pair ULTWGS_0001:DIP2B_CAC
  python3 read_summary.py --results results/Vincent_ASD --cohort DIP2B_CAC   # one site, all people
  add --plot to also write read_plots/<sample>_<locus>.png (needs matplotlib)
Output: read_summary.tsv (one row per read) and a printed summary per lead.
"""
import argparse, csv, os, sys
from collections import defaultdict, Counter

LEADS = [  # step-4 leads (top tiers) plus a normal person at the same site for comparison
    ("ULTWGS_0012", "DIP2B_CAC"), ("ULTWGS_0033", "DIP2B_CAC"), ("ULTWGS_0001", "DIP2B_CAC"),
    ("ULTWGS_0025", "HTT"), ("ULTWGS_0029", "FXN"), ("ULTWGS_0021", "TBP"), ("ULTWGS_0001", "ZIC2"),
    ("ULTWGS_0024", "FOXL2"), ("ULTWGS_0024", "RUNX2"), ("ULTWGS_0024", "ARX_1"),
    ("ULTWGS_0024", "RUNX2__6-45422750-45422801-GCN"), ("ULTWGS_0040", "ARX_2"),
]
MIN_RATIO = 0.85

ap = argparse.ArgumentParser()
ap.add_argument("--results", required=True, help="folder with <sample>_detailed.csv")
ap.add_argument("--pair", action="append", help="SAMPLE:LOCUS (repeatable); replaces the default list")
ap.add_argument("--out", default="read_summary.tsv")
ap.add_argument("--plot", action="store_true")
ap.add_argument("--cohort", metavar="LOCUS",
                help="check one site in EVERY person: how many reads point past --long repeats")
ap.add_argument("--long", type=int, default=20, help="with --cohort: what counts as a long read (default 20)")
a = ap.parse_args()
if a.cohort:
    import glob
    paths = sorted(glob.glob(os.path.join(a.results, "*_detailed.csv")))
    pairs = [(os.path.basename(p).replace("_detailed.csv", ""), a.cohort) for p in paths]
else:
    pairs = [tuple(p.split(":", 1)) for p in a.pair] if a.pair else LEADS

by_sample = defaultdict(set)
for s, l in pairs: by_sample[s].add(l)

rows_out = []
for sample, loci in by_sample.items():
    path = os.path.join(a.results, f"{sample}_detailed.csv")
    if not os.path.exists(path):
        print(f"!! missing {path}", file=sys.stderr); continue
    reads = defaultdict(list)   # (locus, read_id, read_type) -> [(count, score, ratio)]
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            if r["variant_id"] in loci:
                reads[(r["variant_id"], r["read_id"], r["read_type"])].append(
                    (int(r["repeat_count_in_aux"]), float(r["alignment_score"]), float(r["score_ratio"])))
    for (locus, rid, rtype), cands in reads.items():
        best = max(c[1] for c in cands)
        tied = sorted(c[0] for c in cands if c[1] == best)
        ratio = max(c[2] for c in cands)
        tested = sorted({c[0] for c in cands})
        if rtype == "spanning" and len(tied) == 1:
            call = str(tied[0])
        elif rtype == "spanning":
            call = f"{tied[0]}-{tied[-1]} (tie)"
        else:
            # flanking: the read can't see the far end, so every long-enough count scores the same.
            # Lower bound = start of that flat tail (smallest count from which all longer ones tie).
            sc = sorted((c[0], c[1]) for c in cands)
            tail = sc[-1][1]; lb = sc[-1][0]
            for cnt, scr in reversed(sc):
                if scr != tail: break
                lb = cnt
            tied = [lb]
            call = f">={lb}"
        rows_out.append(dict(sample=sample, locus=locus, read_id=rid, read_type=rtype,
                             call=call, best_count=tied[0], score_ratio=f"{ratio:.3f}",
                             kept="yes" if ratio >= MIN_RATIO else "dropped",
                             tested_counts=",".join(map(str, tested))))

rows_out.sort(key=lambda r: (r["sample"], r["locus"], r["read_type"] != "spanning", r["best_count"]))
cols = ["sample", "locus", "read_type", "call", "score_ratio", "kept", "tested_counts", "read_id"]
with open(a.out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore", lineterminator="\n")
    w.writeheader(); w.writerows(rows_out)

groups = defaultdict(list)
for r in rows_out: groups[(r["sample"], r["locus"])].append(r)

if a.cohort:
    print(f"\n{a.cohort}: reads per person; 'long' = kept read whose best count / lower bound is > {a.long}")
    print("sample\treads\tkept\tspanning (top counts)\tlong_spanning\tlong_flanking\tlong_reads")
    for (s, l) in pairs:
        g = groups.get((s, l), []); kept = [r for r in g if r["kept"] == "yes"]
        sp = Counter(r["best_count"] for r in kept if r["read_type"] == "spanning")
        top = ", ".join(f"{k}x{n}" for k, n in sp.most_common(3))
        ls = [r for r in kept if r["read_type"] == "spanning" and r["best_count"] > a.long]
        lf = [r for r in kept if r["read_type"] != "spanning" and r["best_count"] > a.long]
        print(f"{s}\t{len(g)}\t{len(kept)}\t{top}\t{len(ls)}\t{len(lf)}\t"
              + " ".join(r["call"] for r in ls + lf))
    print(f"\nwrote {a.out} ({len(rows_out)} reads)")
    sys.exit(0)
for (s, l) in pairs:
    g = groups.get((s, l), [])
    kept = [r for r in g if r["kept"] == "yes"]
    sp = Counter(r["call"] for r in kept if r["read_type"] == "spanning")
    fl = Counter(r["call"] for r in kept if r["read_type"] != "spanning")
    print(f"\n{s} {l}: {len(g)} reads, {len(kept)} kept (score_ratio >= {MIN_RATIO}), "
          f"{len(g) - len(kept)} dropped")
    print("  spanning:", ", ".join(f"{k} x{n}" for k, n in sorted(sp.items(), key=lambda x: x[0])) or "none")
    print("  flanking:", ", ".join(f"{k} x{n}" for k, n in sorted(fl.items(), key=lambda x: x[0])) or "none")
    if g: print("  counts the tool tried:", g[0]["tested_counts"])
print(f"\nwrote {a.out} ({len(rows_out)} reads)")

if a.plot:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    os.makedirs("read_plots", exist_ok=True)
    for (s, l), g in groups.items():
        fig, ax = plt.subplots(figsize=(7, 3))
        for rtype, color, off in (("spanning", "#2a78d6", -0.2), ("flanking", "#eb6834", 0.2)):
            c = Counter(r["best_count"] for r in g if r["kept"] == "yes" and r["read_type"] == rtype)
            if c: ax.bar([k + off for k in c], list(c.values()), width=0.4, color=color,
                         label=rtype + (" (lower bound)" if rtype == "flanking" else " (exact)"))
        ax.set_xlabel("repeat count the read supports best"); ax.set_ylabel("reads")
        ax.set_title(f"{s.replace('ULTWGS_', '')}  {l}", loc="left", weight="bold")
        ax.spines[["top", "right"]].set_visible(False); ax.legend(frameon=False, fontsize=8)
        fig.savefig(f"read_plots/{s}_{l}.png", dpi=150, bbox_inches="tight"); plt.close(fig)
    print("wrote read_plots/")
