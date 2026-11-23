#!/usr/bin/env python3
"""Make slide figures from the Vincent_ASD STR results.

Inputs (same folder): vincent_str_calls.filtered.tsv, vincent_str_calls.tiered.tsv,
site_quality.tsv. Outputs PNGs to figures/. The figures show participant IDs, so
keep them inside the lab (figures/ is gitignored).

usage: python3 make_figures.py [--stem vincent_str_calls] [--anon]
  --anon  replace participant IDs with letters (A, B, C...) for wider audiences
"""
import argparse, csv, os, re
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

EXACT, LOWER = "#2a78d6", "#eb6834"          # validated pair (blue / orange)
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
CRIT, WARN = "#d03b3b", "#c98500"
plt.rcParams.update({
    "font.size": 11, "axes.edgecolor": "#c3c2b7", "axes.labelcolor": INK2,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
    "axes.spines.right": False, "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.dpi": 200, "savefig.bbox": "tight",
})

ap = argparse.ArgumentParser()
ap.add_argument("--stem", default="vincent_str_calls")
ap.add_argument("--anon", action="store_true")
a = ap.parse_args()
os.makedirs("figures", exist_ok=True)

def read(p): return list(csv.DictReader(open(p), delimiter="\t"))
rows = read(f"{a.stem}.filtered.tsv")
tiers = read(f"{a.stem}.tiered.tsv")
sites = read("site_quality.tsv")

ids = {}
def name(s):
    if not a.anon: return s.replace("ULTWGS_", "")
    return ids.setdefault(s, chr(ord("A") + len(ids)))

def alleles(g):
    """'(14,>=45)' -> [(14, True), (45, False)]  (value, is_exact)"""
    return [(int(n), not ge) for ge, n in re.findall(r"(>=)?(\d+)", g)]

def short(loc):
    """'RUNX2__6-45422750-45422801-GCN' -> 'RUNX2 (GCN)'; long motif names trimmed"""
    if "__" in loc:
        return f"{loc.split('__')[0]} ({loc.rsplit('-', 1)[-1]})"
    base, _, motif = loc.partition("_")
    return f"{base}_{motif[:6]}\u2026" if len(motif) > 8 else loc

def save(fig, f):
    fig.savefig(f"figures/{f}"); plt.close(fig); print("wrote figures/" + f)

# ---------- Figure 1: from all calls to strong leads ----------
n_people = len({r["sample"] for r in rows})
n_sites = len({r["locus"] for r in rows})
# "flagged" = above a threshold (INTERMEDIATE / PATHOGENIC_RANGE); FLANKING_ONLY alone is not a flag here
flagged = sum(1 for r in rows if "INTERMEDIATE" in r["flag"] or "PATHOGENIC" in r["flag"])
leads = sum(1 for r in rows if r["lead"] == "yes")
t = defaultdict(int)
for r in tiers: t[r["tier"]] += 1
steps = [
    (f"All calls\n{n_people} people x {n_sites} sites", len(rows), "#b7d3f6"),
    ("Flagged\nvs. published thresholds", flagged, "#86b6ef"),
    ("Leads\nflag at site with ≤50% flagged,\nor John's site in cohort top 5%", leads, "#3987e5"),
    ("Tier 1 + 2\nexact or disease-size", t["1"] + t["2"], "#1c5cab"),
    ("Tier 1\n≥3 spanning reads, ≥3 over normal", t["1"], "#0d366b"),
]
fig, ax = plt.subplots(figsize=(9, 4.4))
for i, (lab, n, c) in enumerate(steps):
    ax.barh(i, n, color=c, height=0.62)
    ax.text(n * 1.12, i, f"{n:,}", va="center", color=INK, fontsize=12, weight="bold")
ax.set_yticks(range(len(steps)), [s[0] for s in steps], color=INK2)
ax.invert_yaxis(); ax.set_xscale("log"); ax.set_xlim(1, len(rows) * 3)
ax.set_xlabel("number of calls (log scale)")
ax.xaxis.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)
ax.set_title("From every call to the strongest leads", loc="left", color=INK, weight="bold")
save(fig, "fig1_filter_funnel.png")

# ---------- Figure 2: which sites can be measured ----------
s2 = sorted(sites, key=lambda r: float(r["fully_spanned_pct"] or 0))
vals = [float(r["fully_spanned_pct"] or 0) for r in s2]
fig, ax = plt.subplots(figsize=(11, 4))
cols = [CRIT if v < 50 else (WARN if v < 80 else "#86b6ef") for v in vals]
ax.bar(range(len(vals)), vals, color=cols, width=0.8)
ax.axhline(80, color=INK2, lw=1, ls="--")
ax.text(len(vals) - 0.5, 76, "80% = 'well measured'", ha="right", va="top", color=INK2, fontsize=9)
for i, (r, v) in enumerate(zip(s2, vals)):
    if v < 80:
        ax.text(i, v + 2, short(r["locus"]), rotation=90, ha="center", va="bottom", fontsize=8, color=INK)
ax.set_xticks([]); ax.set_ylim(0, 125); ax.set_xlim(-1, len(vals))
ax.set_ylabel("% of people with every\nallele fully spanned")
good = sum(v >= 80 for v in vals)
ax.set_title(f"{good} of {len(vals)} sites are well measured; long, pure, GC-rich repeats are not",
             loc="left", color=INK, weight="bold")
ax.set_xlabel("96 catalog sites, sorted (red < 50%, amber 50–80%)")
save(fig, "fig2_site_measurability.png")

# ---------- shared: one dot-plot panel per site ----------
q = {r["locus"]: r for r in sites}
def call_alleles(r):
    """[(value, exact)] where exact = no '>=' AND that allele has SPANNING support"""
    sup = r["support"].split("/")
    out = []
    for k, (ge, n) in enumerate(re.findall(r"(>=)?(\d+)", r["genotype"])):
        out.append((int(n), not ge and k < len(sup) and sup[k] == "SPANNING"))
    return out

def cohort_cut(loc):
    v = sorted(x for r in rows if r["locus"] == loc for x, _ in alleles(r["genotype"]))
    return v[int(0.95 * (len(v) - 1))] if v else None

def panel(ax, loc, label_samples, show_cut=False, note=""):
    pts, nm, pm = [], "", ""
    for r in rows:
        if r["locus"] == loc:
            nm, pm = r["normal_max"], r["pathogenic_min"]
            pts += [(v, ex, r["sample"]) for v, ex in call_alleles(r)]
    pts.sort(); stack = defaultdict(int)
    for v, ex, s in pts:
        y = stack[v]; stack[v] += 1
        if ex: ax.scatter(v, y, s=30, color=EXACT, edgecolor="white", lw=0.8, zorder=3)
        else:  ax.scatter(v, y, s=30, facecolor="white", edgecolor=LOWER, lw=1.5, zorder=3)
    hi = max([p[0] for p in pts] or [1]); ymax = max(stack.values() or [1]) * 1.3 + 1
    for val, c, lab, yf in ((nm, INK2, "normal", 0.97), (pm, CRIT, "disease", 0.86)):
        if val and val.isdigit() and int(val) <= hi * 1.6:
            ax.axvline(int(val), color=c, lw=1, ls="--", zorder=1)
            ax.text(int(val), ymax * yf, " " + lab + (" \u2264" if lab == "normal" else " \u2265") + val,
                    color=c, fontsize=8, va="top")
    if show_cut and pts:
        c = cohort_cut(loc)
        ax.axvspan(c + 0.5, hi + 2, color="#f0efec", zorder=0)
        ax.text(c + 0.6, ymax * 0.97, " cohort top 5%", color=INK2, fontsize=8, va="top")
    done = set()
    for v, ex, s in sorted(pts, reverse=True):
        if s in label_samples and s not in done:
            done.add(s)
            ax.annotate(name(s), (v, stack_y(pts, v, s)), xytext=(0, 14 + 12 * (len(done) - 1)),
                        textcoords="offset points", ha="center", fontsize=8.5, color=INK, weight="bold",
                        arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
    ne = sum(p[1] for p in pts)
    sub = f"{ne} exact, {len(pts) - ne} not fully spanned"
    if loc in q: sub += f"; site fully spanned in {float(q[loc]['fully_spanned_pct']):.0f}%"
    off = [f"{lab} {sym}{val}" for val, lab, sym in ((nm, "normal", "\u2264"), (pm, "disease", "\u2265"))
           if val and val.isdigit() and int(val) > hi * 1.6]
    if off: note = (note + "; " if note else "") + "off-scale: " + ", ".join(off)
    ax.text(0, 1.0, sub + (f"\n{note}" if note else ""), transform=ax.transAxes, fontsize=8,
            color=INK2, va="bottom")
    ax.set_title(loc, loc="left", fontsize=11, color=INK, weight="bold", pad=24 if note else 14)
    ax.set_yticks([]); ax.spines["left"].set_visible(False)
    ax.set_ylim(-1, ymax); ax.set_xlabel("repeat count", fontsize=9)

def stack_y(pts, v, s):
    same = [p for p in pts if p[0] == v]
    return [p[2] for p in same].index(s)

def legend(fig, y=1.02):
    fig.legend(handles=[
        plt.Line2D([], [], marker="o", ls="", color=EXACT, label="exact count (spanning reads)"),
        plt.Line2D([], [], marker="o", ls="", markerfacecolor="white", markeredgecolor=LOWER,
                   markeredgewidth=1.5, label="not fully spanned (\u2265N or flanking estimate)")],
        loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, y))

# ---------- Figure 3: cohort distributions at key sites ----------
top12 = defaultdict(set)
for r in tiers:
    if r["tier"] in ("1", "2"): top12[r["locus"]].add(r["sample"])
fig, axes = plt.subplots(2, 3, figsize=(13, 7))
for ax, loc in zip(axes.flat, ["HTT", "FXN", "ZIC2", "DIP2B_CAC", "CA10", "FMR1"]):
    panel(ax, loc, top12[loc], show_cut=loc in ("DIP2B_CAC", "CA10"))
legend(fig)
fig.suptitle("Each dot is one allele from one person (labels = tier 1/2 leads)", x=0.01, ha="left",
             y=1.07, color=INK, weight="bold", fontsize=13)
fig.tight_layout()
save(fig, "fig3_key_sites.png")

# ---------- Figure 5: John's regions of interest ----------
JOHN = [("FMR1", "Fragile X; CGG"), ("AFF2", "FRAXE; CCG"), ("PTCHD1", "no thresholds"),
        ("CA10", "no thresholds"), ("TCF4", "FECD; CTG"), ("ATXN8OS", "SCA8; CTG"),
        ("HTT", "Huntington; CAG"), ("DIP2B_CAC", "no thresholds; CAC")]
nocall = {loc: n_people - len({r["sample"] for r in rows if r["locus"] == loc}) for loc, _ in JOHN}
fig, axes = plt.subplots(2, 4, figsize=(16, 7.4))
for ax, (loc, what) in zip(axes.flat, JOHN):
    note = what + (f"; {nocall[loc]} people with no call" if nocall[loc] else "")
    panel(ax, loc, top12[loc], show_cut=what.startswith("no thresholds"), note=note)
legend(fig, 1.03)
fig.suptitle("John's regions of interest in Vincent_ASD (labels = tier 1/2 leads)", x=0.01, ha="left",
             y=1.08, color=INK, weight="bold", fontsize=13)
fig.tight_layout()
save(fig, "fig5_john_sites.png")

# ---------- Figure 4: the tier 1 + 2 leads ----------
lead_rows = sorted([r for r in tiers if r["tier"] in ("1", "2")],
                   key=lambda r: (r["tier"], r["locus"]))
fig, ax = plt.subplots(figsize=(11, 0.42 * len(lead_rows) + 1.2))
ax.axis("off")
hdr = ["tier", "person", "site", "genotype", "spanning reads", "normal max", "disease ≥", "why"]
cells = [[r["tier"], name(r["sample"]), short(r["locus"]), r["genotype"], r["spanning_reads"],
          r["normal_max"], r["pathogenic_min"], r["tier_reason"]] for r in lead_rows]
tb = ax.table(cellText=cells, colLabels=hdr, loc="center", cellLoc="left", colLoc="left")
tb.auto_set_font_size(False); tb.set_fontsize(9); tb.auto_set_column_width(range(len(hdr)))
for (i, j), c in tb.get_celld().items():
    c.set_edgecolor(GRID); c.set_height(0.9 / (len(cells) + 1))
    if i == 0: c.set_text_props(weight="bold", color=INK2); c.set_facecolor("#f9f9f7")
    elif cells[i - 1][0] == "1": c.set_facecolor("#e7f0fc")
ax.set_title(f"Tier 1 (shaded) and tier 2 leads: {t['1']} + {t['2']}", loc="left",
             color=INK, weight="bold")
save(fig, "fig4_top_leads.png")
