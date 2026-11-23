"""Read-level STR pileup (REViewer-style) from a slice BAM made by slice_str.sh.
Each row = one real read. Flank bases grey, each repeat unit coloured (alternating shades),
other bases inside the repeat dark grey. Reads are classed as spanning or flanking
(anchored left / right) and labelled with their unit count.

Usage (Trillium: module load samtools scipy-stack):
  python make_read_pileup.py --bam slices/ULTWGS_0002_HTT.bam --bed slices/HTT.bed --motif CAG \
      --title "HTT, person A" --out figures/pileup_HTT_A.png
Contains participant reads: keep outputs in figures/ (gitignored), lab-internal only.
"""
import argparse, re, subprocess
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

ap = argparse.ArgumentParser()
ap.add_argument("--bam"); ap.add_argument("--sam", help="text SAM instead of --bam (testing)")
ap.add_argument("--bed", required=True); ap.add_argument("--motif", required=True)
ap.add_argument("--flank", type=int, default=10, help="aligned bases needed past the repeat (tool default 10)")
ap.add_argument("--show-flank", type=int, default=40); ap.add_argument("--max-rows", type=int, default=25)
ap.add_argument("--title", default=""); ap.add_argument("--out", required=True)
a = ap.parse_args()

chrom, start, end = open(a.bed).read().split()[:3]
start, end = int(start), int(end)          # 0-based half-open, like the catalog
K = len(a.motif); motifs = {a.motif}

if a.sam:
    lines = [l for l in open(a.sam) if not l.startswith("@")]
else:
    lines = subprocess.run(["samtools", "view", "-F", "0xF04", a.bam, f"{chrom}:{start - a.flank}-{end + a.flank}"],
                           capture_output=True, text=True, check=True).stdout.splitlines()

def walk(pos0, cigar):
    """map ref position -> query index for aligned (M/=/X) bases; also ref end"""
    m, q, r = {}, 0, pos0
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        n = int(n)
        if op in "M=X":
            for i in range(n): m[r + i] = q + i
            q += n; r += n
        elif op in "IS": q += n
        elif op in "DN": r += n
    return m, r

def unit_mask(seq):
    """per base: 2/3 = alternating repeat units, 4 = other base"""
    out, i, u = [], 0, 0
    while i < len(seq):
        if seq[i:i + K] in motifs and len(seq[i:i + K]) == K:
            out += [2 + u % 2] * K; i += K; u += 1
        else:
            out.append(4); i += 1
    return out, u

reads = {"spanning": [], "left": [], "right": []}
for l in lines:
    f = l.split("\t")
    pos0, cigar, seq = int(f[3]) - 1, f[5], f[9]
    if cigar == "*" or seq == "*": continue
    m, rend = walk(pos0, cigar)
    lead = int(re.match(r"^(\d+)S", cigar).group(1)) if re.match(r"^\d+S", cigar) else 0
    left_ok = pos0 <= start - a.flank and (start - 1) in m
    right_ok = rend >= end + a.flank and end in m
    if left_ok and right_ok:
        qa, qb = m[start - 1] + 1, m[end]
        reads["spanning"].append((seq, qa, qb))
    elif left_ok and rend > start - a.flank:
        reads["left"].append((seq, m[start - 1] + 1, len(seq)))
    elif right_ok and (pos0 < end + a.flank or lead > 0):
        reads["right"].append((seq, 0, m[end]))

rows = []   # (label, codes, offset) offset = x of first code relative to repeat edge
for kind in ("spanning", "left"):
    for seq, qa, qb in reads[kind]:
        lf = seq[max(0, qa - a.show_flank):qa]
        rep, n = unit_mask(seq[qa:qb])
        rf = seq[qb:qb + a.show_flank] if kind == "spanning" else ""
        codes = [1] * len(lf) + rep + [1] * len(rf)
        lab = f"spanning: {n}" if kind == "spanning" else f"flanking: ≥{n}"
        rows.append((kind, n, lab, codes, -len(lf)))
right_rows = []
for seq, qa, qb in reads["right"]:
    rep, n = unit_mask(seq[:qb])
    rf = seq[qb:qb + a.show_flank]
    codes = rep + [1] * len(rf)
    right_rows.append(("right", n, f"flanking: ≥{n}", codes, -len(rep)))

rows.sort(key=lambda r: (r[0] != "spanning", r[1]))
rows, right_rows = rows[:a.max_rows], sorted(right_rows, key=lambda r: r[1])[:a.max_rows // 2]
cmap = ListedColormap(["white", "#d4d4d4", "#2a78d6", "#8fb8ec", "#555555"])

def draw(ax, rr, edge_label):
    if not rr:
        ax.text(0.5, 0.5, "no reads", transform=ax.transAxes, ha="center"); ax.axis("off"); return
    lo = min(r[4] for r in rr); hi = max(r[4] + len(r[3]) for r in rr)
    M = np.zeros((len(rr), hi - lo), dtype=int)
    for i, r in enumerate(rr):
        M[i, r[4] - lo:r[4] - lo + len(r[3])] = r[3]
    ax.imshow(M, cmap=cmap, vmin=0, vmax=4, aspect="auto", interpolation="nearest",
              extent=(lo, hi, len(rr), 0))
    ax.hlines(range(1, len(rr)), lo, hi, color="white", lw=2)   # gap between reads
    ax.axvline(0, color="#2a78d6", ls=":", lw=1)
    for i, r in enumerate(rr):
        ax.text(hi + 2, i + 0.5, r[2], va="center", fontsize=8, color="#222222")
    ax.set_yticks([]); ax.set_xlabel(f"bases from {edge_label}")
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)

fig, axs = plt.subplots(2, 1, figsize=(12, 0.22 * (len(rows) + len(right_rows)) + 3),
                        gridspec_kw={"height_ratios": [max(len(rows), 1), max(len(right_rows), 1)]})
draw(axs[0], rows, "repeat start (reads anchored on the left flank)")
draw(axs[1], right_rows, "repeat end (reads anchored on the right flank)")
axs[0].set_title(f"{a.title}  ({len(reads['spanning'])} spanning, "
                 f"{len(reads['left']) + len(reads['right'])} flanking reads; {a.motif} units coloured)", loc="left")
from matplotlib.patches import Patch
fig.legend(handles=[Patch(color="#d4d4d4", label="flanking DNA (outside the repeat)"),
                    Patch(color="#2a78d6", label=f"one {a.motif} unit"), Patch(color="#8fb8ec", label=f"next {a.motif} unit"),
                    Patch(color="#555555", label="base that breaks the pattern (error or interruption)")],
           loc="lower center", ncol=4, frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.04))
import os; os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
fig.tight_layout(); fig.savefig(a.out, dpi=200, bbox_inches="tight")
print(f"{a.out}: {len(reads['spanning'])} spanning, {len(reads['left'])} left-flanking, {len(reads['right'])} right-flanking")
