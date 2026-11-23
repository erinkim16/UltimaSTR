"""Schematic slide figures for the STR talk (no participant data).
Writes figures/talk_*.png: STR schematic, pipeline, read types, background-noise table."""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle

OUT = sys.argv[1] if len(sys.argv) > 1 else "figures"
os.makedirs(OUT, exist_ok=True)
INK, MUTED, FLANK = "#222222", "#666666", "#d4d4d4"
UNIT_A, UNIT_B, ORANGE = "#2a78d6", "#8fb8ec", "#eb6834"
plt.rcParams.update({"font.size": 12, "font.family": "DejaVu Sans"})

def units(ax, x, y, n, h=0.5, w=0.32):
    for i in range(n):
        ax.add_patch(Rectangle((x + i * w, y), w * 0.96, h, color=UNIT_A if i % 2 == 0 else UNIT_B, lw=0))
    return x + n * w

# 1. STR schematic: normal vs expanded
fig, ax = plt.subplots(figsize=(11, 3.6))
rows = [("Normal (HTT: 17 CAG)", 17, "healthy"), ("Expanded (HTT: 40+ CAG)", 42, "disease range")]
for k, (lab, n, note) in enumerate(rows):
    y = 1.6 - k * 1.2
    ax.add_patch(Rectangle((0, y), 3, 0.5, color=FLANK, lw=0))
    end = units(ax, 3, y, n)
    ax.add_patch(Rectangle((end, y), 3, 0.5, color=FLANK, lw=0))
    ax.text(-0.3, y + 0.25, lab, ha="right", va="center", color=INK)
    ax.text(end + 3.3, y + 0.25, note, ha="left", va="center", color=MUTED)
ax.text(3, 2.35, "CAG CAG CAG ...  (one block = one repeat unit)", color=UNIT_A, fontsize=11)
ax.text(1.5, -0.25, "flanking DNA", ha="center", color=MUTED, fontsize=10)
ax.set_xlim(-9, 26); ax.set_ylim(-0.5, 2.7); ax.axis("off")
ax.set_title("A short tandem repeat (STR): a short motif repeated in a row", loc="left", color=INK)
fig.savefig(f"{OUT}/talk_str_schematic.png", dpi=200, bbox_inches="tight"); plt.close(fig)

# 2. Pipeline flow
steps = [("Ultima WGS\nCRAMs", "40 people"), ("Catalog", "86 loci, 96 sites\n(Broad 83 + John 3)\nmale / female"),
         ("str_genotyper", "WDL + miniwdl\nApptainer, Trillium"), ("Merge calls", "3,830 calls"),
         ("Filters + tiers", "cohort filter,\nread support"), ("Read checks", "read_summary.py\nJBrowse")]
fig, ax = plt.subplots(figsize=(13, 3))
W, G = 1.75, 0.45
for i, (t, sub) in enumerate(steps):
    x = i * (W + G)
    ax.add_patch(FancyBboxPatch((x, 0), W, 1.2, boxstyle="round,pad=0.02,rounding_size=0.12",
                                fc="#eef4fc" if i < 3 else "#fdf0e9", ec=UNIT_A if i < 3 else ORANGE, lw=1.5))
    ax.text(x + W / 2, 0.85, t, ha="center", va="center", fontweight="bold", color=INK, fontsize=11)
    ax.text(x + W / 2, 0.35, sub, ha="center", va="center", color=MUTED, fontsize=9)
    if i < len(steps) - 1:
        ax.annotate("", xy=(x + W + G - 0.04, 0.6), xytext=(x + W + 0.04, 0.6),
                    arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.5))
ax.text(0, -0.35, "Running the tool (ported to Trillium)", color=UNIT_A, fontsize=10)
ax.text(3 * (W + G), -0.35, "Separating signal from noise (my analysis)", color=ORANGE, fontsize=10)
ax.set_xlim(-0.1, len(steps) * (W + G)); ax.set_ylim(-0.6, 1.4); ax.axis("off")
fig.savefig(f"{OUT}/talk_pipeline.png", dpi=200, bbox_inches="tight"); plt.close(fig)

# 3. Read types
fig, axs = plt.subplots(1, 2, figsize=(13, 4.2), sharey=True)
for ax, (title, n) in zip(axs, [("Short allele: reads can cross the repeat", 12), ("Long allele: the repeat is longer than a read", 34)]):
    ax.add_patch(Rectangle((0, 0), 10, 0.35, color=FLANK, lw=0))
    end = units(ax, 10, 0, n, h=0.35, w=0.5)
    ax.add_patch(Rectangle((end, 0), 10, 0.35, color=FLANK, lw=0))
    ax.text(10, -0.35, "repeat", color=UNIT_A, fontsize=10)
    for xe in (10, end):
        ax.axvline(xe, ymin=0.05, ymax=0.85, color=UNIT_A, lw=1, ls=":")
    L = 15  # read length (drawing units)
    if n == 12:
        reads = [(4, "spanning", INK), (7, "spanning", INK), (2, "spanning", INK)]
        reads = [(s, lab, c) for s, lab, c in reads]
        ys = [0.9, 1.4, 1.9]
        for (s, lab, c), y in zip(reads, ys):
            ax.annotate("", xy=(s + L + 4, y), xytext=(s, y), arrowprops=dict(arrowstyle="-|>", color=c, lw=2.5))
        ax.text(0, 2.35, "Spanning reads: anchored on BOTH sides → exact count", color=INK, fontsize=11)
    else:
        for s, y in [(2, 0.9), (5, 1.4)]:
            ax.annotate("", xy=(s + L, y), xytext=(s, y), arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=2.5))
        for s, y in [(end + 7, 1.9)]:
            ax.annotate("", xy=(s - L, y), xytext=(s, y), arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=2.5))
        ax.annotate("", xy=(22, 2.3), xytext=(14, 2.3), arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=2.5, ls="--"))
        ax.text(0, 2.75, "Flanking reads: anchored on ONE side → only a minimum (≥N)", color=ORANGE, fontsize=11)
        ax.text(23, 2.25, "in-repeat read (no anchor)", color=MUTED, fontsize=9, va="center")
    ax.set_title(title, loc="left", fontsize=12, color=INK)
    ax.set_xlim(-1, 48); ax.set_ylim(-0.6, 3.0); ax.axis("off")
fig.savefig(f"{OUT}/talk_read_types.png", dpi=200, bbox_inches="tight"); plt.close(fig)

# 4. Background-noise table (slide 9; cohort counts only, no IDs)
cols = ["Site", "People with stray\nlong reads", "The flagged call rests on", "Verdict"]
rows = [["DIP2B_CAC", "32 / 40", "4 flanking reads (others have up to 8)", "Noise"],
        ["FOXL2", "9 / 40", "1–2 flanking reads (others have up to 3)", "Noise"],
        ["HTT", "4 / 40", "1 spanning + 4 flanking, all agree on 29", "Supported"]]
fig, ax = plt.subplots(figsize=(12, 2.6)); ax.axis("off")
t = ax.table(cellText=rows, colLabels=cols, loc="center", cellLoc="left", colLoc="left", colWidths=[0.14, 0.2, 0.5, 0.16])
t.auto_set_font_size(False); t.set_fontsize(11); t.scale(1, 2.6)
for (r, c), cell in t.get_celld().items():
    cell.set_edgecolor("#e0e0e0"); cell.set_linewidth(0.8)
    if r == 0:
        cell.set_text_props(weight="bold", color=INK); cell.set_facecolor("#f4f4f4")
    elif c == 3:
        cell.set_text_props(weight="bold", color=ORANGE if rows[r - 1][3] == "Noise" else UNIT_A)
ax.set_title("A long call only counts if the site is otherwise clean", loc="left", color=INK, fontsize=13)
fig.savefig(f"{OUT}/talk_noise_table.png", dpi=200, bbox_inches="tight"); plt.close(fig)
print("wrote", sorted(f for f in os.listdir(OUT) if f.startswith("talk_")))
