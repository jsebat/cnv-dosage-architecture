# Schematic (no data) of the expected relationship between a region's summed height TWAS z and CNV height effects:
# z > 0 (higher expression -> taller): DEL shorter, DUP taller; z < 0 (higher expression -> shorter): DEL taller, DUP shorter.
# (MS 2026-11-08) Outputs: fig_twas_sign_schematic.{png,pdf}
import os, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
INK = '#1a1a1a'; LIGHT = '#dcdcdc'; UP = '#b03a2e'; DN = '#1f5c8b'; CDEL = '#d17a00'; CDUP = '#6a3d9a'
plt.rcParams.update({"font.family": "DejaVu Sans", "axes.linewidth": .6})
fig, axs = plt.subplots(1, 2, figsize=(6, 3), sharey=True); fig.subplots_adjust(left=.1, right=.98, top=.8, bottom=.16, wspace=.15)
for ax, s, col, title, sub in [(axs[0], 1, UP, "TWAS z > 0", "higher expression → taller"),
                               (axs[1], -1, DN, "TWAS z < 0", "higher expression → shorter")]:
    ax.axhline(0, color=INK, lw=.6, zorder=1)
    ax.plot([1, 2, 3], [-s, 0, s], color=col, lw=2, zorder=2)
    ax.scatter([1], [-s], s=90, color=CDEL, edgecolor="white", lw=1, zorder=3)
    ax.scatter([3], [s], s=90, color=CDUP, marker="s", edgecolor="white", lw=1, zorder=3)
    ax.scatter([2], [0], s=40, color=INK, zorder=3)
    ax.annotate("DEL " + ("shorter" if s > 0 else "taller"), (1, -s), xytext=(1.12, -s), va="center", size=8, color=INK)
    ax.annotate("DUP " + ("taller" if s > 0 else "shorter"), (3, s), xytext=(2.88, s), va="center", ha="right", size=8, color=INK)
    ax.set_title(f"{title}\n", size=10, weight="bold", color=col, pad=2)
    ax.text(.5, 1.02, sub, transform=ax.transAxes, ha="center", va="bottom", size=8, color=INK)
    ax.set_xticks([1, 2, 3]); ax.set_xticklabels(["Deletion\n(1 copy)", "Typical\n(2 copies)", "Duplication\n(3 copies)"], size=7.5)
    ax.set_xlim(.6, 3.4); ax.set_ylim(-1.5, 1.5); ax.set_yticks([]); ax.tick_params(length=0)
    ax.spines[["top", "right"]].set_visible(False)
axs[0].set_ylabel("Height effect\n← shorter     taller →", size=8)
fig.savefig(HERE + "/fig_twas_sign_schematic.png", dpi=300); fig.savefig(HERE + "/fig_twas_sign_schematic.pdf")
