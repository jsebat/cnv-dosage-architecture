# Shows that the sign of a CNV region's summed height TWAS z predicts the direction of its DEL and DUP height effects
# (z > 0: DEL lowers height, DUP raises it; z < 0: the reverse). Input: twas_dosage_panel.csv (summed_twas_z over the contiguous
# TWAS span, DEL / DUP height effects). Top: summed TWAS z per CNV, sorted. Bottom: DEL and DUP effect for the same CNV, joined by a line.
# * = asymmetric CNV. (MS 2026-11-08) Outputs: fig_twas_sign_del_dup.{png,pdf}
import os, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
d = pd.read_csv(HERE + "/twas_dosage_panel.csv").sort_values("summed_twas_z").reset_index(drop=True)
x = np.arange(len(d)); z = d.summed_twas_z.values
INK = '#1a1a1a'; LIGHT = '#dcdcdc'; UP = '#b03a2e'; DN = '#1f5c8b'; CDEL = '#d17a00'; CDUP = '#6a3d9a'
plt.rcParams.update({"font.family": "DejaVu Sans", "axes.linewidth": .6, "xtick.major.width": .6, "ytick.major.width": .6})
fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.5, 5), sharex=True, gridspec_kw=dict(height_ratios=[1, 1.6], hspace=.08))
fig.subplots_adjust(left=.1, right=.98, top=.92, bottom=.2)
split = (z < 0).sum() - .5
for a in (a1, a2):
    a.axvline(split, color=INK, lw=.6, ls=":"); a.spines[["top", "right"]].set_visible(False); a.tick_params(labelsize=6, length=2)
a1.bar(x, z, width=.7, color=np.where(z > 0, UP, DN), zorder=3); a1.axhline(0, color=INK, lw=.5)
a1.set_ylabel("Summed height\nTWAS z", size=7)
a1.text(split / 2, a1.get_ylim()[1] * .92, "TWAS z < 0", ha="center", va="top", size=7.5, weight="bold", color=DN)
a1.text((split + len(d)) / 2, a1.get_ylim()[1] * .92, "TWAS z > 0", ha="center", va="top", size=7.5, weight="bold", color=UP)
a2.axhline(0, color=INK, lw=.5)
for i, r in d.iterrows():
    a2.plot([i, i], [r.DEL, r.DUP], color=LIGHT if (r.DUP > r.DEL) != (r.summed_twas_z > 0) else '#9a9a9a', lw=1.2, zorder=2)
a2.scatter(x, d.DEL, s=22, color=CDEL, zorder=3, label="Deletion", edgecolor="white", lw=.6)
a2.scatter(x, d.DUP, s=22, color=CDUP, zorder=3, label="Duplication", marker="s", edgecolor="white", lw=.6)
a2.set_ylabel("Height effect (SD)", size=7); a2.legend(frameon=False, fontsize=6.5, loc="lower right", handletextpad=.2)
neg, pos = d[z < 0], d[z > 0]
k_neg = int((neg.DEL > neg.DUP).sum()); k_pos = int((pos.DUP > pos.DEL).sum())
a2.text(split / 2, a2.get_ylim()[1] * .95, f"DEL > DUP in {k_neg} of {len(neg)}", ha="center", va="top", size=7, color=INK)
a2.text((split + len(d)) / 2, a2.get_ylim()[1] * .95, f"DUP > DEL in {k_pos} of {len(pos)}", ha="center", va="top", size=7, color=INK)
a2.set_xticks(x); a2.set_xticklabels([c + ("*" if k == "asymmetric" else "") for c, k in zip(d.CNV_ID, d.cls)], rotation=45, ha="right", rotation_mode="anchor", size=5.5)
a2.set_xlim(-.6, len(d) - .4)
fig.suptitle("Regional height TWAS z vs. deletion and duplication height effects", size=9, weight="bold", y=.98)
print("z<0: DEL>DUP", k_neg, "/", len(neg), "| z>0: DUP>DEL", k_pos, "/", len(pos))
for cls in ["mirror", "asymmetric"]:
    s = d[d.cls == cls]; print(cls, "concordant", int(((s.DUP > s.DEL) == (s.summed_twas_z > 0)).sum()), "/", len(s))
fig.savefig(HERE + "/fig_twas_sign_del_dup.png", dpi=300); fig.savefig(HERE + "/fig_twas_sign_del_dup.pdf")
