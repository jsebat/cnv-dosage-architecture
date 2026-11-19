# Summed TWAS z (height) by three CNV classes: mirror, asymmetric DEL-driven, asymmetric DUP-driven (MS 2026-11-06).
# Asymmetric CNVs (Fig. 3 class) split by the allele with the larger |height effect| (twas_per_cnv_sign_test.csv):
#   DEL-driven: 7q11.23 WBS, 22q11.21 A-D, 22q11.21 A-B (added by MS; nested in A-D), 15q13.1-13.2 BP3-BP4
#   DUP-driven: 17q12 (set by hand, MS: DEL -0.22 P 0.10, DUP -0.17 P 3e-4), 10q11.21-11.23, 15q13.3 BP4-BP5, 16p11.2 BP1-BP3, 17p12, 2q11.2, 7q11.23 distal
# Primary non-overlapping set without Williams (drawn as an excluded outlier) plus 22q11.21 A-B. Left: summed z. Right: size-matched z (vs gene-count-matched circular
# null, as panel b). P per class: class mean vs circular null (10,000 draws, seed 2026), two-sided; reuses
# fig_panel_b_with_without_22qGH.py's computation. Expected sign if concordant: dosage effect (DUP - DEL)/2, i.e. negative for
# DUP-driven (all DUP effects < 0) and positive for DEL-driven where the DEL lowers height.
# Outputs: fig_three_class_box.{png,pdf}, three_class_box_stats.csv
import os, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
src = open(HERE + "/fig_panel_b_with_without_22qGH.py").read().split("# check against")[0]; ns = {"__file__": HERE + "/x.py"}; exec(src, ns)
zs, Zp, C, cls, keep, xsum, S = ns["zs"], ns["Zp"], ns["C"], ns["cls"], ns["keep"], ns["x"], ns["S"]
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
p2 = lambda o, nl: min(1, 2 * min(1 + np.sum(nl <= o), 1 + np.sum(nl >= o)) / (len(nl) + 1))
st = ns["st"]; Ls = ns["Ls"]; cs = ns["cs"]; Xraw = cs[st + Ls[None, :]] - cs[st]
dom = np.array(["DEL" if abs(S.loc[c, "DEL_effect"]) > abs(S.loc[c, "DUP_effect"]) else "DUP" for c in C])
dom[C == "17q12"] = "DUP"   # MS 2026-11-06: 17q12 is DUP-driven (only the DUP effect is significant)
keep = keep | (C == "22q11.21_A_B")   # MS 2026-11-06: add 22q11.21 A-B to the DEL-driven class (nested in A-D; shares its first 18 TWAS genes)
WBS = C == "7q11.23_WBS"; keep = keep & ~WBS   # MS 2026-11-06: Williams excluded from boxes and stats, drawn as an outlier (open marker)
grp = np.where(cls == "mirror", "mirror", np.where(dom == "DEL", "asym DEL", "asym DUP"))
ORDER = ["mirror", "asym DEL", "asym DUP"]; rows = []
for g in ORDER:
    m = keep & (grp == g)
    rows.append(dict(group=g, n=int(m.sum()), cnvs="; ".join(C[m]), mean_summed_z=xsum[m].mean(), median_summed_z=np.median(xsum[m]),
                     p_summed=p2(xsum[m].mean(), Xraw[:, m].mean(1)), mean_size_matched_z=zs[m].mean(), p_size_matched=p2(zs[m].mean(), Zp[:, m].mean(1))))
R = pd.DataFrame(rows); R.to_csv(HERE + "/three_class_box_stats.csv", index=False); print(R.drop(columns="cnvs").round(3).to_string(index=False))

INK = '#1a1a1a'; GREY = '#8a8a8a'; LIGHT = '#dcdcdc'; COL = {"mirror": '#1f5c8b', "asym DEL": '#b03a2e', "asym DUP": '#d98a3a'}
LAB = {"mirror": "mirror", "asym DEL": "asymmetric\nDEL-driven", "asym DUP": "asymmetric\nDUP-driven"}
plt.rcParams.update({"font.family": "DejaVu Sans", "axes.linewidth": .6, "xtick.major.width": .6, "ytick.major.width": .6})
SHOW = {"22q11.21_A-D": "22q11.21 A–D", "22q11.21_A_B": "22q11.21 A–B", "16p11.2_BP1-BP3": "16p11.2 BP1–BP3", "2q11.2": "2q11.2",
        "22q11.21_G_H": "22q11.21 G–H", "3q29": "3q29", "17q12": "17q12", "15q13.1-13.2_BP3_BP4": "15q13.1–13.2"}
def draw(ax, v, ylab, pk):
    rng = np.random.default_rng(3)
    data = [v[keep & (grp == g)] for g in ORDER]
    bp = ax.boxplot(data, positions=range(3), widths=.5, showfliers=False, patch_artist=True, medianprops=dict(color=INK, lw=1.4),
                    whiskerprops=dict(color=GREY, lw=.7), capprops=dict(color=GREY, lw=.7), boxprops=dict(lw=.7, edgecolor=GREY))
    for b, g in zip(bp["boxes"], ORDER): b.set_facecolor(COL[g]); b.set_alpha(.18)
    for i, g in enumerate(ORDER):
        m = keep & (grp == g); jit = rng.uniform(-.12, .12, m.sum())
        ax.scatter(i + jit, v[m], s=16, color=COL[g], edgecolor="white", lw=.4, zorder=4)
        if g == "asym DEL":
            ax.scatter([i], v[WBS], s=18, facecolor="white", edgecolor=COL[g], lw=.9, zorder=4)
            ax.annotate("Williams (excluded)", (i, v[WBS][0]), xytext=(5, 0), textcoords="offset points", size=4.8, color=GREY, va="center", zorder=6)
        for xj, yv, c in zip(i + jit, v[m], C[m]):
            if c in SHOW: ax.annotate(SHOW[c], (xj, yv), xytext=(5, 0), textcoords="offset points", size=4.8, color=INK, va="center", zorder=6, bbox=dict(fc="white", ec="none", pad=.4, alpha=.85))
        r = R.set_index("group").loc[g]
        ax.text(i, 1.0, f"n = {r.n}\nP = {r[pk]:.2g}", transform=ax.get_xaxis_transform(), size=5.4, ha="center", va="bottom", color=INK)
    ax.axhline(0, color=LIGHT, lw=.6, zorder=1); ax.set_xticks(range(3)); ax.set_xticklabels([LAB[g] for g in ORDER], size=6.2)
    ax.set_xlim(-.6, 2.75); ax.set_ylabel(ylab, size=6.6); ax.tick_params(labelsize=5.8, length=2); ax.spines[["top", "right"]].set_visible(False)
fig, axs = plt.subplots(1, 2, figsize=(6.6, 3.2), gridspec_kw=dict(wspace=.32)); fig.subplots_adjust(left=.09, right=.98, top=.84, bottom=.23)
for ax, (v, ylab, pk) in zip(axs, [(xsum, "summed TWAS z (height)", "p_summed"), (zs, "summed z vs gene-count-matched\nnull (SD)", "p_size_matched")]): draw(ax, v, ylab, pk)
fig.suptitle("Height TWAS z by CNV class", size=8, weight="bold", y=.98)
import textwrap
FOOT = ("Box: median and IQR; whiskers 1.5 IQR. P: class mean vs circular genomic null (10,000), two-sided. Asymmetric CNVs split by the allele with the larger "
         "height effect, except 17q12, assigned DUP-driven (only its DUP effect is significant: DUP −0.17, P 3e-4; DEL −0.22, P 0.10). Non-overlapping CNV set "
         "(Williams excluded from boxes and P, shown as an open circle) plus 22q11.21 A–B, which is nested in A–D (overlapping genes; the circular null moves them together).")
fig.text(.09, .01, textwrap.fill(FOOT, 175), size=4.8, color=GREY, va="bottom")
fig.savefig(HERE + "/fig_three_class_box.png", dpi=400); fig.savefig(HERE + "/fig_three_class_box.pdf")


# single-panel version: gene-count-corrected (size-matched) z only (MS 2026-11-06)
fig, ax = plt.subplots(figsize=(3.4, 3.3)); fig.subplots_adjust(left=.17, right=.97, top=.84, bottom=.27)
draw(ax, zs, "summed TWAS z vs gene-count-matched\nnull (SD)", "p_size_matched")
fig.suptitle("Height TWAS z by CNV class (gene-count corrected)", size=7.5, weight="bold", y=.98)
fig.text(.03, .01, textwrap.fill(FOOT.replace("P: class mean", "Each CNV: summed z standardised against random windows with the same number of TWAS genes. P: class mean"), 95), size=4.5, color=GREY, va="bottom")
fig.savefig(HERE + "/fig_three_class_box_sizematched.png", dpi=400); fig.savefig(HERE + "/fig_three_class_box_sizematched.pdf")