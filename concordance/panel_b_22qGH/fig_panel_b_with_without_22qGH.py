# Panel b (height; summed TWAS z per CNV vs its gene-count-matched circular null, mirror vs asymmetric) with and without
# 22q11.21 G-H (MS 2026-11-03: G-H is a weak mirror CNV, neither allele nominally significant: DEL +0.36 SD P 0.41, DUP -0.16 SD
# P 0.07). Recomputes twas_class_null_panel.py's strip-plot values (same gene list, same 10,000 circular draws, seed 2026, same 23
# non-overlapping CNVs) and checks them against twas_class_strip_points.csv / twas_class_strip_stats.csv.
# Outputs: fig_panel_b_with_without_22qGH.{png,pdf}, panel_b_with_without_22qGH_stats.csv
import os, numpy as np, pandas as pd, matplotlib
from scipy import stats
matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); CONC = os.path.dirname(HERE); ROOT = os.path.dirname(os.path.dirname(CONC))
genes = pd.read_csv(ROOT + "/pathway_analysis/Final/cnv_genes_annotated.csv")
S = pd.read_csv(CONC + "/twas_per_cnv_sign_test.csv"); S = S[S.trait == "height"].set_index("CNV_ID")
chrnum = lambda c: int(str(c).replace("chr", "").replace("X", "23"))
tw = pd.read_csv(CONC + "/twas30/height.txt", sep=r"\s+")
g = tw.groupby("GENE").agg(z=("TWAS.Z", "mean"), chr=("CHR", "first"), pos=("TxS", "first")).reset_index()
g["c"] = g["chr"].map(chrnum); g = g.sort_values(["c", "pos"]).reset_index(drop=True)
N = len(g); cs = np.concatenate([[0], np.cumsum(np.concatenate([g.z.values, g.z.values]))]); pos_of = dict(zip(g.GENE, g.index))
spans = {}
for cid, sub in genes.groupby("CNV_ID"):
    ix = sorted(pos_of[x] for x in sub.Gene if x in pos_of)
    if ix and cid in S.index and not np.isnan(S.loc[cid, "DEL_effect"]) and not np.isnan(S.loc[cid, "DUP_effect"]): spans[cid] = (ix[0], ix[-1] - ix[0] + 1)
cids = list(spans); order = sorted(cids, key=lambda c: spans[c][0]); clusters, cur, end = [], [], -1
for c in order:
    s0, L = spans[c]
    if cur and s0 <= end: cur.append(c); end = max(end, s0 + L - 1)
    else:
        if cur: clusters.append(cur)
        cur, end = [c], s0 + L - 1
clusters.append(cur)
anchor = {c: min(spans[x][0] for x in cl) for cl in clusters for c in cl}; cl_id = {c: i for i, cl in enumerate(clusters) for c in cl}
rel = np.array([spans[c][0] - anchor[c] for c in cids]); Ls = np.array([spans[c][1] for c in cids]); ci = np.array([cl_id[c] for c in cids]); anc = np.array([anchor[c] for c in cids])
x = np.array([cs[spans[c][0] + Ls[i]] - cs[spans[c][0]] for i, c in enumerate(cids)]); cls = np.array([S.loc[c, "cls"] for c in cids])
rng = np.random.default_rng(2026); K = rng.integers(0, N, size=(10000, len(clusters)))
st = (anc[None, :] + K[:, ci] + rel[None, :]) % N; Xp = cs[st + Ls[None, :]] - cs[st]
mu, sd = Xp.mean(0), Xp.std(0); zs = (x - mu) / sd; Zp = (Xp - mu[None, :]) / sd[None, :]
C = np.array(cids); keep = np.isin(C, pd.read_csv(CONC + "/twas_dosage_panel_nonoverlap.csv").CNV_ID)

# check against the published panel b values
ref = pd.read_csv(CONC + "/twas_class_strip_points.csv").set_index("CNV_ID")
assert np.allclose(pd.Series(zs, index=C)[ref.index], ref.size_matched_z), "does not reproduce twas_class_strip_points.csv"

p2 = lambda o, nl: min(1, 2 * min(1 + np.sum(nl <= o), 1 + np.sum(nl >= o)) / (len(nl) + 1))
rows, sets = [], {"with 22q11.21 G-H": keep, "without 22q11.21 G-H": keep & (C != "22q11.21_G_H")}
for name, k in sets.items():
    for grp in ["mirror", "asymmetric"]:
        m = k & (cls == grp)
        rows.append(dict(cnv_set=name, group=grp, n_cnv=int(m.sum()), median=float(np.median(zs[m])), mean=zs[m].mean(), p_two_sided=p2(zs[m].mean(), Zp[:, m].mean(1))))
    ma, mm = k & (cls == "asymmetric"), k & (cls == "mirror"); o = zs[ma].mean() - zs[mm].mean()
    mw = stats.mannwhitneyu(zs[ma], zs[mm]).pvalue
    rows.append(dict(cnv_set=name, group="asymmetric minus mirror", mean=o, p_two_sided=p2(o, Zp[:, ma].mean(1) - Zp[:, mm].mean(1)), p_mannwhitney=mw))
R = pd.DataFrame(rows); R.to_csv(HERE + "/panel_b_with_without_22qGH_stats.csv", index=False)
chk = pd.read_csv(CONC + "/twas_class_strip_stats.csv").set_index("group"); r0 = R[R.cnv_set == "with 22q11.21 G-H"].set_index("group")
assert np.allclose(r0["mean"], chk.loc[r0.index, "mean"]), "class means do not reproduce twas_class_strip_stats.csv"
print(R.round(4).to_string(index=False))

# ---- figure (Fig. 3 style, as fig_twas_fig3style.py panel b)
INK = '#1a1a1a'; GREY = '#8a8a8a'; LIGHT = '#dcdcdc'; RED = '#b03a2e'; BLUE = '#1f5c8b'
plt.rcParams.update({"font.family": "DejaVu Sans", "axes.linewidth": .6, "xtick.major.width": .6, "ytick.major.width": .6})
fp = lambda p: f"{p:.2g}"
fig, axs = plt.subplots(1, 2, figsize=(4.4, 2.7), sharey=True, gridspec_kw=dict(wspace=.12))
fig.subplots_adjust(left=.15, right=.98, top=.8, bottom=.12)
for j, (bx, (name, k)) in enumerate(zip(axs, sets.items())):
    st_ = R[R.cnv_set == name].set_index("group"); rng = np.random.default_rng(1)
    for i, (grp, c, mk) in enumerate([("mirror", BLUE, "o"), ("asymmetric", RED, "s")]):
        m = k & (cls == grp); v = zs[m]; jit = rng.uniform(-.13, .13, len(v))
        gh = C[m] == "22q11.21_G_H"
        bx.scatter(i + jit[~gh], v[~gh], s=12, color=c, marker=mk, edgecolor="white", lw=.4, zorder=4)
        if gh.any():
            bx.scatter(i + jit[gh], v[gh], s=16, marker=mk, facecolor="white", edgecolor=c, lw=.9, zorder=4)
            bx.annotate("22q11.21 G–H", (i + jit[gh][0], v[gh][0]), xytext=(6, 0), textcoords="offset points", size=5, color=c, va="center")
        bx.plot([i - .22, i + .22], [np.median(v)] * 2, color=INK, lw=2, zorder=5, solid_capstyle="butt")
        bx.text(i, 1.0, f"n = {int(m.sum())}\nP = {fp(st_.loc[grp, 'p_two_sided'])}", size=5.2, ha="center", va="bottom", color=INK, linespacing=1.2)
    bx.axhline(0, color=LIGHT, lw=.6, zorder=1)
    bx.text(.5, 1.5, f"asy vs mir: P = {fp(st_.loc['asymmetric minus mirror', 'p_two_sided'])}", size=5.2, ha="center", va="bottom", color=GREY)
    bx.set_xticks([0, 1]); bx.set_xticklabels(["mir", "asy"], size=6.4); bx.set_xlim(-.6, 1.6); bx.set_ylim(-1.9, 1.75)
    bx.tick_params(labelsize=5.6, length=2, pad=1.5); bx.spines[["top", "right"]].set_visible(False)
    bx.set_title(name[0].upper() + name[1:].replace("G-H", "G\u2013H"), size=7, weight="bold", pad=4)
axs[0].set_ylabel("summed TWAS z vs gene-count-\nmatched null (SD)", size=6.4)
fig.suptitle("CNV genes vs random windows with equal gene count", size=7, weight="bold", y=.97)
fig.text(.15, .005, "Bar, median. P: two-sided, class mean vs circular null (10,000); asy vs mir: difference in class means.", size=4.6, color=GREY)
fig.savefig(HERE + "/fig_panel_b_with_without_22qGH.png", dpi=400); fig.savefig(HERE + "/fig_panel_b_with_without_22qGH.pdf")
print("written")
