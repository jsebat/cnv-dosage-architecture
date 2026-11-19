# All-CNV version of fig_genecount_vs_summed_z.py (MS 2026-11-07): every CNV in twas_per_cnv_sign_test.csv (32), i.e. the 23 of the
# non-overlapping set plus the 9 overlapping sub-intervals (16p11.2 BP2-BP3, 1q21.1 TAR BP1-BP2 / BP2-BP3, 16p13.11 BP1-BP2 / BP2-BP3,
# 22q11.21 A-B / B-C / C-D / B-D). Not here: 15q13.3 CHRNA7 and 2q13 (not in the per-CNV TWAS table).
# No-effect loci (neither DEL nor DUP height P < 0.05) drawn as open markers. Overlapping CNVs share genes, so points are not independent.
# Gene count / summed z: contiguous run of TWAS30 genes spanning the CNV's genes, z = mean over panels (as panel A); checked against
# twas_class_strip_points_all.csv for height. Class = height class for both traits. Grey band: genome-wide null as before.
# Outputs: fig_genecount_vs_summed_z_all.{png,pdf}, genecount_vs_summed_z_all_points.csv, genecount_vs_summed_z_all_stats.csv
import os, numpy as np, pandas as pd, matplotlib
from scipy import stats
matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); CONC = os.path.dirname(HERE); ROOT = os.path.dirname(os.path.dirname(CONC))
chrnum = lambda c: int(str(c).replace("chr", "").replace("X", "23"))
genes = pd.read_csv(ROOT + "/pathway_analysis/Final/cnv_genes_annotated.csv")
H = pd.read_csv(CONC + "/twas_per_cnv_sign_test.csv"); H = H[H.trait == "height"].set_index("CNV_ID")
noeff = (H.DEL_p >= .05) & (H.DUP_p >= .05)
DROP = ["16p11.2_BP2_BP3", "1q21.1_TAR_BP1_BP2", "1q21.1_TAR_BP2_BP3", "16p13.11_BP1_BP2", "16p13.11_BP2_BP3",
        "22q11.21_A_B", "22q11.21_B_C", "22q11.21_C_D", "22q11.21_B-D"]
NMAX = 32; band, pts, rows = {}, [], []
for tr in ["Height", "BMI"]:
    tw = pd.read_csv(CONC + f"/twas30/{tr.lower()}.txt", sep=r"\s+")
    g = tw.groupby("GENE").agg(z=("TWAS.Z", "mean"), chr=("CHR", "first"), pos=("TxS", "first")).reset_index()
    g["c"] = g["chr"].map(chrnum); g = g.sort_values(["c", "pos"]).reset_index(drop=True); z = g.z.values; N = len(z); pos_of = dict(zip(g.GENE, g.index))
    cs = np.concatenate([[0], np.cumsum(np.concatenate([z, z]))]); ns = np.arange(1, NMAX + 1)
    Q = np.array([np.percentile(cs[np.arange(N) + n] - cs[np.arange(N)], [2.5, 25, 50, 75, 97.5]) for n in ns]); band[tr] = (ns, Q)
    for cid, sub in genes.groupby("CNV_ID"):
        if cid not in H.index: continue
        ix = sorted(pos_of[x] for x in sub.Gene if x in pos_of)
        if not ix: continue
        L = ix[-1] - ix[0] + 1
        pts.append(dict(trait=tr, CNV_ID=cid, cls=H.loc[cid, "cls"], n_genes_twas=L, z=cs[ix[0] + L] - cs[ix[0]], no_effect=bool(noeff[cid]), overlapping=cid in DROP))
P = pd.DataFrame(pts)
ref = pd.read_csv(CONC + "/twas_class_strip_points_all.csv").set_index("CNV_ID"); ph = P[P.trait == "Height"].set_index("CNV_ID")
assert np.allclose(ph.z, ref.summed_z.reindex(ph.index)) and (ph.n_genes_twas == ref.n_genes.reindex(ph.index)).all(), "height does not match twas_class_strip_points_all.csv"
P.to_csv(HERE + "/genecount_vs_summed_z_all_points.csv", index=False)
for tr in ["Height", "BMI"]:
    for grp in ["mirror", "asymmetric"]:
        for lab, s in [("all", P[(P.trait == tr) & (P.cls == grp)]), ("effect loci only", P[(P.trait == tr) & (P.cls == grp) & ~P.no_effect])]:
            r, p = stats.spearmanr(s.n_genes_twas, s.z)
            rows.append(dict(trait=tr, cls=grp, set=lab, n_cnv=len(s), n_no_effect=int(s.no_effect.sum()), median_genes=s.n_genes_twas.median(),
                             mean_z=s.z.mean(), spearman_rho_genes_vs_z=r, p_asymptotic=p))
R = pd.DataFrame(rows); R.to_csv(HERE + "/genecount_vs_summed_z_all_stats.csv", index=False); print(R.round(3).to_string(index=False))
print("no-effect loci:", ", ".join(sorted(noeff[noeff].index)))

INK = '#1a1a1a'; GREY = '#8a8a8a'; LIGHT = '#dcdcdc'; RED = '#b03a2e'; BLUE = '#1f5c8b'
plt.rcParams.update({"font.family": "DejaVu Sans", "axes.linewidth": .6, "xtick.major.width": .6, "ytick.major.width": .6})
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.5), gridspec_kw=dict(wspace=.22)); fig.subplots_adjust(left=.08, right=.98, top=.87, bottom=.23)
for ax, tr in zip(axs, ["Height", "BMI"]):
    ns, Q = band[tr]
    ax.fill_between(ns, Q[:, 0], Q[:, 4], color=LIGHT, alpha=.6, lw=0, zorder=1); ax.fill_between(ns, Q[:, 1], Q[:, 3], color=LIGHT, lw=0, zorder=1)
    ax.plot(ns, Q[:, 2], color=GREY, lw=.6, zorder=2); ax.axhline(0, color=GREY, lw=.4, ls=":", zorder=1)
    p = P[P.trait == tr]
    for k, (grp, c, mk) in enumerate([("mirror", BLUE, "o"), ("asymmetric", RED, "s")]):
        s = p[p.cls == grp]; e, ne = s[~s.no_effect], s[s.no_effect]
        ax.scatter(e.n_genes_twas, e.z, s=22, color=c, marker=mk, edgecolor="white", lw=.5, zorder=4)
        ax.scatter(ne.n_genes_twas, ne.z, s=22, facecolor="white", edgecolor=c, marker=mk, lw=.9, zorder=4)
        st = R[(R.trait == tr) & (R.cls == grp) & (R.set == "all")].iloc[0]; st2 = R[(R.trait == tr) & (R.cls == grp) & (R.set == "effect loci only")].iloc[0]
        ax.text(.02, .97 - .1 * k, f"{grp}: n = {st.n_cnv} ({st.n_no_effect} no-effect), ρ = {st.spearman_rho_genes_vs_z:+.2f} (P = {st.p_asymptotic:.2g})\n"
                f"   effect loci only: n = {st2.n_cnv}, ρ = {st2.spearman_rho_genes_vs_z:+.2f} (P = {st2.p_asymptotic:.2g})",
                transform=ax.transAxes, size=5.2, color=c, va="top", linespacing=1.3)
    for _, r in p[(p.z.abs() >= 4) | (p.n_genes_twas >= 10)].iterrows():
        off = {("BMI", "2q11.2"): (-5, 5), ("BMI", "10q11.21_11.23"): (-5, -5), ("Height", "22q11.21_A-D"): (-4, 7), ("BMI", "22q11.21_A-D"): (-4, -7), ("BMI", "16p13.11_BP1-BP3"): (4, 7), ("BMI", "22q11.21_B-D"): (2, -7)}.get((tr, r.CNV_ID), (4, 0))
        ax.annotate(r.CNV_ID.replace("_", " "), (r.n_genes_twas, r.z), xytext=off, textcoords="offset points", size=4.6, color=INK, va="center", ha="right" if off[0] < 0 else "left")
    ax.set_xlim(0, NMAX + 1); ax.set_ylim(*{"Height": (-15, 13), "BMI": (-15, 11)}[tr]); ax.set_title(tr, size=8, weight="bold"); ax.set_xlabel("TWAS-modelled genes in CNV", size=7)
    ax.tick_params(labelsize=6, length=2); ax.spines[["top", "right"]].set_visible(False)
axs[0].set_ylabel("summed TWAS z", size=7)
fig.suptitle("Gene count vs summed TWAS z per CNV, all CNVs (incl. overlapping and no-effect loci)", size=8.5, weight="bold", y=.97)
import textwrap
fig.text(.08, .01, textwrap.fill("Open markers: no-effect loci (neither DEL nor DUP height P < 0.05). All 32 CNVs with TWAS genes, including the 9 overlapping sub-intervals, so points share genes and "
         "are not independent. Grey: genome-wide null for n consecutive TWAS genes (light 95%, dark 50%, line median). Class = height class for both traits. "
         "Labelled: |z| ≥ 4 or ≥ 10 genes. ρ: Spearman (descriptive). No-effect loci: " + ", ".join(sorted(x.replace("_", " ") for x in noeff[noeff].index)) + ".", 185), size=5, color=GREY, va="bottom")
fig.savefig(HERE + "/fig_genecount_vs_summed_z_all.png", dpi=300); fig.savefig(HERE + "/fig_genecount_vs_summed_z_all.pdf")
