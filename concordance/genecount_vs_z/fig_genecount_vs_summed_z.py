# Number of TWAS-modelled genes vs summed TWAS z per CNV, height and BMI, mirror vs asymmetric (MS 2026-11-06).
# Points: primary non-overlapping CNV set (height: twas_dosage_panel_nonoverlap.csv, n = 23; BMI: twas_dosage_bmi.csv, n = 22,
# no 7q11.23 distal). Gene count = contiguous run of TWAS30 genes spanning the CNV; z = mean over TWAS30 panels (as panel A).
# Class = HEIGHT class for both traits (BMI has no asymmetric class of its own; every CNV is mirror for BMI).
# Grey band: null for a window of n consecutive TWAS genes anywhere in the genome (all N circular windows; 2.5-97.5% and 25-75%).
# Outputs: fig_genecount_vs_summed_z.{png,pdf}, genecount_vs_summed_z_points.csv, genecount_vs_summed_z_stats.csv
import os, numpy as np, pandas as pd, matplotlib
from scipy import stats
matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); CONC = os.path.dirname(HERE)
chrnum = lambda c: int(str(c).replace("chr", "").replace("X", "23"))
h = pd.read_csv(CONC + "/twas_dosage_panel_nonoverlap.csv").rename(columns={"summed_twas_z": "z"})
b = pd.read_csv(CONC + "/twas_dosage_bmi.csv").rename(columns={"summed_bmi_twas_z": "z"}).merge(h[["CNV_ID", "cls"]], on="CNV_ID")
D = {"Height": h, "BMI": b}; NMAX = 32; band, pts, rows = {}, [], []
for tr, d in D.items():
    tw = pd.read_csv(CONC + f"/twas30/{tr.lower()}.txt", sep=r"\s+")
    g = tw.groupby("GENE").agg(z=("TWAS.Z", "mean"), chr=("CHR", "first"), pos=("TxS", "first")).reset_index()
    g["c"] = g["chr"].map(chrnum); z = g.sort_values(["c", "pos"]).z.values; N = len(z)
    cs = np.concatenate([[0], np.cumsum(np.concatenate([z, z]))]); ns = np.arange(1, NMAX + 1)
    Q = np.array([np.percentile(cs[np.arange(N) + n] - cs[np.arange(N)], [2.5, 25, 50, 75, 97.5]) for n in ns]); band[tr] = (ns, Q)
    lo, hi = Q[d.n_genes_twas.values - 1, 0], Q[d.n_genes_twas.values - 1, 4]
    d = d.assign(trait=tr, outside_95=(d.z < lo) | (d.z > hi)); pts.append(d[["trait", "CNV_ID", "cls", "n_genes_twas", "z", "outside_95"]])
    for grp in ["mirror", "asymmetric"]:
        s = d[d.cls == grp]; r, p = stats.spearmanr(s.n_genes_twas, s.z)
        rows.append(dict(trait=tr, cls=grp, n_cnv=len(s), median_genes=s.n_genes_twas.median(), spearman_rho_genes_vs_z=r, p_asymptotic=p,
                         n_outside_null_95=int(s.outside_95.sum()), n_below_null_95=int((s.z < lo[d.cls == grp]).sum())))
P = pd.concat(pts); P.to_csv(HERE + "/genecount_vs_summed_z_points.csv", index=False)
R = pd.DataFrame(rows); R.to_csv(HERE + "/genecount_vs_summed_z_stats.csv", index=False); print(R.round(3).to_string(index=False))

INK = '#1a1a1a'; GREY = '#8a8a8a'; LIGHT = '#dcdcdc'; RED = '#b03a2e'; BLUE = '#1f5c8b'
plt.rcParams.update({"font.family": "DejaVu Sans", "axes.linewidth": .6, "xtick.major.width": .6, "ytick.major.width": .6})
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.3), gridspec_kw=dict(wspace=.22)); fig.subplots_adjust(left=.08, right=.98, top=.86, bottom=.2)
for ax, (tr, d) in zip(axs, D.items()):
    ns, Q = band[tr]
    ax.fill_between(ns, Q[:, 0], Q[:, 4], color=LIGHT, alpha=.6, lw=0, zorder=1); ax.fill_between(ns, Q[:, 1], Q[:, 3], color=LIGHT, lw=0, zorder=1)
    ax.plot(ns, Q[:, 2], color=GREY, lw=.6, zorder=2); ax.axhline(0, color=GREY, lw=.4, ls=":", zorder=1)
    p = P[P.trait == tr]
    for grp, c, mk in [("mirror", BLUE, "o"), ("asymmetric", RED, "s")]:
        s = p[p.cls == grp]; ax.scatter(s.n_genes_twas, s.z, s=22, color=c, marker=mk, edgecolor="white", lw=.5, zorder=4, label=grp)
        st = R[(R.trait == tr) & (R.cls == grp)].iloc[0]
        ax.text(.02, .97 - .07 * (grp == "asymmetric"), f"{grp}: n = {st.n_cnv}, ρ = {st.spearman_rho_genes_vs_z:+.2f} (P = {st.p_asymptotic:.2g})",
                transform=ax.transAxes, size=6, color=c, va="top")
    for _, r in p[(p.z.abs() >= 4) | (p.n_genes_twas >= 15)].iterrows():
        off = {("BMI", "2q11.2"): (-5, 5), ("BMI", "10q11.21_11.23"): (-5, -5), ("Height", "22q11.21_A-D"): (-4, 7), ("BMI", "22q11.21_A-D"): (-4, -7)}.get((tr, r.CNV_ID), (4, 0))
        ax.annotate(r.CNV_ID.replace("_", " "), (r.n_genes_twas, r.z), xytext=off, textcoords="offset points", size=5, color=INK, va="center",
                    ha="right" if off[0] < 0 else "left")
    ax.set_xlim(0, NMAX); ax.set_title(tr, size=8, weight="bold"); ax.set_xlabel("TWAS-modelled genes in CNV", size=7)
    ax.tick_params(labelsize=6, length=2); ax.spines[["top", "right"]].set_visible(False)
axs[0].set_ylabel("summed TWAS z", size=7)
fig.suptitle("Gene count vs summed TWAS z per CNV", size=8.5, weight="bold", y=.97)
fig.text(.08, .01, "Grey: null for n consecutive TWAS genes anywhere in the genome (light 95%, dark 50%, line median). Class = height class for both traits "
         "(no CNV is asymmetric for BMI).\nLabelled: |summed z| ≥ 4 or ≥ 15 genes; no CNV lies outside the 95% band. ρ: Spearman, gene count vs summed z within class (asymptotic P, descriptive).",
         size=5.3, color=GREY, va="bottom")
fig.savefig(HERE + "/fig_genecount_vs_summed_z.png", dpi=300); fig.savefig(HERE + "/fig_genecount_vs_summed_z.pdf")
print(P[P.outside_95].to_string(index=False))
