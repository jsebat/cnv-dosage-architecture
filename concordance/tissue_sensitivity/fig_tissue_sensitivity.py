# Tissue sensitivity of Fig. panel A (height; primary non-overlapping CNV set, all CNVs without Williams, n = 22).
# Top: z from a single TWAS30 panel only. Bottom: z averaged over all panels except that one; rightmost bar = all panels (as in
# the figure). Bars: Spearman rho; * one-sided circular-permutation P < 0.05, ** P < 0.05 / n tissues with models (Bonferroni).
# Inputs: tissue_drivers_rho.csv, tissue_drivers_cov_decomp.csv (from tissue_drivers.py). Run: python fig_tissue_sensitivity.py
import os, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
HERE = os.path.dirname(os.path.abspath(__file__))
R = pd.read_csv(f"{HERE}/tissue_drivers_rho.csv"); D = pd.read_csv(f"{HERE}/tissue_drivers_cov_decomp.csv", index_col=0)
COL = {"brain": "#2a78d6", "blood": "#eb6834", "adipose": "#1baf7a", "other_GTEx": "#eda100"}
GLAB = {"brain": "Brain", "blood": "Blood", "adipose": "Adipose", "other_GTEx": "Other GTEx"}
ALL = "#52514e"; INK = "#1f1f1e"; MUTED = "#6b6a66"; GRID = "#e4e3df"
def nice(t):
    s = {"YFS": "YFS (blood)", "NTR": "NTR (blood)", "CMC": "CMC (DLPFC)", "MET": "METSIM (adipose)"}.get(t)
    return s or t.replace("GTEX_", "").replace("_basal_ganglia", "").replace("_", " ")
get = lambda a: R.set_index("analysis").loc[a]
D["single_rho"] = [get(f"only {t}").rho_all_noW for t in D.index]; D["single_p"] = [get(f"only {t}").p_all_noW for t in D.index]
D["drop_rho"] = [get(f"drop {t}").rho_all_noW for t in D.index]; D["drop_p"] = [get(f"drop {t}").p_all_noW for t in D.index]
D.loc[D.single_rho.isna(), "single_p"] = np.nan   # no gene models in any panel-A CNV: undefined
GC = pd.read_csv(f"{HERE}/tissue_gene_counts.csv", index_col=0)   # from tissue_gene_counts.py
D["n_genes_single"] = GC.n_genes_single.reindex(D.index).astype(int); D["n_genes_all_but"] = GC.n_genes_all_but.reindex(D.index).astype(int); N_ALL = int(GC.loc["ALL", "n_genes_single"])
gorder = ["blood", "brain", "adipose", "other_GTEx"]
D["g"] = D.group.map(gorder.index); D = D.sort_values(["g", "n_gene_models_in_cnvs"], ascending=[True, False])
base = get("baseline (all panels)")
# x positions with a gap between tissue groups and before the all-tissues bar
xs, x, prev = [], 0.0, None
for g in D.group:
    if prev is not None and g != prev: x += 0.8
    xs.append(x); x += 1; prev = g
xs = np.array(xs); x_all = x + 0.8
labels = [nice(t) for t in D.index]
cols = [COL[g] for g in D.group]
plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.edgecolor": MUTED, "axes.linewidth": 0.6, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.labelcolor": INK, "text.color": INK})
NTEST = int(D.single_rho.notna().sum()); BONF = 0.05 / NTEST   # tissues with gene models in these CNVs
star = lambda p: "" if np.isnan(p) else "**" if p < BONF else "*" if p < 0.05 else ""
fig, axs = plt.subplots(2, 1, figsize=(11, 6.4), sharex=True, gridspec_kw=dict(hspace=0.22))
def bars(ax, v, p, extra=None):
    ax.bar(xs, np.nan_to_num(v), width=0.78, color=cols, linewidth=0, zorder=2)
    for xi, vi, pi in zip(xs, v, p):
        if np.isnan(vi): ax.text(xi, 0.02, "n/a", ha="center", va="bottom", fontsize=6, color=MUTED, rotation=90); continue
        if star(pi): ax.text(xi, vi + (0.01 if vi >= 0 else -0.04), star(pi), ha="center", va="bottom" if vi >= 0 else "top", fontsize=9, color=INK)
    if extra is not None:
        ax.bar([x_all], [extra[0]], width=0.78, color=ALL, linewidth=0, zorder=2)
        ax.text(x_all, extra[0] + 0.01, star(extra[1]), ha="center", va="bottom", fontsize=9, color=INK)
    ax.axhline(0, color=MUTED, lw=0.6, zorder=3); ax.grid(axis="y", color=GRID, lw=0.5, zorder=0)
    for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
    ax.set_ylim(-1.0, 0.95); ax.set_ylabel("Spearman ρ")
    ax.axhline(base.rho_all_noW, color=ALL, ls="--", lw=0.8, zorder=1)
def counts(ax, n, extra=None):   # (number of CNV genes with an expression model) along the bottom of the panel
    for xi, ni in zip(list(xs) + ([x_all] if extra is not None else []), list(n) + ([extra] if extra is not None else [])):
        ax.text(xi, -0.97, f"({ni})", ha="center", va="bottom", fontsize=6, color=MUTED, rotation=90)
a = axs[0]; bars(a, D.single_rho.values, D.single_p.values); counts(a, D.n_genes_single)
a.text(xs[-1] + 0.6, base.rho_all_noW + 0.03, f"all tissues ρ = {base.rho_all_noW:.2f}", ha="right", va="bottom", fontsize=7, color=ALL)
a.set_title("A   z-scores from a single tissue only", loc="left", fontsize=10, fontweight="bold", pad=18)
a = axs[1]; bars(a, D.drop_rho.values, D.drop_p.values, (base.rho_all_noW, base.p_all_noW)); counts(a, D.n_genes_all_but, N_ALL)
a.set_title("B   z-scores averaged over all tissues except that one", loc="left", fontsize=10, fontweight="bold", pad=6)
a.set_xticks(list(xs) + [x_all]); a.set_xticklabels(labels + ["All tissues"], rotation=60, ha="right", rotation_mode="anchor", fontsize=7)
for lab in a.get_xticklabels():
    if lab.get_text() == "All tissues": lab.set_fontweight("bold"); lab.set_color(INK)
a.set_xlim(-0.8, x_all + 0.8)
# group brackets above top panel
for g in gorder:
    gx = xs[D.group.values == g]; axs[0].annotate("", xy=(gx.min() - 0.4, 0.93), xytext=(gx.max() + 0.4, 0.93), xycoords="data",
                                                  arrowprops=dict(arrowstyle="-", color=COL[g], lw=2))
    axs[0].text(gx.mean(), 0.95, GLAB[g], ha="center", va="bottom", fontsize=7.5, color=INK)
fig.legend(handles=[Patch(color=COL[g], label=GLAB[g]) for g in gorder] + [Patch(color=ALL, label="All tissues (as in figure)")],
           loc="upper right", ncol=5, frameon=False, fontsize=7.5, bbox_to_anchor=(0.99, 0.995))
fig.text(0.01, 0.003, "Height; non-overlapping CNVs excluding Williams (n = 22), as in panel A. ρ: Spearman, summed TWAS z (TWAS30) vs dosage effect (DUP − DEL)/2. "
         f"P: one-sided circular permutation (10,000).\n* P < 0.05; ** P < 0.05/{NTEST} = {BONF:.4f} (Bonferroni over the {NTEST} tissues with gene models in these CNVs). (n): CNV genes with an expression model in that tissue (A) or in any tissue but that one (B); {N_ALL} in any tissue. n/a: tissue has no gene models in these CNVs.",
         fontsize=6.5, color=MUTED, va="bottom")
fig.subplots_adjust(left=0.06, right=0.99, top=0.94, bottom=0.31)
fig.savefig(f"{HERE}/fig_tissue_sensitivity.png", dpi=300); fig.savefig(f"{HERE}/fig_tissue_sensitivity.pdf")
D[["group", "n_genes_single", "n_genes_all_but", "single_rho", "single_p", "drop_rho", "drop_p", "share_cov"]].to_csv(f"{HERE}/tissue_sensitivity_table.csv")
print("done")
