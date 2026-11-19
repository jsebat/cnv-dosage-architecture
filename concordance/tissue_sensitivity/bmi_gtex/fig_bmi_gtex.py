# Figures for bmi_gtex_analysis.py (MS 2026-11-03). Usage: python fig_bmi_gtex.py [gtex|gtex_pec]
# Tissues with < MIN_GENES CNV genes with a model are left out of the sensitivity figure and of the Bonferroni count
# (they still contribute to the all-tissue mean).
#  fig_bmi_gtex_scatter: BMI dosage effect (DUP - DEL)/2 vs summed GTEx-only BMI TWAS z per CNV (panel A style; all CNVs mirror
#    for BMI; Theil-Sen line through all CNVs).
#  fig_bmi_gtex_tissue_sensitivity: Spearman rho with z from one GTEx tissue only (A) and from all GTEx tissues but one (B);
#    * one-sided circular P < 0.05, ** P < 0.05 / n tissues with models in these CNVs; (n) = CNV genes with an expression model.
import os, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import sys
MODE = sys.argv[1] if len(sys.argv) > 1 else "gtex"; MIN_GENES = 5
SRC = {"gtex": "GTEx", "gtex_pec": "GTEx + PsychENCODE"}[MODE]
HERE = os.path.dirname(os.path.abspath(__file__)); pre = f"{HERE}/bmi_{MODE}_"; fpre = f"{HERE}/fig_bmi_{MODE}_"
P = pd.read_csv(pre + "points.csv"); ST = pd.read_csv(pre + "stats.csv").set_index("group")
R = pd.read_csv(pre + "tissue_rho.csv"); GC = pd.read_csv(pre + "gene_counts.csv", index_col=0)
NMISS = 22 - len(P)
INK = "#1f1f1e"; MUTED = "#6b6a66"; GRID = "#e4e3df"; ALL = "#52514e"; BLUE = "#1f5c8b"; LIGHT = "#dcdcdc"
fp = lambda p: f"{p:.2g}"
plt.rcParams.update({"font.family": "Arial", "axes.edgecolor": MUTED, "axes.linewidth": 0.6, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.labelcolor": INK, "text.color": INK})

# ---------------- scatter
fig, ax = plt.subplots(figsize=(3.6, 3.0)); fig.subplots_adjust(left=.19, right=.96, top=.84, bottom=.2)
ax.scatter(P.summed_bmi_twas_z, P.dosage_effect_bmi, s=18, color=BLUE, edgecolor="white", lw=.4, zorder=4)
xr = np.array([P.summed_bmi_twas_z.min(), P.summed_bmi_twas_z.max()]); a = ST.loc["all"]
ax.plot(xr, a.theilsen_intercept + a.theilsen_slope * xr, color=INK, lw=.9, zorder=3)
ax.axhline(0, color=LIGHT, lw=.6, zorder=1); ax.axvline(0, color=LIGHT, lw=.6, zorder=1)
LAB_GTEX = {"16p11.2_BP1-BP3": ("16p11.2 BP1\u2013BP3", 4, 3, "left"), "16p11.2_BP4_BP5": ("16p11.2 BP4\u2013BP5", 4, 3, "left"),
       "22q11.21_A-D": ("22q11.2 A\u2013D", -5, 3, "right"), "7q11.23_WBS": ("Williams", -5, -3, "right"), "3q29": ("3q29", -5, 0, "right"),
       "17q12": ("17q12", 5, 0, "left"), "2q11.2": ("2q11.2", -5, 0, "right"), "1q21.1": ("1q21.1", 0, 6, "center"),
       "10q11.21_11.23": ("10q11.21\u201311.23", 0, -9, "center")}
LAB_PEC = {"16p11.2_BP1-BP3": ("16p11.2 BP1\u2013BP3", -6, 0, "right"), "16p11.2_BP4_BP5": ("16p11.2 BP4\u2013BP5", 6, 0, "left"),
           "3q29": ("3q29", -6, 0, "right"), "17q12": ("17q12", 6, 0, "left"), "22q11.21_A-D": ("22q11.2 A\u2013D", 4, 7, "left"),
           "2q11.2": ("2q11.2", -6, 0, "right"), "7q11.23_WBS": ("Williams", 6, 0, "left"), "10q11.21_11.23": ("10q11.21\u201311.23", 0, -8, "center"),
           "1q21.1": ("1q21.1", 6, 0, "left")}
lab = LAB_PEC if MODE == "gtex_pec" else LAB_GTEX
for cid, (t, dx, dy, ha) in lab.items():
    r = P[P.CNV_ID == cid]
    if len(r): r = r.iloc[0]; ax.annotate(t, (r.summed_bmi_twas_z, r.dosage_effect_bmi), xytext=(dx, dy), textcoords="offset points", size=5.2, color=BLUE, ha=ha, va="center")
w = ST.loc["all without Williams"]
txt = f"All CNVs (n = {int(a.n_cnv)}): ρ = {a.spearman_rho:.2f}, P = {fp(a.p_one_sided)}\nwithout Williams: ρ = {w.spearman_rho:.2f}, P = {fp(w.p_one_sided)}"
if MODE == "gtex_pec":
    rr = R.set_index(["analysis", "tissue"])
    txt += f"\nPsychENCODE only: ρ = {rr.loc[('only', 'PsychENCODE'), 'rho_all']:.2f}, P = {fp(rr.loc[('only', 'PsychENCODE'), 'p1_all'])}" \
           f"\nGTEx only: ρ = {rr.loc[('drop', 'PsychENCODE'), 'rho_all']:.2f}, P = {fp(rr.loc[('drop', 'PsychENCODE'), 'p1_all'])}"
ax.text(.02, .98, txt, transform=ax.transAxes, size=5.4, va="top", color=INK, linespacing=1.35)
ax.set_ylim(P.dosage_effect_bmi.min() - 0.08, P.dosage_effect_bmi.max() + (0.38 if MODE == "gtex_pec" else 0.18))
ax.set_xlabel(f"summed BMI TWAS z of genes in the CNV\n({SRC}; + = higher expression, higher BMI)", size=6.4)
ax.set_ylabel("CNV dosage effect on BMI (SD)\n(DUP − DEL) / 2", size=6.4); ax.tick_params(labelsize=5.6, length=2, pad=1.5)
ax.spines[["top", "right"]].set_visible(False)
ax.set_title(f"BMI: {SRC}-predicted expression of CNV genes\nvs CNV dosage effect", size=7, weight="bold")
fig.text(.01, .01, "P: one-sided circular permutation (10,000)." + (f" {NMISS} of 22 CNVs have no gene with a {SRC} model and are not shown." if NMISS else "")
         + (" PsychENCODE z from Tan et al. 2024 BMI GWAS; GTEx z from TWAS30." if MODE == "gtex_pec" else ""), size=4.4, color=MUTED)
fig.savefig(fpre + "scatter.png", dpi=400); fig.savefig(fpre + "scatter.pdf"); plt.close(fig)

# ---------------- tissue sensitivity
COL = {"brain": "#2a78d6", "blood": "#eb6834", "adipose": "#1baf7a", "other_GTEx": "#eda100"}
GLAB = {"brain": "Brain", "blood": "Blood", "adipose": "Adipose", "other_GTEx": "Other GTEx"}
GROUP = lambda s: "brain" if ("Brain" in s or s == "PsychENCODE" or s == "GTEX_Pituitary") else "blood" if s == "GTEX_Whole_Blood" else "adipose" if "Adipose" in s else "other_GTEx"
nice = lambda t: "PsychENCODE cortex" if t == "PsychENCODE" else t.replace("GTEX_", "").replace("_basal_ganglia", "").replace("_", " ")
base = R[R.analysis == "all"].iloc[0]
D = R[R.analysis == "only"].set_index("tissue")[["rho_all", "p1_all"]].rename(columns={"rho_all": "single_rho", "p1_all": "single_p"})
D = D.join(R[R.analysis == "drop"].set_index("tissue")[["rho_all", "p1_all"]].rename(columns={"rho_all": "drop_rho", "p1_all": "drop_p"}))
D = D.join(GC.drop(index="ALL").astype(int)); D["group"] = [GROUP(t) for t in D.index]
N_EXCL = int((D.n_genes_single < MIN_GENES).sum()); D = D[D.n_genes_single >= MIN_GENES]
gorder = [g for g in ["blood", "brain", "adipose", "other_GTEx"] if g in set(D.group)]; D["g"] = D.group.map(gorder.index)
D = D.sort_values(["g", "n_genes_single"], ascending=[True, False]); N_ALL = int(GC.loc["ALL", "n_genes_single"])
NTEST = int(D.single_rho.notna().sum()); BONF = 0.05 / NTEST
star = lambda p: "" if np.isnan(p) else "**" if p < BONF else "*" if p < 0.05 else ""
xs, x, prev = [], 0.0, None
for g in D.group:
    if prev is not None and g != prev: x += 0.8
    xs.append(x); x += 1; prev = g
xs = np.array(xs); x_all = x + 0.8; cols = [COL[g] for g in D.group]
plt.rcParams.update({"font.size": 8})
fig, axs = plt.subplots(2, 1, figsize=(max(7.0, 0.3 * (len(D) + 1) + 1.6), 6.6), sharex=True, gridspec_kw=dict(hspace=0.22))
def bars(ax, v, p, n, extra=None):
    ax.bar(xs, np.nan_to_num(v), width=0.78, color=cols, linewidth=0, zorder=2,
           hatch=["////" if t == "PsychENCODE" else "" for t in D.index], edgecolor="white")
    for xi, vi, pi in zip(xs, v, p):
        if np.isnan(vi): ax.text(xi, 0.02, "n/a", ha="center", va="bottom", fontsize=6, color=MUTED, rotation=90); continue
        if star(pi): ax.text(xi, vi + (0.01 if vi >= 0 else -0.04), star(pi), ha="center", va="bottom" if vi >= 0 else "top", fontsize=9, color=INK)
    xx, nn = list(xs), list(n)
    if extra is not None:
        ax.bar([x_all], [extra[0]], width=0.78, color=ALL, linewidth=0, zorder=2)
        ax.text(x_all, extra[0] + 0.01, star(extra[1]), ha="center", va="bottom", fontsize=9, color=INK); xx.append(x_all); nn.append(N_ALL)
    for xi, ni in zip(xx, nn): ax.text(xi, -0.97, f"({ni})", ha="center", va="bottom", fontsize=6, color=MUTED, rotation=90)
    ax.axhline(0, color=MUTED, lw=0.6, zorder=3); ax.grid(axis="y", color=GRID, lw=0.5, zorder=0)
    for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
    ax.set_ylim(-1.0, 0.95); ax.set_ylabel("Spearman ρ"); ax.axhline(base.rho_all, color=ALL, ls="--", lw=0.8, zorder=1)
a = axs[0]; bars(a, D.single_rho.values, D.single_p.values, D.n_genes_single)
a.text(x_all, base.rho_all + 0.03, f"all\nρ = {base.rho_all:.2f}", ha="center", va="bottom", fontsize=7, color=ALL)
a.set_title(f"A   BMI: z-scores from a single tissue only ({SRC})", loc="left", fontsize=10, fontweight="bold", pad=18)
a = axs[1]; bars(a, D.drop_rho.values, D.drop_p.values, D.n_genes_all_but, (base.rho_all, base.p1_all))
a.set_title(f"B   BMI: z-scores averaged over all tissues except that one ({SRC})", loc="left", fontsize=10, fontweight="bold", pad=6)
a.set_xticks(list(xs) + [x_all]); a.set_xticklabels([nice(t) for t in D.index] + ["All tissues"], rotation=60, ha="right", rotation_mode="anchor", fontsize=7)
for lab_ in a.get_xticklabels():
    if lab_.get_text() == "All tissues": lab_.set_fontweight("bold"); lab_.set_color(INK)
a.set_xlim(-0.8, x_all + 0.8)
for g in gorder:
    gx = xs[D.group.values == g]
    axs[0].annotate("", xy=(gx.min() - 0.4, 0.93), xytext=(gx.max() + 0.4, 0.93), arrowprops=dict(arrowstyle="-", color=COL[g], lw=2))
    axs[0].text(gx.mean(), 0.95, GLAB[g], ha="center", va="bottom", fontsize=7.5, color=INK)
import textwrap
fig.text(0.01, 0.003, "\n".join(textwrap.fill(l, 175) for l in (f"BMI; non-overlapping CNVs with ≥ 1 modelled gene (n = {int(ST.loc['all'].n_cnv)}; all mirror for BMI). ρ: Spearman, summed BMI TWAS z ({SRC}; GTEx from TWAS30" + (", PsychENCODE (hatched) from Tan et al. 2024 GWAS" if MODE == "gtex_pec" else "") + ") vs dosage effect (DUP − DEL)/2. "
         f"P: one-sided circular permutation (10,000).\n* P < 0.05; ** P < 0.05/{NTEST} = {BONF:.4f} (Bonferroni over the {NTEST} tissues shown). {N_EXCL} tissues with < {MIN_GENES} CNV genes with a model not shown (still in the all-tissue mean).\n"
         f"(n): CNV genes with an expression model in that tissue (A) or in any tissue but that one (B); {N_ALL} in any tissue.").split("\n")),
         fontsize=6.5, color=MUTED, va="bottom")
fig.subplots_adjust(left=0.08, right=0.99, top=0.91, bottom=0.33)
fig.savefig(fpre + "tissue_sensitivity.png", dpi=300); fig.savefig(fpre + "tissue_sensitivity.pdf")
D.drop(columns="g").to_csv(pre + "tissue_sensitivity_table.csv"); print("done")
