# Main-text ConGA figure of the CNV GWAS paper (UK Biobank): signal groups by LC/HC x DEL/DUP/dose-aligned (a), HC dose-aligned z
# with its range over the 14 locus-removal runs (b), and the top HC dose-aligned cell types, tissues and brain regions (c).
# Draws from finished aggregate tables only (conga_contrasts.py stage B output); no new tests.
# null model: genome shift of the gene annotation (2,000 shared shifts), matched in every locus-removal run (from the input run)
"""python fig_conga_main.py --run <conga_fig16 dir> --groups <signal groups dir> --out <png>
Marks: ** locus-robust (FDR < 5% and holds in all 14 removal runs), * FDR < 5% only. Signal groups are named after pathways only."""
import argparse
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from conga_contrasts import label as ctx_label   # same context names as Supplementary Fig. 1

TN = {"height": "Height", "bmi": "BMI"}
TC = {"height": "#D55E00", "bmi": "#0072B2"}
SG = {"height": {1: "Ossification; osteoblast differentiation", 2: "Ser/Thr kinase receptor (TGF-β family) signalling",
                 3: "Nuclear receptor transcription; hormone signalling", 4: "Viral carcinogenesis; cellular senescence",
                 5: "Respiratory and endocrine system development", 6: "Regulation of translation",
                 7: "Response to misfolded protein; stem cell differentiation", 8: "Cortisol synthesis and secretion",
                 9: "Chromatin organisation; cellular senescence", 10: "Muscle tissue and organ development",
                 11: "Negative regulation of cell development and growth"},
      "bmi": {1: "Second-messenger signalling", 2: "MicroRNAs in cancer; cytomegalovirus infection", 3: "β-catenin degradation; S phase",
              4: "Axon development; synapse assembly", 5: "Transcription factor binding", 6: "Glutamatergic synapse; synaptic membrane",
              7: "Synapse organisation", 8: "Phosphoprotein phosphatase activity", 9: "PTEN regulation; oocyte meiosis",
              10: "Neuron differentiation", 11: "Insulin secretion; dopaminergic synapse", 12: "Protein heterodimerisation; ER stress",
              13: "Oncogenic MAPK signalling", 14: "Learning; cognition", 15: "TP53 metabolic genes; longevity"}}
COLS = [("LC", "DEL"), ("LC", "DUP"), ("LC", "aligned"), ("HC", "DEL"), ("HC", "DUP"), ("HC", "aligned")]


def mark(r):
    if r.FDR_BH < 0.05:
        return "**" if bool(r.locus_robust) else "*"
    return ""


def panel_a(ax, T, t, lim):
    k = sorted(SG[t]); units = [f"SG:{t}:{i}" for i in k]
    Z = np.full((len(k), len(COLS)), np.nan); M = np.full(Z.shape, "", object)
    idx = T[T.trait == t].set_index(["unit", "constraint", "stat"])
    for i, u in enumerate(units):
        for j, (c, s) in enumerate(COLS):
            r = idx.loc[(u, c, s)]; Z[i, j] = r.z_emp; M[i, j] = mark(r)
    im = ax.imshow(Z, cmap="RdBu_r", norm=TwoSlopeNorm(0, -lim, lim), aspect="auto")
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            if M[i, j]:
                ax.text(j, i + 0.08, M[i, j], ha="center", va="center", fontsize=8, fontweight="bold")
    ax.set_yticks(range(len(k))); ax.set_yticklabels([f"SG{i}  {SG[t][i]}" for i in k], fontsize=7)
    SL = {"DEL": "DEL", "DUP": "DUP", "aligned": "dose-\naligned"}
    ax.set_xticks(range(len(COLS))); ax.set_xticklabels([c + "\n" + SL[s] for c, s in COLS], fontsize=7)
    ax.axvline(2.5, color="black", lw=1.2); ax.tick_params(length=0)
    ax.set_title(f"{TN[t]}: CNVs in {TN[t].lower() if t == 'height' else TN[t]} signal-group genes", fontsize=9, loc="left")
    return im


def panel_b(ax, T, H, t):
    k = sorted(SG[t]); y = np.arange(len(k))
    idx = T[(T.trait == t) & (T.constraint == "HC") & (T.stat == "aligned")].set_index("unit")
    h = H[(H.trait == t) & (H.stat == "aligned")]
    for i, kk in enumerate(k):
        u = f"SG:{t}:{kk}"; r = idx.loc[u]; sig = r.FDR_BH < 0.05
        zr = h[h.set_id == f"{u}|HC"].z_emp_run
        col = TC[t] if sig else "#9a9a9a"
        if len(zr):
            ax.plot([zr.min(), zr.max()], [i, i], color=col, lw=2.2, alpha=0.45, solid_capstyle="round")
        ax.scatter([r.z_emp], [i], s=36, color=col if bool(r.locus_robust) else "white", edgecolors=col, lw=1.3, zorder=3)
        if not len(zr):
            ax.text(r.z_emp + (0.25 if r.z_emp >= 0 else -0.25), i, "no removal runs", fontsize=5.5, va="center",
                    ha="left" if r.z_emp >= 0 else "right", color="#777777")
    ax.axvline(0, color="black", lw=0.6)
    ax.set_yticks(y); ax.set_yticklabels([f"SG{i}" for i in k], fontsize=7); ax.invert_yaxis()
    ax.set_xlabel("HC dose-aligned z (bar: range over 14 locus-removal runs)", fontsize=7); ax.tick_params(labelsize=7)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title(f"{TN[t]}: HC genes, deletions and duplications combined", fontsize=9, loc="left")


def panel_c(ax, T, names, t, n_per=5):
    x = T[(T.trait == t) & (T.category != "signal group")]
    hc = x[(x.constraint == "HC") & (x.stat == "aligned")]
    top = pd.concat([hc[hc.category == c].sort_values("P_emp").head(n_per) for c in ["cell type", "tissue", "brain region"]])
    lc = x[(x.constraint == "LC") & (x.stat == "aligned")].set_index("unit")
    y = np.arange(len(top))
    ax.barh(y, top.z_emp, color=[TC[t] if f < 0.05 else "#bbbbbb" for f in top.FDR_BH], height=0.7, label="HC dose-aligned")
    ax.scatter(lc.loc[top.unit].z_emp, y, marker="o", s=16, facecolors="white", edgecolors="black", lw=0.7, zorder=3,
               label="LC dose-aligned")
    for yy, r in zip(y, top.itertuples()):
        m = mark(r)
        if m:
            ax.text(r.z_emp + (0.15 if r.z_emp > 0 else -0.15), yy, m, va="center", ha="left" if r.z_emp > 0 else "right", fontsize=7)
    for b in (n_per - 0.5, 2 * n_per - 0.5):
        ax.axhline(b, color="black", lw=0.5, ls=":")
    ax.set_yticks(y); ax.set_yticklabels([ctx_label(u, names) for u in top.unit], fontsize=6.3); ax.invert_yaxis()
    ax.axvline(0, color="black", lw=0.6); ax.set_xlabel("dose-aligned z vs genome-shift null", fontsize=7, labelpad=1); ax.tick_params(labelsize=7)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(fontsize=6.5, frameon=False, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.09))
    ax.set_title(f"{TN[t]}: top cell types, tissues and brain regions (HC dose-aligned)", fontsize=9, loc="left")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True); ap.add_argument("--groups", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    T = pd.read_csv(f"{a.run}/conga_contrasts_full.tsv.gz", sep="\t"); T = T[T.constraint.isin(["LC", "HC"])]
    T["locus_robust"] = T.locus_robust.astype(str).str.lower().eq("true")
    H = pd.read_csv(f"{a.run}/conga_contrasts_removal_runs.tsv.gz", sep="\t")
    names = {}
    for t in ["height", "bmi"]:
        g = pd.read_csv(f"{a.groups}/{t}.groups.tsv", sep="\t"); names.update(dict(zip(g.set_group, g.name)))
    fig = plt.figure(figsize=(14, 16))
    gs = fig.add_gridspec(3, 2, height_ratios=[15, 15, 15], hspace=0.45, wspace=0.95, left=0.22, right=0.97, top=0.94, bottom=0.065)
    lim = 6.2
    for j, t in enumerate(["height", "bmi"]):
        im = panel_a(fig.add_subplot(gs[0, j]), T, t, lim)
        panel_b(fig.add_subplot(gs[1, j]), T, H, t)
        panel_c(fig.add_subplot(gs[2, j]), T, names, t)
    cax = fig.add_axes([0.975, 0.70, 0.008, 0.2]); cb = fig.colorbar(im, cax=cax); cb.ax.tick_params(labelsize=6)
    cb.set_label("z (red: trait higher)", fontsize=7)
    for y, l in [(0.955, "a"), (0.64, "b"), (0.325, "c")]:
        fig.text(0.01, y, l, fontsize=14, fontweight="bold")
    fig.suptitle("Person-level ConGA of CNV dosage in signal groups, UK Biobank (All of Us replication pending)", fontsize=11, x=0.5, y=0.985)
    fig.text(0.01, 0.005, "** locus-robust (FDR < 5% and holds in all 14 locus-removal runs); * FDR < 5% only. b: filled, locus-robust; "
             "open, not; grey, FDR ≥ 5%.", fontsize=7)
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True); fig.savefig(a.out, dpi=250); plt.close(fig)


if __name__ == "__main__":
    main()
