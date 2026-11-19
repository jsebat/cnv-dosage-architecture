# TWAS z (TWAS30, mean over tissues) for 16p11.2 BP4-BP5, height only (MS 2026-11-08).
# Shows only genes with a TWAS model, over the same contiguous span used for the TWAS dosage panel (twas_dosage_panel_stats.py):
# every TWAS-modelled gene on chr16 from QPRT to CORO1A, ordered by TxS. This includes RP11-455F5.3 (lncRNA, not in the CNV gene list),
# so the total matches summed_twas_z in twas_dosage_panel.csv (+2.45). lncRNAs are kept.
# Count of CNV-list genes with no TWAS model (cnv_genes_annotated.csv, 46 genes) is printed on the plot, split into lncRNAs
# (AC... clone IDs, *-AS1) and other genes.
# Outputs: fig_16p11.2_BP4_BP5_twas_z.{png,pdf}, 16p11.2_BP4_BP5_twas_z.csv
import os, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); CONC = os.path.dirname(HERE); ROOT = os.path.dirname(os.path.dirname(CONC))
cnv = pd.read_csv(ROOT + "/pathway_analysis/Final/cnv_genes_annotated.csv"); cnv = cnv[cnv.CNV_ID == "16p11.2_BP4_BP5"].Gene
tw = pd.read_csv(f"{CONC}/twas30/height.txt", sep=r"\s+")
tw = tw[tw.CHR.astype(str).str.replace("chr", "") == "16"]
g = tw.groupby("GENE").agg(TxS=("TxS", "min"), z_height=("TWAS.Z", "mean"), n_tissues=("TWAS.Z", "size")).reset_index().rename(columns={"GENE": "Gene"})
listed = g[g.Gene.isin(cnv)]
g = g[g.TxS.between(listed.TxS.min(), listed.TxS.max())].sort_values("TxS").reset_index(drop=True)
g["in_cnv_gene_list"] = g.Gene.isin(cnv)
g.to_csv(HERE + "/16p11.2_BP4_BP5_twas_z.csv", index=False)
none = cnv[~cnv.isin(g.Gene)]; is_lnc = none.str.match(r"^AC\d+\.\d+$") | none.str.endswith("-AS1")
n_lnc, n_oth = int(is_lnc.sum()), int((~is_lnc).sum())

INK = '#1a1a1a'; LIGHT = '#dcdcdc'; UP = '#b03a2e'; DN = '#1f5c8b'
plt.rcParams.update({"font.family": "DejaVu Sans", "axes.linewidth": .6, "xtick.major.width": .6, "ytick.major.width": .6})
fig, ax = plt.subplots(figsize=(4, 3.2)); fig.subplots_adjust(left=.15, right=.98, top=.87, bottom=.2)
x = np.arange(len(g)); z = g.z_height.values; tot, mean = z.sum(), z.mean()
ylim = max(3, np.abs(z).max() * 1.3)
ax.bar(x, z, width=.7, color=np.where(z > 0, UP, DN), zorder=3)
ax.axhline(0, color=INK, lw=.5, zorder=4)
for v in (1.96, -1.96): ax.axhline(v, color=LIGHT, lw=.5, zorder=1)
ax.text(.02, .97, f"Height: total z = {tot:+.2f} ({len(g)} genes)", transform=ax.transAxes, size=7.5, weight="bold", va="top", color=INK)
ax.text(.02, .03, f"No TWAS z: {n_oth} genes, {n_lnc} lncRNAs", transform=ax.transAxes, size=6.5, va="bottom", color=INK)
ax.set_ylim(-ylim, ylim); ax.set_xlim(-.6, len(g) - .4); ax.set_ylabel("TWAS z (mean over tissues)", size=7)
ax.tick_params(labelsize=6, length=2); ax.spines[["top", "right"]].set_visible(False)
ax.set_xticks(x); ax.set_xticklabels(g.Gene, rotation=45, ha="right", rotation_mode="anchor", size=6)
fig.suptitle("16p11.2 BP4–BP5: TWAS z (TWAS30)", size=9, weight="bold", y=.97)
print(g.to_string()); print("total", round(tot, 3), "| no TWAS:", n_oth, "genes", n_lnc, "lncRNAs")
fig.savefig(HERE + "/fig_16p11.2_BP4_BP5_twas_z.png", dpi=300); fig.savefig(HERE + "/fig_16p11.2_BP4_BP5_twas_z.pdf")
