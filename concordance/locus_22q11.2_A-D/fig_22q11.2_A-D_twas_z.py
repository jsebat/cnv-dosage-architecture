# TWAS z (TWAS30, mean over tissues) for every gene of 22q11.21 A-D (cnv_genes_annotated.csv, CNV_ID 22q11.21_A-D), height only (MS 2026-11-08).
# Same layout as ../locus_16p11.2_BP4_BP5/fig_16p11.2_BP4_BP5_twas_z.py. Genes ordered by GRCh38 start (Ensembl REST overlap chr22:18.9-21.2Mb,
# saved as ensembl_GRCh38_chr22_18.9-21.2Mb.json); SEPT5 = SEPTIN5, USP41 = USP41P. SNORA15 (the Ensembl symbol now maps to chr7; no TWAS model)
# is put at the right end. Long non-coding RNAs are excluded: AC... clone IDs, LINC*, and THAP7-AS1 (THAP7-AS1 has a TWAS model, mean height z
# = -0.31 over 20 tissues, so dropping it changes the total by +0.31). TMEM191A (HGNC: pseudogene) is kept.
# Genes with no TWAS model: empty bar. Mean z over genes with a model; total = sum over listed genes.
# Outputs: fig_22q11.2_A-D_twas_z.{png,pdf}, 22q11.2_A-D_twas_z.csv
import os, json, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); CONC = os.path.dirname(HERE); ROOT = os.path.dirname(os.path.dirname(CONC))
E = pd.DataFrame(json.load(open(HERE + "/ensembl_GRCh38_chr22_18.9-21.2Mb.json"))).dropna(subset=["external_name"])
pos = E.groupby("external_name").start.min().to_dict(); pos["SEPT5"] = pos["SEPTIN5"]; pos["USP41"] = pos["USP41P"]
g = pd.read_csv(ROOT + "/pathway_analysis/Final/cnv_genes_annotated.csv"); g = g[g.CNV_ID == "22q11.21_A-D"][["Gene"]].copy()
lnc = g.Gene.str.match(r"^AC\d+\.\d+$") | g.Gene.str.startswith("LINC") | g.Gene.str.endswith("-AS1")
print("dropped lncRNAs:", " ".join(g.Gene[lnc])); g = g[~lnc]
g["pos"] = g.Gene.map(pos); g = pd.concat([g[g.pos.notna()].sort_values("pos"), g[g.pos.isna()].sort_values("Gene")]).reset_index(drop=True)
tw = pd.read_csv(f"{CONC}/twas30/height.txt", sep=r"\s+"); m = tw.groupby("GENE")["TWAS.Z"].agg(["mean", "size"])
g["z_height"] = g.Gene.map(m["mean"]); g["n_tissues_height"] = g.Gene.map(m["size"]).fillna(0).astype(int)
g.to_csv(HERE + "/22q11.2_A-D_twas_z.csv", index=False)
n_pos = int(g.pos.notna().sum())
INK = '#1a1a1a'; GREY = '#8a8a8a'; LIGHT = '#dcdcdc'; UP = '#b03a2e'; DN = '#1f5c8b'
plt.rcParams.update({"font.family": "DejaVu Sans", "axes.linewidth": .6, "xtick.major.width": .6, "ytick.major.width": .6})
fig, ax = plt.subplots(figsize=(8, 3.4)); fig.subplots_adjust(left=.08, right=.99, top=.88, bottom=.2)
x = np.arange(len(g)); z = g.z_height.values; has = ~np.isnan(z); zz = np.nan_to_num(z)
ylim = max(3, np.nanmax(np.abs(z)) * 1.25)
ax.bar(x[has], zz[has], width=.75, color=np.where(zz[has] > 0, UP, DN), zorder=3)
ax.bar(x[~has], np.full((~has).sum(), ylim * .12), bottom=-ylim * .06, width=.75, facecolor="none", edgecolor=GREY, lw=.5, ls=(0, (2, 1.5)), zorder=2)
mean, tot, n = np.nanmean(z), np.nansum(z), int(has.sum())
ax.axhline(0, color=INK, lw=.5, zorder=4); ax.axhline(mean, color=INK, lw=.9, ls="--", zorder=4)
ax.text(-.5, mean, f" mean z = {mean:+.2f}", va="bottom", ha="left", size=6.5, color=INK)
for v in (1.96, -1.96): ax.axhline(v, color=LIGHT, lw=.5, zorder=1)
ax.text(.005, .97, f"Height:  total z = {tot:+.2f}  ({n} of {len(g)} genes with a TWAS model; mean {mean:+.2f})",
        transform=ax.transAxes, size=7.5, weight="bold", va="top", color=INK)
ax.set_ylim(-ylim, ylim); ax.set_ylabel("TWAS z\n(mean over tissues)", size=7); ax.tick_params(labelsize=6, length=2)
ax.spines[["top", "right"]].set_visible(False); ax.set_xlim(-.7, len(g) - .3)
ax.set_xticks(x); ax.set_xticklabels(g.Gene, rotation=45, ha="right", rotation_mode="anchor", size=4)
for lab, hz in zip(ax.get_xticklabels(), has):
    if not hz: lab.set_color(GREY)
print("height total", round(tot, 3), "mean", round(mean, 3), "n", n, "of", len(g))
fig.suptitle("22q11.2 A–D: TWAS z of each gene (TWAS30)", size=9, weight="bold", y=.97)
fig.savefig(HERE + "/fig_22q11.2_A-D_twas_z.png", dpi=300); fig.savefig(HERE + "/fig_22q11.2_A-D_twas_z.pdf")
