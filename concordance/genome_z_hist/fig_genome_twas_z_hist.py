# Genome-wide distribution of TWAS z (TWAS30, Mancuso 2017), height and BMI (MS 2026-11-07).
# Per-gene z = mean over the panels with a model (as used in panel A / the circular nulls). Grey curve: N(0, 1) scaled to the counts.
# Outputs: fig_genome_twas_z_hist.{png,pdf}, genome_twas_z_summary.csv
import os, numpy as np, pandas as pd, matplotlib
from scipy import stats
matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); CONC = os.path.dirname(HERE)
INK = '#1a1a1a'; GREY = '#8a8a8a'; BLUE = '#1f5c8b'
plt.rcParams.update({"font.family": "DejaVu Sans", "axes.linewidth": .6, "xtick.major.width": .6, "ytick.major.width": .6})
fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.8), gridspec_kw=dict(wspace=.25)); fig.subplots_adjust(left=.09, right=.98, top=.85, bottom=.2)
rows = []
for ax, tr in zip(axs, ["height", "bmi"]):
    tw = pd.read_csv(f"{CONC}/twas30/{tr}.txt", sep=r"\s+"); z = tw.groupby("GENE")["TWAS.Z"].mean().values
    rows.append(dict(trait=tr, n_genes=len(z), n_gene_tissue=len(tw), mean=z.mean(), median=np.median(z), sd=z.std(), min=z.min(), max=z.max(),
                     frac_abs_gt_1_96=np.mean(abs(z) > 1.96), frac_abs_gt_5=np.mean(abs(z) > 5), raw_gene_tissue_sd=tw["TWAS.Z"].std()))
    lim = 10; bins = np.arange(-lim, lim + .25, .25); zc = np.clip(z, -lim, lim)
    ax.hist(zc, bins=bins, color=BLUE, alpha=.85, lw=0)
    xs = np.linspace(-lim, lim, 400); ax.plot(xs, stats.norm.pdf(xs) * len(z) * .25, color=GREY, lw=.9, ls="--")
    r = rows[-1]
    ax.text(.98, .97, f"{r['n_genes']:,} genes\nmean {r['mean']:+.2f}, SD {r['sd']:.2f}\n|z| > 1.96: {100*r['frac_abs_gt_1_96']:.0f}%\n|z| > 5: {100*r['frac_abs_gt_5']:.1f}%\nrange {r['min']:.0f} to {r['max']:.0f}",
            transform=ax.transAxes, ha="right", va="top", size=5.8, color=INK, linespacing=1.3)
    ax.set_title("Height" if tr == "height" else "BMI", size=8, weight="bold"); ax.set_xlabel("TWAS z (mean over tissues)", size=7)
    ax.set_xlim(-lim, lim); ax.tick_params(labelsize=6, length=2); ax.spines[["top", "right"]].set_visible(False)
axs[0].set_ylabel("genes", size=7)
fig.suptitle("Genome-wide TWAS z distribution (TWAS30)", size=8.5, weight="bold", y=.97)
fig.text(.09, .01, f"Per-gene z averaged over the expression panels with a model. Dashed: N(0, 1) scaled to the gene count. Values beyond ±{lim} piled into the end bins.", size=5.2, color=GREY, va="bottom")
fig.savefig(HERE + "/fig_genome_twas_z_hist.png", dpi=300); fig.savefig(HERE + "/fig_genome_twas_z_hist.pdf")
R = pd.DataFrame(rows); R.to_csv(HERE + "/genome_twas_z_summary.csv", index=False); print(R.round(3).to_string(index=False))
