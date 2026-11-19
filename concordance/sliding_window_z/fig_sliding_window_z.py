# Sliding windows of K = 10 adjacent TWAS genes vs random sets of 10 genes (MS 2026-11-07): are sums of adjacent genes' TWAS z
# inflated (neighbouring genes share eQTLs / LD, so their z are correlated) relative to sums of 10 genes drawn at random?
# Genes: TWAS30 per-gene z (mean over panels), ordered by chromosome and TxS; windows do not cross chromosomes.
# Random: 100,000 sets of 10 distinct genes drawn genome-wide (seed 2026). Variance ratio = var(adjacent sums) / var(random sums);
# under independence both are ~10 var(z). Right column: the same ratio for window sizes 1-40.
# Outputs: fig_sliding_window_z.{png,pdf}, sliding_window_z_summary.csv, sliding_window_z_by_k.csv, sliding_window_z_windows.csv
import os, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); CONC = os.path.dirname(HERE)
chrnum = lambda c: int(str(c).replace("chr", "").replace("X", "23"))
K = 10; NR = 100_000; rng = np.random.default_rng(2026)
def load(tr):
    tw = pd.read_csv(f"{CONC}/twas30/{tr}.txt", sep=r"\s+")
    g = tw.groupby("GENE").agg(z=("TWAS.Z", "mean"), chr=("CHR", "first"), pos=("TxS", "first")).reset_index()
    g["c"] = g["chr"].map(chrnum); return g.sort_values(["c", "pos"]).reset_index(drop=True)
def windows(g, k):   # sums of k consecutive genes within each chromosome; returns sums, start index
    out, st = [], []
    for c, d in g.groupby("c", sort=True):
        z = d.z.values
        if len(z) < k: continue
        cs = np.concatenate([[0], np.cumsum(z)]); out.append(cs[k:] - cs[:-k]); st.append(d.index.values[:len(z) - k + 1])
    return np.concatenate(out), np.concatenate(st)
def rand_sums(z, k, n):
    idx = np.argsort(rng.random((n, len(z))), axis=1)[:, :k] if len(z) * n < 5e7 else np.array([rng.choice(len(z), k, replace=False) for _ in range(n)])
    return z[idx].sum(1)
G = {tr: load(tr) for tr in ["height", "bmi"]}; summ, byk, wins = [], [], []
for tr, g in G.items():
    z = g.z.values; adj, st = windows(g, K); rnd = rand_sums(z, K, NR)
    same_chr = g.c.values[1:] == g.c.values[:-1]; lag1 = np.corrcoef(z[:-1][same_chr], z[1:][same_chr])[0, 1]
    summ.append(dict(trait=tr, n_genes=len(z), n_windows=len(adj), sd_gene_z=z.std(), sd_adjacent=adj.std(), sd_random=rnd.std(),
                     var_ratio=adj.var() / rnd.var(), mean_adjacent=adj.mean(), mean_random=rnd.mean(), lag1_corr_adjacent_genes=lag1,
                     frac_adjacent_abs_gt_random_p975=np.mean(np.abs(adj) > np.quantile(np.abs(rnd), .975))))
    for k in range(1, 41):
        a, _ = windows(g, k); r = rand_sums(z, k, 20_000)
        byk.append(dict(trait=tr, k=k, var_adjacent=a.var(), var_random=r.var(), var_ratio=a.var() / r.var()))
    wins.append(pd.DataFrame(dict(trait=tr, start_gene=g.GENE.values[st], chr=g.c.values[st], pos=g.pos.values[st], sum_z=adj)))
S = pd.DataFrame(summ); B = pd.DataFrame(byk); W = pd.concat(wins)
S.to_csv(HERE + "/sliding_window_z_summary.csv", index=False); B.to_csv(HERE + "/sliding_window_z_by_k.csv", index=False); W.to_csv(HERE + "/sliding_window_z_windows.csv", index=False)
print(S.round(3).T.to_string()); print(B[B.k.isin([1, 2, 5, 10, 20, 40])].round(2).to_string(index=False))

INK = '#1a1a1a'; GREY = '#8a8a8a'; LIGHT = '#dcdcdc'; BLUE = '#1f5c8b'; ORANGE = '#d98a3a'
plt.rcParams.update({"font.family": "DejaVu Sans", "axes.linewidth": .6, "xtick.major.width": .6, "ytick.major.width": .6})
fig = plt.figure(figsize=(9, 5)); gs = fig.add_gridspec(2, 3, width_ratios=[2.4, 1, 1], wspace=.32, hspace=.5)
fig.subplots_adjust(left=.07, right=.98, top=.9, bottom=.13)
for i, tr in enumerate(["height", "bmi"]):
    g = G[tr]; adj, st = windows(g, K); rnd = rand_sums(g.z.values, K, NR); s = S[S.trait == tr].iloc[0]; T = "Height" if tr == "height" else "BMI"
    ax = fig.add_subplot(gs[i, 0]); x = np.arange(len(adj)); cc = g.c.values[st]
    lo, hi = np.quantile(rnd, [.025, .975]); ax.axhspan(lo, hi, color=LIGHT, lw=0, zorder=0)
    ax.scatter(x, adj, s=.6, c=np.where(cc % 2, BLUE, "#7fa6c8"), lw=0, rasterized=True, zorder=2)
    mids = [x[cc == c].mean() for c in np.unique(cc)]; ax.set_xticks(mids); ax.set_xticklabels([str(c) if c < 23 else "X" for c in np.unique(cc)], size=4.6)
    ax.axhline(0, color=GREY, lw=.4); ax.set_xlim(-20, len(adj) + 20); ax.set_ylabel(f"sum of {K} adjacent z", size=6.5)
    ax.set_title(f"{T}: sliding windows of {K} adjacent TWAS genes", size=7.5, weight="bold", loc="left"); ax.tick_params(labelsize=5.5, length=2)
    ax.spines[["top", "right"]].set_visible(False); ax.set_xlabel("chromosome (window position)" if i == 1 else "", size=6)
    ax2 = fig.add_subplot(gs[i, 1]); lim = np.quantile(np.abs(np.concatenate([adj, rnd])), .995); bins = np.linspace(-lim, lim, 60)
    ax2.hist(rnd, bins=bins, density=True, color=GREY, alpha=.5, lw=0, label="random 10 genes")
    ax2.hist(adj, bins=bins, density=True, histtype="step", color=BLUE, lw=1.1, label=f"{K} adjacent genes")
    ax2.text(.03, .97, f"SD adjacent {s.sd_adjacent:.1f}\nSD random {s.sd_random:.1f}\nvariance ratio {s.var_ratio:.2f}", transform=ax2.transAxes, size=5, va="top", color=INK)
    ax2.set_ylim(0, ax2.get_ylim()[1] * 1.55); ax2.set_title("distribution of sums", size=7, weight="bold"); ax2.tick_params(labelsize=5.5, length=2); ax2.set_yticks([])
    ax2.spines[["top", "right", "left"]].set_visible(False); ax2.set_xlabel(f"sum of {K} z", size=6)
    if i == 0: ax2.legend(frameon=False, fontsize=5, loc="upper right", handlelength=1.2)
    ax3 = fig.add_subplot(gs[i, 2]); b = B[B.trait == tr]
    ax3.plot(b.k, b.var_ratio, color=BLUE, lw=1.1); ax3.axhline(1, color=GREY, lw=.6, ls="--"); ax3.axvline(K, color=LIGHT, lw=.6)
    ax3.set_title("adjacent / random variance", size=7, weight="bold"); ax3.set_xlabel("genes per window", size=6); ax3.set_ylim(0.8, max(b.var_ratio.max() * 1.1, 1.5))
    ax3.text(.97, .05, f"lag-1 corr of\nneighbouring genes {s.lag1_corr_adjacent_genes:+.2f}", transform=ax3.transAxes, ha="right", size=5.3, color=INK)
    ax3.tick_params(labelsize=5.5, length=2); ax3.spines[["top", "right"]].set_visible(False)
fig.suptitle("Do adjacent genes' TWAS z add up more than random genes? (TWAS30)", size=8.5, weight="bold", y=.98)
fig.text(.07, .01, f"Left: sum of z for each run of {K} consecutive TWAS-modelled genes (windows within chromosomes); grey band, central 95% of sums of {K} randomly drawn genes. "
         "Middle: adjacent windows (line) vs 100,000 random sets (grey). Right: var(adjacent sums) / var(random sums) by window size; 1 = no excess from adjacency.", size=5, color=GREY, va="bottom", wrap=True)
fig.savefig(HERE + "/fig_sliding_window_z.png", dpi=300); fig.savefig(HERE + "/fig_sliding_window_z.pdf")
