# Which TWAS30 tissues drive the Fig panel A correlation (height; primary non-overlapping set, all CNVs without Williams)?
# Reuses the exact pipeline from twas_dosage_panel_stats.py / twas_nonoverlapping.py, swapping the per-gene z.
import numpy as np, pandas as pd, sys
from scipy import stats
HERE = "/Users/mollysacks/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly/CNV_ASD_growth_analysis/concordance"
ROOT = "/Users/mollysacks/sebatlab Dropbox/Manuscripts/PGC_joint_submission/Molly"
genes = pd.read_csv(ROOT + "/pathway_analysis/Final/cnv_genes_annotated.csv")
S = pd.read_csv(HERE + "/twas_per_cnv_sign_test.csv"); S = S[S.trait == "height"].set_index("CNV_ID")
chrnum = lambda c: int(str(c).replace("chr", "").replace("X", "23"))
tw = pd.read_csv(HERE + "/twas30/height.txt", sep=r"\s+")
g = tw.groupby("GENE").agg(z=("TWAS.Z", "mean"), chr=("CHR", "first"), pos=("TxS", "first")).reset_index()
g["c"] = g["chr"].map(chrnum); g = g.sort_values(["c", "pos"]).reset_index(drop=True)
N = len(g); pos_of = dict(zip(g.GENE, g.index))
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
s0 = np.array([spans[c][0] for c in cids])
y = np.array([(S.loc[c, "DUP_effect"] - S.loc[c, "DEL_effect"]) / 2 for c in cids]); cls = np.array([S.loc[c, "cls"] for c in cids])
rng = np.random.default_rng(2026); K = rng.integers(0, N, size=(10000, len(clusters)))
st = (anc[None, :] + K[:, ci] + rel[None, :]) % N
DROP = ["16p11.2_BP2_BP3", "1q21.1_TAR_BP1_BP2", "1q21.1_TAR_BP2_BP3", "16p13.11_BP1_BP2", "16p13.11_BP2_BP3",
        "22q11.21_A_B", "22q11.21_B_C", "22q11.21_C_D", "22q11.21_B-D"]
C = np.array(cids); keep = ~np.isin(C, DROP)
kAll = keep & (C != "7q11.23_WBS"); kMir = keep & (cls == "mirror"); kAsy = keep & (cls == "asymmetric") & (C != "7q11.23_WBS")

def run(zvec):
    cs = np.concatenate([[0], np.cumsum(np.concatenate([zvec, zvec]))])
    x = cs[s0 + Ls] - cs[s0]; Xp = cs[st + Ls[None, :]] - cs[st]; out = {}
    for lab, k in [("all_noW", kAll), ("mirror", kMir), ("asym_noW", kAsy)]:
        r = stats.spearmanr(x[k], y[k])[0]; yr = stats.rankdata(y[k]); yr = (yr - yr.mean()) / yr.std()
        Xr = np.apply_along_axis(stats.rankdata, 1, Xp[:, k]); sd = Xr.std(1, keepdims=True); sd[sd == 0] = np.nan
        Xr = (Xr - Xr.mean(1, keepdims=True)) / sd; null = np.nan_to_num((Xr * yr[None, :]).mean(1), nan=-9)
        out[f"rho_{lab}"] = r; out[f"p_{lab}"] = (1 + np.sum(null >= r)) / 10001
    return out, x

GROUP = lambda s: ("brain" if s == "CMC" or "Brain" in s or s == "GTEX_Pituitary" else
                   "blood" if s in ("YFS", "NTR", "GTEX_Whole_Blood") else
                   "adipose" if s == "MET" or "Adipose" in s else "other_GTEx")
tw["grp"] = tw.STUDY.map(GROUP)
W = tw.pivot_table(index="GENE", columns="STUDY", values="TWAS.Z", aggfunc="mean").reindex(g.GENE)
def zmean(cols):  # per-gene mean over the given panels, 0 where the gene has none of them
    return W[cols].mean(1).fillna(0).values

rows = []
base, xbase = run(g.z.values); rows.append(dict(analysis="baseline (all panels)", **base))
for grp in ["brain", "blood", "adipose", "other_GTEx"]:
    cols = sorted(tw.loc[tw.grp == grp, "STUDY"].unique()); other = [c for c in W.columns if c not in cols]
    rows.append(dict(analysis=f"only {grp} ({len(cols)} panels)", **run(zmean(cols))[0]))
    rows.append(dict(analysis=f"drop {grp}", **run(zmean(other))[0]))
for t in W.columns:
    rows.append(dict(analysis=f"only {t}", n_genes_in_cnvs=None, **run(zmean([t]))[0]))
    rows.append(dict(analysis=f"drop {t}", **run(zmean([c for c in W.columns if c != t]))[0]))
R = pd.DataFrame(rows)

# Additive decomposition of the baseline x: x_cnv = sum_t sum_{genes in span} z_gt / n_g.
# Contribution of panel t to cov(x, y) over the panel-A CNVs (Pearson) and to cov(rank-agnostic) via y ranks.
n_g = W.notna().sum(1).values
ctrib = {}
for t in W.columns:
    v = (W[t].fillna(0).values / np.where(n_g > 0, n_g, 1)); cs = np.concatenate([[0], np.cumsum(v)])
    ctrib[t] = cs[s0 + Ls] - cs[s0]
Ct = pd.DataFrame(ctrib, index=C)
k = kAll; yc = y[k] - y[k].mean(); yrk = stats.rankdata(y[k]); yrk = yrk - yrk.mean()
dec = pd.DataFrame({"cov_y": (Ct[k].sub(Ct[k].mean())).mul(yc, axis=0).sum(), "cov_rank_y": (Ct[k].sub(Ct[k].mean())).mul(yrk, axis=0).sum(),
                    "n_gene_models_in_cnvs": [int(sum(W[t].notna().values[s:s + L].sum() for s, L in zip(s0[k], Ls[k]))) for t in Ct.columns]})
dec["share_cov"] = dec.cov_y / dec.cov_y.sum(); dec["share_cov_rank"] = dec.cov_rank_y / dec.cov_rank_y.sum()
dec["group"] = [GROUP(t) for t in dec.index]
counts = R.copy()
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 200)
print(R.round(4).to_string(index=False))
print(dec.sort_values("share_cov", ascending=False).round(3).to_string())
print(dec.groupby("group")[["share_cov", "share_cov_rank", "n_gene_models_in_cnvs"]].sum().round(3))
R.to_csv(sys.argv[1] + "/tissue_drivers_rho.csv", index=False); dec.to_csv(sys.argv[1] + "/tissue_drivers_cov_decomp.csv")
Ct.assign(dosage_effect=y, in_panel=kAll).to_csv(sys.argv[1] + "/tissue_contrib_per_cnv.csv")
