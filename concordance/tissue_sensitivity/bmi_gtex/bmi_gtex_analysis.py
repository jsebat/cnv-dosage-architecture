# BMI version of panel A using GTEx tissues (+ optionally PsychENCODE cortex) (MS 2026-11-03).
# Per CNV, BMI dosage effect (DUP - DEL)/2 against the summed BMI TWAS z of its genes, each gene's z = mean over the panels used.
#   mode gtex     : GTEx v6 panels of TWAS30 (Mancuso 2017; GIANT BMI GWAS). YFS, NTR, CMC, METSIM excluded.
#   mode gtex_pec : the same GTEx panels plus PsychENCODE cortex (PECfull FUSION weights, Tan et al. 2024 BMI GWAS; Jon's
#                   psychencode_twas/raw/BMI_TAN, as in twas_psychencode.py). NB the PsychENCODE z comes from a different BMI GWAS.
# Genes keyed by current HGNC symbol (TWAS30 / CNV gene names mapped through previous and alias symbols; PsychENCODE Ensembl IDs
# through hgnc_map.tsv). Gene universe = genes with >= 1 model in the panels used, ordered by hg19 position; CNV = contiguous run
# of those genes. Same non-overlapping CNV set as twas_dosage_bmi.py (spanning intervals kept, 22q11.2 A-D); all CNVs mirror for BMI.
# Spearman with circular genomic permutation P (10,000, seed 2026; locus clusters moved as units).
# Sensitivity: z from one panel only, and mean over all panels except one (genes lacking the panel(s) -> 0).
# Usage: python bmi_gtex_analysis.py [gtex|gtex_pec]   Outputs: bmi_<mode>_{points,stats,tissue_rho,gene_counts}.csv
import os, sys, glob, numpy as np, pandas as pd
from scipy import stats
MODE = sys.argv[1] if len(sys.argv) > 1 else "gtex"
HERE = os.path.dirname(os.path.abspath(__file__))
CONC = os.path.dirname(os.path.dirname(HERE)); ROOT = os.path.dirname(os.path.dirname(CONC))
DROP = ["16p11.2_BP2_BP3", "1q21.1_TAR_BP1_BP2", "1q21.1_TAR_BP2_BP3", "16p13.11_BP1_BP2", "16p13.11_BP2_BP3",
        "22q11.21_A_B", "22q11.21_B_C", "22q11.21_C_D", "22q11.21_B-D"]
hg = pd.read_csv(CONC + "/hgnc_map.tsv", sep="\t"); cur = set(hg.symbol); ens2sym = dict(zip(hg.ensembl_gene_id, hg.symbol)); alias = {}
for col in ("prev_symbol", "alias_symbol"):          # previous symbols take priority over aliases
    for s, v in zip(hg.symbol, hg[col]):
        if isinstance(v, str):
            for a in v.split("|"): alias.setdefault(a, s)
canon = lambda s: s if s in cur else alias.get(s, s)
genes = pd.read_csv(ROOT + "/pathway_analysis/Final/cnv_genes_annotated.csv"); genes["sym"] = genes.Gene.map(canon)
S = pd.read_csv(CONC + "/twas_per_cnv_sign_test.csv"); S = S[S.trait == "bmi"].set_index("CNV_ID")
chrnum = lambda c: int(str(c).replace("chr", "").replace("X", "23"))
tw = pd.read_csv(CONC + "/twas30/bmi.txt", sep=r"\s+"); tw = tw[tw.STUDY.str.startswith("GTEX_")].copy()
tw = pd.DataFrame(dict(sym=tw.GENE.map(canon), panel=tw.STUDY, z=tw["TWAS.Z"], chr=tw.CHR.map(chrnum), pos=tw.TxS))
if MODE == "gtex_pec":
    pe = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(glob.glob(CONC + "/psychencode_twas/raw/BMI_TAN/PECfull.*.dat"))])
    for c in ("TWAS.Z", "P0", "CHR"): pe[c] = pd.to_numeric(pe[c], errors="coerce")
    pe = pe.dropna(subset=["TWAS.Z", "P0", "CHR"]); pe["sym"] = pe.ID.str.split(".").str[0].map(ens2sym); pe = pe.dropna(subset=["sym"])
    tw = pd.concat([tw, pd.DataFrame(dict(sym=pe.sym, panel="PsychENCODE", z=pe["TWAS.Z"], chr=pe.CHR.astype(int), pos=pe.P0))], ignore_index=True)
g = tw.groupby("sym").agg(z=("z", "mean"), c=("chr", "first"), pos=("pos", "first")).reset_index()
g = g.sort_values(["c", "pos"]).reset_index(drop=True)
N = len(g); pos_of = dict(zip(g.sym, g.index))
spans = {}
for cid, sub in genes.groupby("CNV_ID"):
    if cid in DROP: continue
    ix = sorted({pos_of[x] for x in sub.sym if x in pos_of})
    if ix and cid in S.index and not np.isnan(S.loc[cid, "DEL_effect"]) and not np.isnan(S.loc[cid, "DUP_effect"]): spans[cid] = (ix[0], ix[-1] - ix[0] + 1)
cids = list(spans); order = sorted(cids, key=lambda c: spans[c][0]); clusters, cur_, end = [], [], -1
for c in order:
    s0_, L = spans[c]
    if cur_ and s0_ <= end: cur_.append(c); end = max(end, s0_ + L - 1)
    else:
        if cur_: clusters.append(cur_)
        cur_, end = [c], s0_ + L - 1
clusters.append(cur_)
anchor = {c: min(spans[x][0] for x in cl) for cl in clusters for c in cl}; cl_id = {c: i for i, cl in enumerate(clusters) for c in cl}
rel = np.array([spans[c][0] - anchor[c] for c in cids]); Ls = np.array([spans[c][1] for c in cids]); ci = np.array([cl_id[c] for c in cids]); anc = np.array([anchor[c] for c in cids])
s0 = np.array([spans[c][0] for c in cids]); C = np.array(cids)
y = np.array([(S.loc[c, "DUP_effect"] - S.loc[c, "DEL_effect"]) / 2 for c in cids])
rng = np.random.default_rng(2026); K = rng.integers(0, N, size=(10000, len(clusters)))
st = (anc[None, :] + K[:, ci] + rel[None, :]) % N
kAll = np.ones(len(C), bool); kNoW = C != "7q11.23_WBS"

def run(zvec, sets=(("all", kAll), ("noW", kNoW))):
    cs = np.concatenate([[0], np.cumsum(np.concatenate([zvec, zvec]))])
    x = cs[s0 + Ls] - cs[s0]; Xp = cs[st + Ls[None, :]] - cs[st]; out = {}
    for lab, k in sets:
        if np.all(x[k] == x[k][0]): out.update({f"rho_{lab}": np.nan, f"p1_{lab}": np.nan, f"p2_{lab}": np.nan}); continue
        r = stats.spearmanr(x[k], y[k])[0]; yr = stats.rankdata(y[k]); yr = (yr - yr.mean()) / yr.std()
        Xr = np.apply_along_axis(stats.rankdata, 1, Xp[:, k]); sd = Xr.std(1, keepdims=True); sd[sd == 0] = np.nan
        null = np.nan_to_num(((Xr - Xr.mean(1, keepdims=True)) / sd * yr[None, :]).mean(1), nan=0)
        out.update({f"rho_{lab}": r, f"p1_{lab}": (1 + np.sum(null >= r)) / 10001, f"p2_{lab}": (1 + np.sum(np.abs(null) >= abs(r))) / 10001})
    return out, x

W = tw.pivot_table(index="sym", columns="panel", values="z", aggfunc="mean").reindex(g.sym)
zmean = lambda cols: W[cols].mean(1).fillna(0).values
base, xb = run(g.z.values)
ts = stats.theilslopes(y, xb); tsW = stats.theilslopes(y[kNoW], xb[kNoW])
xpec = None
if MODE == "gtex_pec": xpec = run(zmean(["PsychENCODE"]))[1]; xgtex = run(zmean([c for c in W.columns if c != "PsychENCODE"]))[1]
pts = pd.DataFrame(dict(CNV_ID=C, summed_bmi_twas_z=xb, DEL=[S.loc[c, "DEL_effect"] for c in C], DUP=[S.loc[c, "DUP_effect"] for c in C],
                        dosage_effect_bmi=y, n_genes_twas=Ls))
if xpec is not None: pts["summed_z_psychencode_only"] = xpec; pts["summed_z_gtex_only"] = xgtex
pre = f"{HERE}/bmi_{MODE}_"; pts.to_csv(pre + "points.csv", index=False)
pd.DataFrame([dict(group="all", n_cnv=int(kAll.sum()), spearman_rho=base["rho_all"], p_one_sided=base["p1_all"], p_two_sided=base["p2_all"],
                   theilsen_slope=ts[0], theilsen_intercept=ts[1]),
              dict(group="all without Williams", n_cnv=int(kNoW.sum()), spearman_rho=base["rho_noW"], p_one_sided=base["p1_noW"], p_two_sided=base["p2_noW"],
                   theilsen_slope=tsW[0], theilsen_intercept=tsW[1])]).to_csv(pre + "stats.csv", index=False)
rows = [dict(analysis="all", tissue="ALL", **base)]
for t in W.columns:
    rows.append(dict(analysis="only", tissue=t, **run(zmean([t]))[0]))
    rows.append(dict(analysis="drop", tissue=t, **run(zmean([c for c in W.columns if c != t]))[0]))
R = pd.DataFrame(rows); R.to_csv(pre + "tissue_rho.csv", index=False)
idx = np.unique(np.concatenate([np.arange(a, a + L) for a, L in zip(s0, Ls)])); H = W.iloc[idx].notna()
GC = pd.DataFrame({"n_genes_single": H.sum(), "n_genes_all_but": [int(H.drop(columns=t).any(axis=1).sum()) for t in H.columns]})
GC.loc["ALL"] = [int(H.any(axis=1).sum()), np.nan]; GC.to_csv(pre + "gene_counts.csv")
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 200)
print(MODE, "genes in universe:", N, "CNVs:", len(C))
print(pd.read_csv(pre + "stats.csv").round(4).to_string(index=False))
print(R.round(3).sort_values(["analysis", "p1_all"]).to_string(index=False))
