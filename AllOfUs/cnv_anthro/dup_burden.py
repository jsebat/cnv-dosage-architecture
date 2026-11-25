"""Genome-wide duplication burden and height: is the burden -> short stature association driven by monogenic-like duplications,
with duplicated gene count raising height once they are removed? (Jon 2026-11-15; Sacks Fig. 2D vs 2M.)

Per person: n genes duplicated (DUP spans the whole transcript; each gene once), n genes deleted (any exon), and whether any
duplication contains a 'hit' gene (top-5% GeneBayes s_het or ClinGen haploinsufficiency 2-3). Models (OLS, covariate-residualised):
  M1 all people:                          trait ~ n_DUP + n_DEL
  M2 asymmetric-interval carriers removed (Molly's 13 intervals, either dosage, Molly's recurrent CNV genotyping)
  M3 carriers of a hit-gene duplication removed
  M4 all people:                          trait ~ n_DUP + n_DEL + any hit-gene DUP
  M5 carriers of a hit-gene duplication only vs non-carriers: n_DUP slope within the hit group
Also the duplicated-gene-count dose-response by bin in the M3 sample. Aggregate outputs only.
"""
import argparse, json, os
import numpy as np, pandas as pd
from scipy import stats
from common import read_plink_cnv, design_matrix, rcg
from conga_axes_qt import load_genes, incidence
from additivity_constraint import qualifying_calls

BINS = [(0, 0, "0"), (1, 2, "1-2"), (3, 5, "3-5"), (6, 10, "6-10"), (11, 20, "11-20"), (21, 10 ** 6, "21+")]


def fit(yr, X, Q, names):
    # Rank-safe (as ols() in additivity_constraint.py): a column is kept only if it raises the rank of the residualised design, in
    # order; aliased columns (e.g. no duplicated gene in a small WGS stratum) are reported as NaN and df comes from the kept terms.
    # Identical to the old full-rank solve whenever the design is full rank.
    Xr = X - Q @ (Q.T @ X); keep = []
    for j in range(Xr.shape[1]):
        if np.linalg.matrix_rank(Xr[:, keep + [j]]) == len(keep) + 1:
            keep.append(j)
    out = {n: (np.nan, np.nan, np.nan) for n in names}
    df = len(yr) - Q.shape[1] - len(keep)
    if not keep or df < 1:
        return out
    K = Xr[:, keep]; A = K.T @ K; b = np.linalg.solve(A, K.T @ yr)
    se = np.sqrt(np.diag(np.linalg.inv(A)) * (yr @ yr - b @ (K.T @ yr)) / df)
    for i, j in enumerate(keep):
        out[names[j]] = (b[i], se[i], 2 * stats.t.sf(abs(b[i] / se[i]), df))
    return out


def fit_guarded(yr, X, Q, names, min_n):
    """fit() unless the group is too small (N < min_n) -> NaN rows; terms with no variation come back NaN from fit()."""
    if len(yr) < min_n:
        return {n: (np.nan, np.nan, np.nan) for n in names}
    return fit(yr, X, Q, names)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pheno", required=True); ap.add_argument("--trait", required=True); ap.add_argument("--cnv", required=True)
    ap.add_argument("--exons", required=True); ap.add_argument("--gene_cov", required=True); ap.add_argument("--clingen", required=True)
    ap.add_argument("--asym_loci", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--numeric", default="age,age2,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"); ap.add_argument("--factors", default="sex,batch")
    ap.add_argument("--min_n", type=int, default=50, help="model groups with fewer people are not fitted (NaN rows)")
    ap.add_argument("--single_def", choices=["all_calls", "rare_gene"], default="all_calls",
                    help="signed-burden panel: who counts as deletion-only / duplication-only and what is counted. all_calls = every call "
                         "(arrays, UKBB); rare_gene = rare gene-hitting calls only, the WGS reference rule of the 'only' runs (no_rare_gene)")
    ap.add_argument("--rare_freq", type=float, default=0.01)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    ph = pd.read_csv(a.pheno, sep="\t", dtype={"IID": str})
    num = [c for c in a.numeric.split(",") if c]; fac = [c for c in a.factors.split(",") if c]
    ph = ph[np.isfinite(ph[a.trait])].dropna(subset=num + fac).reset_index(drop=True); N = len(ph)
    cnv = read_plink_cnv(a.cnv, keep_iids=ph.IID); cnv = cnv[cnv.CHR <= 22].reset_index(drop=True)
    cnv["row"] = pd.Series(np.arange(N), index=ph.IID).loc[cnv.IID].to_numpy()
    U, ex, tx = load_genes(a.exons); M = len(U)
    gc = pd.read_csv(a.gene_cov, sep="\t", dtype={"entrez": str}); gc.loc[gc.shet_imputed == 1, "shet"] = np.nan
    shet = U[["entrez"]].merge(gc[["entrez", "shet"]], on="entrez", how="left").shet.to_numpy()
    q95 = float(np.nanquantile(gc.shet.to_numpy(float), 0.95))
    cg = pd.read_csv(a.clingen, sep="\t", skiprows=5, dtype=str, quoting=3); cg.columns = [c.strip().lstrip("#") for c in cg.columns]
    hi = U.symbol.isin(set(cg.loc[cg["Haploinsufficiency Score"].isin(["2", "3"]), cg.columns[0]])).to_numpy()
    hitw = ((np.nan_to_num(shet, nan=-1) >= q95) | hi).astype(float)
    Gd = incidence(cnv[cnv.TYPE == 1], ex, "overlap", N, M); Gu = incidence(cnv[cnv.TYPE == 3], tx, "within", N, M)
    nDUP = np.asarray(Gu.sum(1)).ravel(); nDEL = np.asarray(Gd.sum(1)).ravel(); hitDUP = (np.asarray(Gu @ hitw).ravel() > 0)
    hcDUP = (np.asarray(Gu @ (np.nan_to_num(shet, nan=-1) > 0.029).astype(float)).ravel() > 0)   # any HC (s_het > 0.029) gene duplicated
    asym = set()
    asym = rcg.carriers(rcg.genotype_calls(cnv), pd.read_csv(a.asym_loci, sep="\t").NAME, "both")   # Molly's genotyping (Jon 2026-11-24)
    asymm = ph.IID.isin(asym).to_numpy()
    Zc, _ = design_matrix(ph, num, fac); y = ph[a.trait].to_numpy(float)
    rows = []
    def run(name, keep, extra=None):
        Q, _ = np.linalg.qr(Zc[keep]) if keep.sum() else (np.zeros((0, Zc.shape[1])), None); yr = y[keep] - Q @ (Q.T @ y[keep])
        cols = [nDUP[keep] / 10, nDEL[keep] / 10]; names = ["DUP genes (per 10)", "DEL genes (per 10)"]
        if extra is not None:
            cols.append(extra[keep]); names.append("any hit-gene DUP")
        r = fit_guarded(yr, np.column_stack(cols), Q, names, a.min_n)
        for k, (b, se, p) in r.items():
            rows.append(dict(model=name, term=k, beta_SD=b, SE=se, P=p, N=int(keep.sum()), dup_carriers=int((nDUP[keep] > 0).sum())))
    allp = np.ones(N, bool)
    run("M1 all people", allp)
    run("M2 asymmetric-interval carriers removed", ~asymm)
    run("M3 hit-gene DUP carriers removed", ~hitDUP)
    run("M4 all people + any hit-gene DUP", allp, hitDUP.astype(float))
    run("M5 hit-gene DUP carriers + non-DUP carriers", hitDUP | (nDUP == 0))
    R = pd.DataFrame(rows); R.to_csv(os.path.join(a.out, "dup_burden_models.tsv"), sep="\t", index=False, float_format="%.6g")
    Q, _ = np.linalg.qr(Zc); yr = y - Q @ (Q.T @ y); b = []
    for grp, m in [("all", allp), ("hit-gene DUP carriers removed", ~hitDUP), ("hit-gene DUP carriers", hitDUP)]:
        for lo, hi_, lab in BINS:
            k = m & (nDUP >= lo) & (nDUP <= hi_)
            if k.sum() > 1:
                b.append(dict(group=grp, dup_genes=lab, n=int(k.sum()), mean_SD=yr[k].mean(), SE=yr[k].std(ddof=1) / np.sqrt(k.sum())))
    pd.DataFrame(b).to_csv(os.path.join(a.out, "dup_burden_bins.tsv"), sep="\t", index=False, float_format="%.6g")
    # ---- signed genome-wide burden (Jon 2026-11-15: U-shaped dose response caused by monogenic asymmetric effects)
    # people with deletions only, duplications only, or neither (carriers of both dosages left out of this display);
    # x = -(genes deleted) to the left, +(genes duplicated) to the right; same bins on each side
    if a.single_def == "rare_gene":
        q = qualifying_calls(cnv, N, a.rare_freq, a.exons); cq = cnv[q]
        sDEL = np.asarray(incidence(cq[cq.TYPE == 1], ex, "overlap", N, M).sum(1)).ravel()
        Gq = incidence(cq[cq.TYPE == 3], tx, "within", N, M); sDUP = np.asarray(Gq.sum(1)).ravel()
        sHIT = np.asarray(Gq @ hitw).ravel() > 0; sHC = np.asarray(Gq @ (np.nan_to_num(shet, nan=-1) > 0.029).astype(float)).ravel() > 0
    else:
        sDEL, sDUP, sHIT, sHC = nDEL, nDUP, hitDUP, hcDUP
    # signed-burden section only (restored below): counts, and the hit/HC-gene DUP removals, use the same calls
    saved = (nDEL, nDUP, hitDUP, hcDUP); nDEL, nDUP, hitDUP, hcDUP = sDEL, sDUP, sHIT, sHC
    single = (nDEL == 0) | (nDUP == 0)
    sb = []; sm = []
    for grp, m in [("all", single), ("hit-gene DUP carriers removed", single & ~hitDUP), ("asymmetric-interval carriers removed", single & ~asymm),
                   ("HC-gene DUP carriers removed", single & ~hcDUP)]:
        for side, cnt, sgn in [("DEL", nDEL, -1), ("DUP", nDUP, 1)]:
            for lo, hi_, lab in BINS[1:]:
                k = m & (cnt >= lo) & (cnt <= hi_) & ((nDUP == 0) if side == "DEL" else (nDEL == 0))
                if k.sum() > 1:
                    sb.append(dict(group=grp, bin=f"{side} {lab}", side=side, order=sgn * BINS.index((lo, hi_, lab)), n=int(k.sum()),
                                   mean_SD=yr[k].mean(), SE=yr[k].std(ddof=1) / np.sqrt(k.sum())))
        k = m & (nDEL == 0) & (nDUP == 0)
        if k.sum() > 1:
            sb.append(dict(group=grp, bin="none", side="none", order=0, n=int(k.sum()), mean_SD=yr[k].mean(), SE=yr[k].std(ddof=1) / np.sqrt(k.sum())))
        # piecewise slopes and the quadratic (U-shape) term of the signed burden s = (genes duplicated - genes deleted) / 10
        Qm, _ = np.linalg.qr(Zc[m]) if m.sum() else (np.zeros((0, Zc.shape[1])), None); ym = y[m] - Qm @ (Qm.T @ y[m])
        r = fit_guarded(ym, np.column_stack([nDEL[m] / 10, nDUP[m] / 10]), Qm, ["DEL genes (per 10)", "DUP genes (per 10)"], a.min_n)
        s = (nDUP[m] - nDEL[m]) / 10.0
        r.update(fit_guarded(ym, np.column_stack([s, s ** 2]), Qm, ["signed burden s", "s^2 (U-shape)"], a.min_n))
        for kk, (bb, se, p) in r.items():
            sm.append(dict(group=grp, term=kk, beta_SD=bb, SE=se, P=p, N=int(m.sum())))
    nDEL, nDUP, hitDUP, hcDUP = saved
    pd.DataFrame(sb, columns=["group", "bin", "side", "order", "n", "mean_SD", "SE"]).sort_values(["group", "order"]).to_csv(os.path.join(a.out, "signed_burden_bins.tsv"), sep="\t", index=False, float_format="%.6g")
    pd.DataFrame(sm).to_csv(os.path.join(a.out, "signed_burden_models.tsv"), sep="\t", index=False, float_format="%.6g")
    print(pd.DataFrame(sm).to_string())
    info = dict(trait=a.trait, N=N, single_def=a.single_def, min_n=a.min_n, single_dosage_N=int(single.sum()),
                rare_freq=a.rare_freq if a.single_def == "rare_gene" else None, dup_carriers=int((nDUP > 0).sum()), hit_dup_carriers=int(hitDUP.sum()), asym_carriers=int(asymm.sum()), hc_dup_carriers=int(hcDUP.sum()), q95=q95)
    json.dump(info, open(os.path.join(a.out, "run_info.json"), "w"), indent=1)
    pd.set_option("display.width", 200); print(json.dumps(info)); print(R.to_string()); print(pd.DataFrame(b).to_string())


if __name__ == "__main__":
    main()
