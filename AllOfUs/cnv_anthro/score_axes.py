"""Person-level fitted scores for the CNV gene-dosage axes (mirror, constraint, burden) and what they explain.

A design variant (as in conga_grid_qt.py: DEL terms + DUP terms) gives each person counts summed over the genes hit by ALL of
their CNVs of that dosage (a gene hit twice in one dosage counts once). Terms are grouped into axes:
  burden     n_all                                   (total genes hit)
  mirror     g:... TWAS measures (any call-size prefix)
  constraint none_H, n95, pres:...                   (constrained-gene measures)
Fit mode (default): K-fold cross-fitting. OLS of the covariate-residualised trait on the residualised counts in K-1 folds gives
weights; each held-out person gets axis scores = sum over that axis's columns of weight x count. Out-of-fold scores are then
tested: trait ~ covariates + score (per axis alone, and all axes jointly), with beta per SD of score, P and partial R2 (share
of trait variance explained). Additivity: among carriers, does the score slope differ between carriers of one gene-hitting CNV
and of two or more (interaction test)?
Full-sample weights are written (aggregate) so other cohorts can apply them: --weights <file> scores with fixed weights
(no refit) and reports the same statistics, i.e. an out-of-sample test.

  python score_axes.py --pheno P --cnv C --trait height_rint --twas T --design designs/joint_2026-11-15.json --variant height_selected \\
      --exons E --gene_cov G --clingen CG --out OUT [--folds 5] [--weights other_cohort_weights.tsv]
"""
import argparse
import json
import os
import numpy as np
import pandas as pd
from scipy import stats
from common import read_plink_cnv, design_matrix
from conga_axes_qt import load_genes, incidence
from conga_grid_qt import twas_z_both, measure


def axis_of(term):
    base = term.split("~")[-1]
    if base == "n_all":
        return "burden"
    if base.startswith("g:"):
        return "mirror"
    return "constraint"


def ols_score_test(yr, s, Q, df_extra=0):
    """y_resid ~ s (s residualised on covariates); returns beta per SD, SE, P, partial R2."""
    sr = s - Q @ (Q.T @ s)
    ss = sr @ sr
    if ss <= 0:
        return dict(beta_per_SD=np.nan, SE=np.nan, P=np.nan, partial_R2=np.nan)
    b = (sr @ yr) / ss; n = len(yr); dfr = n - Q.shape[1] - 1 - df_extra
    rss = yr @ yr - b * (sr @ yr); se = np.sqrt(rss / dfr / ss); t = b / se
    sd = np.sqrt(ss / n)
    return dict(beta_per_SD=b * sd, SE=se * sd, T=t, P=2 * stats.t.sf(abs(t), dfr), partial_R2=t ** 2 / (t ** 2 + dfr))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pheno", required=True); ap.add_argument("--trait", required=True); ap.add_argument("--cnv", required=True)
    ap.add_argument("--twas", required=True); ap.add_argument("--panels", default="gtex_nonbrain")
    ap.add_argument("--exons", required=True); ap.add_argument("--gene_cov", required=True); ap.add_argument("--clingen", required=True)
    ap.add_argument("--design", required=True); ap.add_argument("--variant", required=True)
    ap.add_argument("--folds", type=int, default=5); ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--weights", default=None, help="apply these fixed weights (term, dosage, weight) instead of fitting")
    ap.add_argument("--numeric", default="age,age2,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10")
    ap.add_argument("--factors", default="sex,batch"); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    V = [v for v in json.load(open(a.design)) if v["name"] == a.variant][0]
    ph = pd.read_csv(a.pheno, sep="\t", dtype={"IID": str})
    num = [c for c in a.numeric.split(",") if c]; fac = [c for c in a.factors.split(",") if c]
    ph = ph[np.isfinite(ph[a.trait])].dropna(subset=num + fac).reset_index(drop=True)
    Zc, _ = design_matrix(ph, num, fac); Q, _ = np.linalg.qr(Zc); N = len(ph)
    y = ph[a.trait].to_numpy(float); yr = y - Q @ (Q.T @ y)
    cnv = read_plink_cnv(a.cnv, keep_iids=ph.IID); cnv = cnv[cnv.CHR <= 22].reset_index(drop=True)
    cnv["row"] = pd.Series(np.arange(N), index=ph.IID).loc[cnv.IID].to_numpy()
    U, ex, tx = load_genes(a.exons); M = len(U)
    gc = pd.read_csv(a.gene_cov, sep="\t", dtype={"entrez": str}); gc.loc[gc.shet_imputed == 1, "shet"] = np.nan
    shet = U[["entrez"]].merge(gc[["entrez", "shet"]], on="entrez", how="left").shet.to_numpy()
    q95 = float(np.nanquantile(gc.shet.to_numpy(float), 0.95))
    cg = pd.read_csv(a.clingen, sep="\t", skiprows=5, dtype=str, quoting=3); cg.columns = [c.strip().lstrip("#") for c in cg.columns]
    hi = U.symbol.isin(set(cg.loc[cg["Haploinsufficiency Score"].isin(["2", "3"]), cg.columns[0]])).to_numpy()
    Z, _ = twas_z_both(a.twas, a.panels, U, gc)
    # per-call gene counts (for call-size classes and for the number of gene-hitting CNVs per person)
    cc_ = cnv.copy(); cc_["row"] = np.arange(len(cc_)); ng = np.zeros(len(cnv))
    for typ, feats, mode in [(1, ex, "overlap"), (3, tx, "within")]:
        m = (cnv.TYPE == typ).to_numpy(); ng[m] = incidence(cc_[m], feats, mode, len(cnv), M).getnnz(axis=1)[m]
    cnv["ng"] = ng
    n_calls = np.bincount(cnv.row[cnv.ng > 0], minlength=N)
    G = {}
    def G_for(dos, cls):
        key = (dos, cls)
        if key not in G:
            c_ = cnv[cnv.TYPE == (1 if dos == "DEL" else 3)]
            if cls:
                kk = int(cls[1:]); c_ = c_[(c_.ng > kk) if cls[0] == "L" else (c_.ng <= kk)]
            G[key] = incidence(c_, ex if dos == "DEL" else tx, "overlap" if dos == "DEL" else "within", N, M)
        return G[key]
    cols, names = [], []
    for dos in ("DEL", "DUP"):
        for t in V[dos]:
            cls, base = (t.split("~")[0], t.split("~")[1]) if "~" in t else ("", t)
            f = np.asarray(G_for(dos, cls) @ measure(base, Z, shet, q95, hi)).ravel()
            if base.startswith("pres:"):
                f = (f > 0).astype(float)
            cols.append(f); names.append((t, dos))
    F = np.column_stack(cols); Fr = F - Q @ (Q.T @ F)
    axes = np.array([axis_of(t) for t, _ in names])
    info = dict(trait=a.trait, N=int(N), variant=a.variant, terms=[f"{t}|{d}" for t, d in names],
                carriers_with_gene_hit=int((n_calls > 0).sum()), carriers_2plus=int((n_calls >= 2).sum()))
    # full-sample weights (exportable)
    A = Fr.T @ Fr; b_full = np.linalg.solve(A, Fr.T @ yr)
    s2 = (yr @ yr - b_full @ (Fr.T @ yr)) / (N - Q.shape[1] - F.shape[1]); se_full = np.sqrt(np.diag(np.linalg.inv(A)) * s2)
    W = pd.DataFrame({"term": [t for t, _ in names], "dosage": [d for _, d in names], "axis": axes, "weight": b_full, "SE": se_full,
                      "mean_count": F.mean(0), "sd_count": F.std(0)})
    W.to_csv(os.path.join(a.out, "weights_fullsample.tsv"), sep="\t", index=False, float_format="%.6g")
    # scores: fixed weights (another cohort's) or K-fold cross-fitted
    if a.weights:
        w = pd.read_csv(a.weights, sep="\t").set_index(["term", "dosage"]).weight
        beta = np.array([w.get((t, d), 0.0) for t, d in names]); S = F * beta; mode = f"fixed weights from {os.path.basename(a.weights)}"
    else:
        rng = np.random.default_rng(a.seed); fold = rng.integers(0, a.folds, N); S = np.zeros_like(F)
        for k in range(a.folds):
            tr = fold != k; te = ~tr
            bk = np.linalg.solve(Fr[tr].T @ Fr[tr], Fr[tr].T @ yr[tr]); S[te] = F[te] * bk
        mode = f"{a.folds}-fold cross-fitted"
    info["mode"] = mode
    scores = {ax: S[:, axes == ax].sum(1) for ax in ["mirror", "constraint", "burden"] if (axes == ax).any()}
    scores["total"] = S.sum(1)
    rows = []
    for ax, s in scores.items():
        rows.append(dict(axis=ax, model="alone", **ols_score_test(yr, s, Q), carriers_nonzero=int((s != 0).sum())))
    # joint: all axis scores together
    Xs = np.column_stack([scores[ax] for ax in scores if ax != "total"]); Xr = Xs - Q @ (Q.T @ Xs)
    bj = np.linalg.solve(Xr.T @ Xr, Xr.T @ yr); dfr = N - Q.shape[1] - Xs.shape[1]
    s2j = (yr @ yr - bj @ (Xr.T @ yr)) / dfr; sej = np.sqrt(np.diag(np.linalg.inv(Xr.T @ Xr)) * s2j); tj = bj / sej
    for j, ax in enumerate([ax for ax in scores if ax != "total"]):
        sd = Xs[:, j].std()
        rows.append(dict(axis=ax, model="joint (all axis scores)", beta_per_SD=bj[j] * sd, SE=sej[j] * sd, T=tj[j],
                         P=2 * stats.t.sf(abs(tj[j]), dfr), partial_R2=tj[j] ** 2 / (tj[j] ** 2 + dfr), carriers_nonzero=int((Xs[:, j] != 0).sum())))
    R = pd.DataFrame(rows)
    tot = R[(R.axis == "total") & (R.model == "alone")].partial_R2.iloc[0]
    R["share_of_total_R2"] = np.where(R.model == "alone", R.partial_R2 / tot, np.nan)
    R.to_csv(os.path.join(a.out, "axis_scores.tsv"), sep="\t", index=False, float_format="%.6g")
    # additivity: score slope in carriers of 1 vs >= 2 gene-hitting CNVs (interaction), per axis
    add = []
    car = n_calls > 0; multi = (n_calls >= 2).astype(float)
    for ax, s in scores.items():
        X = np.column_stack([s, multi, s * multi])[car]; Qc = Q[car]; yc = yr[car]
        Xr_ = X - Qc @ (Qc.T @ X); yc_ = yc - Qc @ (Qc.T @ yc)
        try:
            bb = np.linalg.solve(Xr_.T @ Xr_, Xr_.T @ yc_)
        except np.linalg.LinAlgError:
            continue
        dfc = car.sum() - Qc.shape[1] - 3; s2c = (yc_ @ yc_ - bb @ (Xr_.T @ yc_)) / dfc; sec = np.sqrt(np.diag(np.linalg.inv(Xr_.T @ Xr_)) * s2c)
        add.append(dict(axis=ax, slope_single=bb[0], SE_single=sec[0], slope_multi_minus_single=bb[2], SE_diff=sec[2], P_interaction=2 * stats.t.sf(abs(bb[2] / sec[2]), dfc),
                        n_single=int(((n_calls == 1)).sum()), n_multi=int((n_calls >= 2).sum())))
    pd.DataFrame(add).to_csv(os.path.join(a.out, "additivity.tsv"), sep="\t", index=False, float_format="%.6g")
    json.dump(info, open(os.path.join(a.out, "run_info.json"), "w"), indent=1)
    print(json.dumps(info, indent=1)); print(R.to_string()); print(pd.DataFrame(add).to_string())


if __name__ == "__main__":
    main()
