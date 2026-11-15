"""Shared I/O, trait transform and fast exact OLS for the CNV anthropometry module.

Portable: no cluster paths in here. Everything cohort-specific comes in through arguments.

OLS strategy (exact, not an approximation): with covariate matrix Z (N x k, incl. intercept) and its thin QR
factor Q, the coefficient of any extra regressor x in  y ~ Z + x  equals the regression of the Z-residualised y
on the Z-residualised x (Frisch-Waugh-Lovell). For a sparse 0/1 carrier vector only the carriers' rows of Q are
needed, so a probe test costs O(n_carriers * k) instead of O(N * k).
"""
import numpy as np
import pandas as pd
from scipy import stats

PLINK_CNV_COLS = ["FID", "IID", "CHR", "BP1", "BP2", "TYPE", "SCORE", "SITES"]


def read_plink_cnv(path, keep_iids=None):
    """PLINK .cnv with header (FID IID CHR BP1 BP2 TYPE SCORE SITES). TYPE 1 = DEL, 3 = DUP. CHR numeric (23 = X)."""
    d = pd.read_csv(path, sep=r"\s+", dtype={"FID": str, "IID": str}, keep_default_na=False)
    d = d[[c for c in PLINK_CNV_COLS if c in d.columns]].copy()
    for c in ["CHR", "BP1", "BP2", "TYPE"]:
        d[c] = pd.to_numeric(d[c]).astype(np.int64)
    if keep_iids is not None:
        d = d[d["IID"].isin(set(keep_iids))]
    d = d[d["TYPE"].isin([1, 3])].reset_index(drop=True)
    d["CNVID"] = np.arange(len(d))
    return d


def rint(x):
    """Rank-based inverse normal transform (Blom offset 3/8); ties get average ranks; NaN kept."""
    x = np.asarray(x, float)
    out = np.full(x.shape, np.nan)
    ok = np.isfinite(x)
    r = stats.rankdata(x[ok])
    n = ok.sum()
    out[ok] = stats.norm.ppf((r - 0.375) / (n + 0.25))
    return out


def rint_within(x, groups):
    x = np.asarray(x, float)
    out = np.full(x.shape, np.nan)
    g = np.asarray(groups)
    for v in pd.unique(g):
        m = g == v
        out[m] = rint(x[m])
    return out


def design_matrix(df, numeric, factors):
    """Intercept + numeric columns + treatment-coded factors (first level dropped). Returns (Z, names)."""
    cols = [np.ones(len(df))]
    names = ["(Intercept)"]
    for c in numeric:
        cols.append(df[c].to_numpy(float))
        names.append(c)
    for f in factors:
        lv = pd.Categorical(df[f].astype(str))
        for i, l in enumerate(lv.categories[1:], start=1):
            cols.append((lv.codes == i).astype(float))
            names.append(f"{f}[{l}]")
    Z = np.column_stack(cols)
    # drop exactly collinear columns (e.g. empty factor levels)
    q, r = np.linalg.qr(Z)
    keep = np.abs(np.diag(r)) > 1e-8 * np.abs(np.diag(r)).max()
    return Z[:, keep], [n for n, k in zip(names, keep) if k]


class Residualiser:
    """Thin-QR of Z; residualises y and evaluates single sparse regressors / small dense blocks exactly."""

    def __init__(self, Z):
        self.Q, _ = np.linalg.qr(Z)
        self.N, self.k = Z.shape

    def resid(self, v):
        v = np.asarray(v, float)
        return v - self.Q @ (self.Q.T @ v)

    def test_carriers(self, idx, yr, yryr):
        """OLS t-test of a 0/1 regressor that is 1 at rows idx, given residualised y (yr) and yr'yr."""
        nc = len(idx)
        if nc == 0 or nc == self.N:
            return np.nan, np.nan, np.nan, np.nan
        qx = self.Q[idx].sum(axis=0)
        xrxr = nc - qx @ qx
        if xrxr <= 1e-10:
            return np.nan, np.nan, np.nan, np.nan
        xy = yr[idx].sum()
        b = xy / xrxr
        df = self.N - self.k - 1
        rss = yryr - b * xy
        se = np.sqrt(rss / df / xrxr)
        t = b / se
        return b, se, t, 2 * stats.t.sf(abs(t), df)


def ols_block(Q, y, F):
    """Joint OLS of y on [Z, F] (Z represented by its QR factor Q); returns beta, se, t for the F columns."""
    yr = y - Q @ (Q.T @ y)
    Fr = F - Q @ (Q.T @ F)
    A = Fr.T @ Fr
    b = np.linalg.solve(A, Fr.T @ yr)
    df = len(y) - Q.shape[1] - F.shape[1]
    s2 = (yr @ yr - b @ (Fr.T @ yr)) / df
    se = np.sqrt(np.diag(np.linalg.inv(A)) * s2)
    return b, se, b / se


def stouffer_n(z, ncase_or_n, ncont=None):
    """Omar's sample-size weighted Stouffer (weight sqrt(N_eff)); for quantitative traits pass N only."""
    z = np.asarray(z, float)
    if ncont is None:
        w = np.sqrt(np.asarray(ncase_or_n, float))
    else:
        w = np.sqrt(4.0 / (1.0 / np.asarray(ncase_or_n, float) + 1.0 / np.asarray(ncont, float)))
    zz = np.nansum(z * w) / np.sqrt(np.sum(w ** 2))
    return zz, 2 * stats.norm.sf(abs(zz))
