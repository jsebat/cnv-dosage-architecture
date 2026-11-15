"""Reference (naive) implementations kept runnable, and the equivalence + timing check of the fast code.

Reference = one full dense OLS per test (numpy lstsq on [covariates, regressor(s)], N x (k+p)), i.e. what a per-probe
lm()/glm() loop does. Fast = cnv_gwas_qt.py (Frisch-Waugh-Lovell with the covariate QR reused, sparse carriers,
incremental carrier sums along the chromosome) and conga_axes_qt.RunData (gene-space sufficient statistics).
Both compute the same OLS Wald t-test; no statistic is changed.

Usage: python reference_check.py --pheno ... --cnv ... --bp_dir ... --gwas_out <fast chr22 output dir> --out <dir>
        [--n_probes 300] [axes options as in conga_axes_qt.py]
"""
import argparse
import json
import os
import time
import numpy as np
import pandas as pd
from scipy import stats
from common import read_plink_cnv, design_matrix


def naive_ols(Z, X, y):
    """Full OLS of y on [Z, X]; returns beta, se, t, p for the X columns."""
    A = np.column_stack([Z, X])
    b, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A @ b
    df = len(y) - A.shape[1]
    s2 = r @ r / df
    cov = s2 * np.linalg.inv(A.T @ A)
    k = Z.shape[1]
    be = b[k:]; se = np.sqrt(np.diag(cov)[k:]); t = be / se
    return be, se, t, 2 * stats.t.sf(np.abs(t), df)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pheno", required=True); ap.add_argument("--trait", default="height_rint")
    ap.add_argument("--cnv", required=True); ap.add_argument("--gwas_out", required=True)
    ap.add_argument("--chrom", type=int, default=22); ap.add_argument("--n_probes", type=int, default=300)
    ap.add_argument("--numeric", default="age,age2,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10")
    ap.add_argument("--factors", default="sex,batch")
    ap.add_argument("--axes", action="store_true"); ap.add_argument("--twas"); ap.add_argument("--exons"); ap.add_argument("--gene_cov")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    ph = pd.read_csv(a.pheno, sep="\t", dtype={"IID": str})
    num = a.numeric.split(","); fac = [f for f in a.factors.split(",") if f]
    ph = ph[np.isfinite(ph[a.trait])].dropna(subset=num + fac).reset_index(drop=True)
    Z, _ = design_matrix(ph, num, fac); y = ph[a.trait].to_numpy(float)
    rows = pd.Series(np.arange(len(ph)), index=ph.IID)
    cnv = read_plink_cnv(a.cnv, keep_iids=ph.IID); cnv = cnv[cnv.CHR == a.chrom]
    rep = {"N": len(ph), "k": Z.shape[1]}
    # ---------------- probe GWAS
    out = []
    for tname, tc in [("del", 1), ("dup", 3)]:
        F = pd.read_csv(os.path.join(a.gwas_out, f"sumstats_{tname}.tsv.gz"), sep="\t")
        F = F[(F.CHR == a.chrom) & np.isfinite(F.BETA)]
        F = F.sample(min(a.n_probes, len(F)), random_state=1).sort_values("BP")
        c = cnv[cnv.TYPE == tc]
        t0 = time.time()
        for r in F.itertuples():
            car = rows.loc[c[(c.BP1 <= r.BP) & (c.BP2 >= r.BP)].IID.unique()].to_numpy()
            x = np.zeros(len(ph)); x[car] = 1
            b, se, t, p = naive_ols(Z, x[:, None], y)
            out.append((tname, r.BP, r.NCNVs, len(car), r.BETA, b[0], r.SE, se[0], r.P, p[0]))
        rep[f"naive_sec_per_probe_{tname}"] = (time.time() - t0) / len(F)
    D = pd.DataFrame(out, columns=["type", "BP", "ncnv_fast", "ncnv_ref", "beta_fast", "beta_ref", "se_fast", "se_ref", "p_fast", "p_ref"])
    D.to_csv(os.path.join(a.out, "gwas_fast_vs_reference.tsv"), sep="\t", index=False)
    rep["gwas_probes_compared"] = len(D)
    rep["gwas_ncnv_mismatch"] = int((D.ncnv_fast != D.ncnv_ref).sum())
    rep["gwas_max_rel_diff_beta"] = float(np.max(np.abs(D.beta_fast - D.beta_ref) / np.abs(D.beta_ref)))
    rep["gwas_max_rel_diff_se"] = float(np.max(np.abs(D.se_fast - D.se_ref) / D.se_ref))
    rep["gwas_max_abs_diff_log10p"] = float(np.max(np.abs(np.log10(D.p_fast) - np.log10(D.p_ref))))
    # ---------------- axes model (one run, observed annotation + 3 shifted annotations)
    if a.axes:
        import conga_axes_qt as cx
        U, ex, tx = cx.load_genes(a.exons); M = len(U)
        gc = pd.read_csv(a.gene_cov, sep="\t", dtype={"entrez": str}); gc.loc[gc.shet_imputed == 1, "shet"] = np.nan
        shet = U[["entrez"]].merge(gc[["entrez", "shet"]], on="entrez", how="left").shet.to_numpy()
        z, _, _ = cx.twas_gene_z(a.twas, "gtex_nonbrain", U, gc); W0 = cx.weights(z, shet)
        call = read_plink_cnv(a.cnv, keep_iids=ph.IID); call = call[call.CHR <= 22].copy(); call["row"] = rows.loc[call.IID].to_numpy()
        Gd = cx.incidence(call[call.TYPE == 1], ex, "overlap", len(ph), M); Gu = cx.incidence(call[call.TYPE == 3], tx, "within", len(ph), M)
        Q, _ = np.linalg.qr(Z); yr = y - Q @ (Q.T @ y)
        t0 = time.time(); RD = cx.RunData(Gd, Gu, Q, yr, float(yr @ yr), len(ph), Z.shape[1]); rep["axes_rundata_setup_sec"] = time.time() - t0
        cols = cx.VARIANTS["primary"]; res = []
        rng = np.random.default_rng(7)
        for rep_i, k in enumerate([0] + list(rng.integers(1000, M - 1000, 3))):
            perm = (np.arange(M) - k) % M
            Wm = np.column_stack([W0[c][perm] for c in cols])
            t0 = time.time(); _, _, tf = RD.fit({"DEL": Wm, "DUP": Wm}); tfast = time.time() - t0
            t0 = time.time()
            Fm = np.column_stack([np.asarray(Gd @ Wm), np.asarray(Gu @ Wm)])
            _, _, tr, _ = naive_ols(Z, Fm, y); tref = time.time() - t0
            res.append(dict(shift=int(k), max_abs_diff_t=float(np.nanmax(np.abs(tf - tr))), sec_fast=tfast, sec_ref=tref))
        rep["axes"] = res
    json.dump(rep, open(os.path.join(a.out, "reference_check.json"), "w"), indent=1, default=float)
    print(json.dumps(rep, indent=1, default=float))


if __name__ == "__main__":
    main()
