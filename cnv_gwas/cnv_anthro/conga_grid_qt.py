"""Genome-wide (person-level) univariate grid search and joint model for the CNV gene-dosage axes of a quantitative trait.

Follows the locus-level two-axis report v3 (growth grid, Table G1; height's duplication-axis hit, G3), applied to person-level
burden in one cohort. Shares genes, CNV->gene rules, loci, LOO runs and the burden-matched genome-shift null with
conga_axes_qt.py (ConGA consensus2), which stays the reference implementation of the ConGA definitions.

Measures (per gene weights, summed over the genes a person's DEL or DUP hits):
  mirror grid  g:<genes>:<form>:<sign>:<tissue>
     genes  all | HC (s_het > 0.029)
     form   sum (sum of z) | c1..c4 (count with |z| > t)
     sign   sig (z, or n(z > t) - n(z < -t)) | uns (|z|, or n(|z| > t))
     tissue scaled (sum z / sqrt(n panels), ConGA) | mean (mean over panels)
  constraint   none_H (s_het > 0.029, ConGA) | n95 (top-5% s_het) | pres:n95 (person carries >= 1 top-5% gene)
               | pres:n95HI (>= 1 gene that is top-5% s_het or ClinGen haploinsufficient, score 2-3)
  n_all        total genes hit
A design is a JSON list of variants {"name", "DEL": [terms], "DUP": [terms]}; one OLS per variant per run:
  trait ~ covariates + DEL terms + DUP terms. Readouts: every term on each arm, and mirror = (t_DUP - t_DEL)/sqrt2 for
  any term on both arms. Null, LOO and the no-loci run as in conga_axes_qt.py (same shifts, matched across runs and variants).
The full-run null t of every readout is saved (null_full.npz) so a max-over-grid (selection-adjusted) P can be computed
across task outputs that used the same --seed / --nrep.
"""
import argparse
import json
import os
import time
import numpy as np
import pandas as pd
from scipy import sparse
from common import read_plink_cnv, design_matrix
from conga_axes_qt import load_genes, incidence, assign_loci, SHET_CUT

TCUT = {"c1": 1.0, "c2": 2.0, "c3": 3.0, "c4": 4.0}


def twas_z_both(twas_file, panels, U, gene_cov):
    tw = pd.read_csv(twas_file, sep=r"\s+")
    tw["TWAS.Z"] = pd.to_numeric(tw["TWAS.Z"], errors="coerce"); tw = tw[np.isfinite(tw["TWAS.Z"])]
    st_ = tw["STUDY"].astype(str); brain = st_.str.startswith("GTEX_Brain") | (st_ == "CMC")
    tw = {"gtex_nonbrain": tw[st_.str.startswith("GTEX_") & ~brain], "brain": tw[brain], "all": tw}[panels]
    g = tw.groupby("GENE")["TWAS.Z"]
    sym = pd.concat([gene_cov[["symbol", "entrez"]], U[["symbol", "entrez"]]]).dropna().drop_duplicates("symbol")
    out = {}
    for nm, s in [("scaled", g.agg(lambda v: v.sum() / np.sqrt(len(v)))), ("mean", g.mean())]:
        zz = pd.DataFrame({"symbol": s.index, "z": s.values}).merge(sym, on="symbol").drop_duplicates("entrez").set_index("entrez").z
        out[nm] = U.entrez.map(zz).to_numpy(float)
    return out, int(tw.STUDY.nunique())


def measure(name, Z, shet, q95, hi):
    """Per-gene weight vector for a term name (see module docstring)."""
    hc = np.nan_to_num(shet, nan=-1) > SHET_CUT
    if name == "n_all":
        return np.ones(len(shet))
    if name == "none_H":
        return hc.astype(float)
    if name in ("n95", "pres:n95"):
        return (np.nan_to_num(shet, nan=-1) >= q95).astype(float)
    if name == "pres:n95HI":
        return ((np.nan_to_num(shet, nan=-1) >= q95) | hi).astype(float)
    p = name.split(":")
    assert p[0] == "g" and len(p) == 5, name
    _, genes, form, sign, tissue = p
    z = np.nan_to_num(Z[tissue], nan=0.0)
    if form == "sum":
        w = z if sign == "sig" else np.abs(z)
    else:
        t = TCUT[form]
        w = ((z > t).astype(float) - (z < -t)) if sign == "sig" else (np.abs(z) > t).astype(float)
    return w * hc if genes == "HC" else w


class RunData:
    """Sufficient statistics for one CNV set. Linear terms use Q'G (gene space); presence terms use the carriers' Q rows."""

    def __init__(self, Gd, Gu, Q, yr, yy, N, k):
        rows = np.flatnonzero((Gd.getnnz(axis=1) + Gu.getnnz(axis=1)) > 0)
        self.G = {"DEL": Gd[rows].tocsr(), "DUP": Gu[rows].tocsr()}
        self.QtG = {"DEL": np.asarray((Gd.T @ Q).T), "DUP": np.asarray((Gu.T @ Q).T)}
        self.Qr = Q[rows]
        self.yr = yr[rows]; self.yy = yy; self.N = N; self.k = k

    def cols(self, dosage, terms, Wfun):
        F, QF = [], []
        for t in terms:
            w = Wfun(t)
            f = np.asarray(self.G[dosage] @ w).ravel()
            if t.startswith("pres:"):
                f = (f > 0).astype(float); QF.append(self.Qr.T @ f)
            else:
                QF.append(self.QtG[dosage] @ w)
            F.append(f)
        return F, QF

    def fit(self, var, Wfun):
        F, QF = [], []
        for d in ("DEL", "DUP"):
            a, b = self.cols(d, var[d], Wfun); F += a; QF += b
        F = np.column_stack(F); QF = np.column_stack(QF)
        A = F.T @ F - QF.T @ QF; Fy = F.T @ self.yr
        p = len(Fy); t = np.full(p, np.nan); bb = t.copy(); se = t.copy()
        keep = np.diag(A) > 1e-8
        if keep.sum() == 0:
            return bb, se, t
        try:
            Ai = np.linalg.inv(A[np.ix_(keep, keep)])
        except np.linalg.LinAlgError:
            return bb, se, t
        b = Ai @ Fy[keep]; s2 = (self.yy - b @ Fy[keep]) / (self.N - self.k - keep.sum())
        bb[keep] = b; se[keep] = np.sqrt(np.diag(Ai) * s2); t[keep] = b / se[keep]
        return bb, se, t


def readout_names(var):
    n = [f"{x}|DEL" for x in var["DEL"]] + [f"{x}|DUP" for x in var["DUP"]]
    return n + [f"{x}|mirror" for x in var["DEL"] if x in var["DUP"] and x != "n_all"]


def readouts(t, var):
    nd = len(var["DEL"]); r = list(t)
    for j, x in enumerate(var["DEL"]):
        if x in var["DUP"] and x != "n_all":
            r.append((t[nd + var["DUP"].index(x)] - t[j]) / np.sqrt(2))
    return np.array(r)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pheno", required=True); ap.add_argument("--trait", required=True); ap.add_argument("--cnv", required=True)
    ap.add_argument("--twas", required=True); ap.add_argument("--panels", default="gtex_nonbrain")
    ap.add_argument("--exons", required=True); ap.add_argument("--gene_cov", required=True); ap.add_argument("--clingen", required=True)
    ap.add_argument("--loci", required=True); ap.add_argument("--loo_runs", required=True)
    ap.add_argument("--design", required=True, help="JSON list of variants"); ap.add_argument("--only", default="", help="comma list of variant names")
    ap.add_argument("--nrep", type=int, default=1000); ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--numeric", default="age,age2,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10")
    ap.add_argument("--factors", default="sex,batch"); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True); t0 = time.time()
    design = json.load(open(a.design))
    if a.only:
        keep = set(a.only.split(",")); design = [v for v in design if v["name"] in keep]
    ph = pd.read_csv(a.pheno, sep="\t", dtype={"IID": str})
    num = [c for c in a.numeric.split(",") if c]; fac = [c for c in a.factors.split(",") if c]
    ph = ph[np.isfinite(ph[a.trait])].dropna(subset=num + fac).reset_index(drop=True)
    Zc, _ = design_matrix(ph, num, fac); Q, _ = np.linalg.qr(Zc)
    y = ph[a.trait].to_numpy(float); yr = y - Q @ (Q.T @ y); yy = float(yr @ yr); N, k = Zc.shape
    cnv = read_plink_cnv(a.cnv, keep_iids=ph.IID); cnv = cnv[cnv.CHR <= 22].reset_index(drop=True)
    cnv["row"] = pd.Series(np.arange(N), index=ph.IID).loc[cnv.IID].to_numpy()
    U, ex, tx = load_genes(a.exons); M = len(U)
    gc = pd.read_csv(a.gene_cov, sep="\t", dtype={"entrez": str}); gc.loc[gc.shet_imputed == 1, "shet"] = np.nan
    shet = U[["entrez"]].merge(gc[["entrez", "shet"]], on="entrez", how="left").shet.to_numpy()
    q95 = float(np.nanquantile(gc.shet.to_numpy(float), 0.95))
    cg = pd.read_csv(a.clingen, sep="\t", skiprows=5, dtype=str, quoting=3); cg.columns = [c.strip().lstrip("#") for c in cg.columns]
    hi_sym = set(cg.loc[cg["Haploinsufficiency Score"].isin(["2", "3"]), cg.columns[0]])
    hi = U.symbol.isin(hi_sym).to_numpy()
    Z, npan = twas_z_both(a.twas, a.panels, U, gc)
    info = dict(trait=a.trait, N=int(N), panels=a.panels, n_panels=npan, q95_shet=q95, clingen_HI_genes_in_circle=int(hi.sum()),
                genes_twas=int(np.isfinite(Z["scaled"]).sum()), variants=[v["name"] for v in design], nrep=a.nrep, seed=a.seed)

    # runs: full, 18 LOO, no_loci (as conga_axes_qt.py)
    loci = pd.read_csv(a.loci, sep="\t"); lab = assign_loci(cnv, loci)
    lr = pd.read_csv(a.loo_runs, sep="\t", keep_default_na=False)
    runs = {"full": np.ones(len(cnv), bool)}; is_del = cnv.TYPE.to_numpy() == 1; union = np.zeros(len(cnv), bool)
    for r in lr.itertuples():
        md = set(x for x in str(r.del_members).split(",") if x); mu = set(x for x in str(r.dup_members).split(",") if x)
        drop = (is_del & np.isin(lab, list(md))) | (~is_del & np.isin(lab, list(mu))); runs[f"loo_{r.run}"] = ~drop; union |= drop
    runs["no_loci"] = ~union

    # null shifts: identical to conga_axes_qt.py (scaled z categories), so every task with the same seed/nrep shares them
    Gd0 = incidence(cnv[cnv.TYPE == 1], ex, "overlap", N, M); Gu0 = incidence(cnv[cnv.TYPE == 3], tx, "within", N, M)
    cc = np.asarray(Gd0.sum(0)).ravel() + np.asarray(Gu0.sum(0)).ravel()
    zs = np.nan_to_num(Z["scaled"], nan=0.0)
    cats = [zs > 3, zs < -3, np.nan_to_num(shet, nan=-1) > SHET_CUT]
    obs_b = np.array([cc[m].mean() if m.any() else 0 for m in cats])
    perms = []
    for rep in range(1, a.nrep + 1):
        rng = np.random.default_rng(a.seed * 100003 + rep); best = None
        for _ in range(200):
            kk = int(rng.integers(1000, M - 1000)); perm = (np.arange(M) - kk) % M
            nb = np.array([cc[m[perm]].mean() if m.any() else 0 for m in cats]); use = obs_b > 0
            dev = abs(np.median(nb[use] / obs_b[use]) - 1)
            if best is None or dev < best[1]:
                best = (perm, dev)
            if dev <= 0.1:
                break
        perms.append(best[0])

    terms = sorted({x for v in design for d in ("DEL", "DUP") for x in v[d]})
    W0 = {x: measure(x, Z, shet, q95, hi) for x in terms}
    obs, nul = [], {}
    for rname, sel in runs.items():
        c_ = cnv[sel]
        Gd = Gd0 if rname == "full" else incidence(c_[c_.TYPE == 1], ex, "overlap", N, M)
        Gu = Gu0 if rname == "full" else incidence(c_[c_.TYPE == 3], tx, "within", N, M)
        RD = RunData(Gd, Gu, Q, yr, yy, N, k)
        for v in design:
            names = readout_names(v)
            b, se, t = RD.fit(v, lambda x: W0[x])
            ro = readouts(t, v)
            for j, nm in enumerate(names):
                obs.append(dict(run=rname, variant=v["name"], readout=nm, beta=b[j] if j < len(b) else np.nan,
                                se=se[j] if j < len(se) else np.nan, t=ro[j]))
            NU = np.empty((a.nrep, len(names)))
            for i, perm in enumerate(perms):
                _, _, tp = RD.fit(v, lambda x: W0[x][perm]); NU[i] = readouts(tp, v)
            nul[(rname, v["name"])] = NU
        print(f"[{time.time() - t0:7.1f}s] run {rname} done", flush=True)
    obs = pd.DataFrame(obs); obs.to_csv(os.path.join(a.out, "observed_terms.tsv"), sep="\t", index=False, float_format="%.6g")

    loo = [r for r in runs if r.startswith("loo_")]; out = []; per = []
    ot = obs.set_index(["run", "variant", "readout"]).t
    for v in design:
        for j, nm in enumerate(readout_names(v)):
            def emp(rn):
                o = ot[(rn, v["name"], nm)]; nu = nul[(rn, v["name"])][:, j]; nu = nu[np.isfinite(nu)]
                if not np.isfinite(o) or len(nu) < 10:
                    return o, np.nan, np.nan
                mu, sd = nu.mean(), nu.std(ddof=1)
                return o, (o - mu) / sd, (1 + np.sum(np.abs(nu - mu) >= abs(o - mu))) / (len(nu) + 1)
            o, ze, pe = emp("full"); s = np.sign(ze)
            lz = {rn: emp(rn)[1] for rn in loo}
            for rn in runs:
                per.append(dict(variant=v["name"], readout=nm, run=rn, t=ot[(rn, v["name"], nm)], z_emp=emp(rn)[1]))
            held = sum(1 for q in lz.values() if np.isfinite(q) and abs(q) > 1.96 and np.sign(q) == s)
            rev = sum(1 for q in lz.values() if np.isfinite(q) and np.sign(q) == -s)
            fin = {rn: q for rn, q in lz.items() if np.isfinite(q)}
            wk = min(fin, key=lambda rn: s * fin[rn]) if fin else ""
            T_obs = np.array([ot[(rn, v["name"], nm)] for rn in loo]); T_nul = np.column_stack([nul[(rn, v["name"])][:, j] for rn in loo])
            mo = np.nanmin(s * T_obs); mn = np.nanmin(s * T_nul, axis=1); pmin = (1 + np.sum(mn >= mo)) / (len(mn) + 1)
            o2, z2, p2 = emp("no_loci")
            out.append(dict(trait=a.trait, panels=a.panels, variant=v["name"], readout=nm, t_analytic=o, z_emp=ze, P_emp=pe,
                            loo_held=held, loo_n=len(loo), loo_reversed=rev, weakest=wk.replace("loo_", ""), weakest_z=fin.get(wk, np.nan),
                            P_minLOO_matched=pmin, noloci_z_emp=z2, noloci_P_emp=p2, nrep=a.nrep))
    pd.DataFrame(out).to_csv(os.path.join(a.out, "grid_summary.tsv"), sep="\t", index=False, float_format="%.6g")
    pd.DataFrame(per).to_csv(os.path.join(a.out, "per_run_z.tsv"), sep="\t", index=False, float_format="%.6g")
    np.savez_compressed(os.path.join(a.out, "null_full.npz"), **{f"{v['name']}": nul[("full", v["name"])] for v in design},
                        **{f"names__{v['name']}": np.array(readout_names(v)) for v in design})
    info["minutes"] = (time.time() - t0) / 60
    json.dump(info, open(os.path.join(a.out, "run_info.json"), "w"), indent=1)
    print(pd.DataFrame(out)[["variant", "readout", "z_emp", "P_emp", "loo_held", "noloci_z_emp"]].to_string())


if __name__ == "__main__":
    main()
