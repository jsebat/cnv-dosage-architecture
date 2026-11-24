"""Person-level burden test of the two ConGA gene-dosage axes (Engchuan et al., companion manuscript) on a
quantitative trait, with the genotype-preserving genome-shift null, leave-one-locus-out (18 runs) and an
all-18-loci-removed run.

Axis definitions are ConGA's consensus2 design (conga_joint.R), unchanged:
  mirror measure     count_A = n(genes hit with TWAS z > 3) - n(genes hit with z < -3)
  constraint measure none_H  = n(genes hit with GeneBayes s_het > 0.029)   (imputed s_het excluded)
  total genes hit    n_all
CNV -> gene (ConGA v2 defaults): DEL hits a gene if it overlaps any refGene exon by >= 1 bp; DUP hits a gene if it
spans an entire refGene transcript (any isoform). Autosomes only.
Model (one fit per run): trait ~ covariates + [n_all, count_A, none_H]_DEL + [n_all, count_A, none_H]_DUP (OLS).
Readouts: count_A DEL arm, DUP arm, mirror net = (z_DUP - z_DEL)/sqrt(2); none_H DEL and DUP arms.
Secondary (height-specific, quantified alongside, never selected over the primary):
  'unsigned' = count_A replaced by n(|z| > 3);  'split' = count_A replaced by n(z > 3) and n(z < -3) jointly.
Null: ConGA true_perm -- one circular shift of the autosomal gene circle per replicate (k ~ U(1000, N-1000)), redrawn
up to 200 times until the median null/observed carrier burden of the annotated gene sets is within 10%; the CNV->gene
incidence is fixed and the annotations (TWAS z and s_het together) move. The same shift is applied to every run
(matched null for LOO and the no-loci run).
Gene-space algebra: with Q = QR factor of the covariates, every replicate needs only G'G, Q'G and G'y_resid, so no
N-dimensional work is repeated.
"""
import argparse
import json
import os
import time
import numpy as np
import pandas as pd
from scipy import sparse, stats
from common import read_plink_cnv, design_matrix, rcg

SHET_CUT = 0.0290
Z_CUT = 3.0


# ---------------------------------------------------------------- genes
def load_genes(exon_file):
    ex = pd.read_csv(exon_file, sep="\t", header=None, names=["chr", "start", "end", "tx", "symbol", "entrez"], dtype={"entrez": str})
    ex = ex[ex.chr.isin([f"chr{i}" for i in range(1, 23)])].copy()
    ex["c"] = ex.chr.str.replace("chr", "").astype(int)
    # gene circle: each gene on the chromosome where it has most exon rows, ordered chr1..22 by first start
    cnt = ex.groupby(["entrez", "c"]).size().rename("n").reset_index().sort_values(["entrez", "n"], ascending=[True, False])
    home = cnt.drop_duplicates("entrez")[["entrez", "c"]]
    pos = ex.merge(home, on=["entrez", "c"]).groupby(["entrez", "c"]).start.min().rename("pos").reset_index()
    U = pos.sort_values(["c", "pos"]).reset_index(drop=True)
    U["gidx"] = np.arange(len(U))
    gi = dict(zip(U.entrez, U.gidx))
    ex["gidx"] = ex.entrez.map(gi)
    tx = ex.groupby(["tx", "c", "gidx"]).agg(ts=("start", "min"), te=("end", "max")).reset_index()
    U = U.merge(ex.drop_duplicates("entrez")[["entrez", "symbol"]], on="entrez", how="left")
    return U, ex[["c", "start", "end", "gidx"]], tx[["c", "ts", "te", "gidx"]]


def incidence(calls, feats, mode, n_rows, n_genes):
    """Sparse person x gene 0/1. mode 'overlap': feature [start,end] overlaps call by >= 1 bp (DEL/exon);
    'within': feature [ts,te] lies inside the call (DUP/transcript)."""
    R, C = [], []
    for c, cc in calls.groupby("CHR"):
        f = feats[feats.c == c]
        if len(f) == 0:
            continue
        s = f.iloc[:, 1].to_numpy(); e = f.iloc[:, 2].to_numpy(); g = f.gidx.to_numpy()
        o = np.argsort(s); s, e, g = s[o], e[o], g[o]
        maxlen = (e - s).max()
        for r, b1, b2 in zip(cc.row.to_numpy(), cc.BP1.to_numpy(), cc.BP2.to_numpy()):
            if mode == "overlap":
                lo = np.searchsorted(s, b1 - maxlen, "left"); hi = np.searchsorted(s, b2, "right")
                m = e[lo:hi] >= b1
            else:
                lo = np.searchsorted(s, b1, "left"); hi = np.searchsorted(s, b2, "right")
                m = e[lo:hi] <= b2
            gg = np.unique(g[lo:hi][m])
            R.append(np.full(len(gg), r)); C.append(gg)
    if not R:
        return sparse.csr_matrix((n_rows, n_genes))
    R = np.concatenate(R); C = np.concatenate(C)
    M = sparse.csr_matrix((np.ones(len(R)), (R, C)), shape=(n_rows, n_genes))
    M.data[:] = 1.0                       # binary: a gene hit twice counts once
    return M


# ---------------------------------------------------------------- loci for LOO (port of assign_loci.awk)
def assign_loci(cnv, loci):
    """Per-call locus label for the locus-removal runs. Carriers come from Molly's recurrent CNV genotyping
    (recurrentCNVgenotyping/aou, rcg.label_calls); this module has no overlap rule of its own (Jon 2026-11-24)."""
    return rcg.label_calls(cnv, loci)


# ---------------------------------------------------------------- TWAS
def twas_gene_z(twas_file, panels, U, gene_cov):
    tw = pd.read_csv(twas_file, sep=r"\s+")
    tw = tw[np.isfinite(pd.to_numeric(tw["TWAS.Z"], errors="coerce"))]
    tw["TWAS.Z"] = tw["TWAS.Z"].astype(float)
    st_ = tw["STUDY"].astype(str)
    brain = st_.str.startswith("GTEX_Brain") | (st_ == "CMC")
    if panels == "gtex_nonbrain":
        tw = tw[st_.str.startswith("GTEX_") & ~brain]
    elif panels == "brain":
        tw = tw[brain]
    elif panels != "all":
        raise SystemExit("--panels gtex_nonbrain | brain | all")
    z = tw.groupby("GENE")["TWAS.Z"].agg(lambda v: v.sum() / np.sqrt(len(v)))          # 'scaled' tissue combination
    sym = pd.concat([gene_cov[["symbol", "entrez"]], U[["symbol", "entrez"]]]).dropna().drop_duplicates("symbol")   # gene_covariates first, then refGene
    zz = pd.DataFrame({"symbol": z.index, "z": z.values}).merge(sym, on="symbol").drop_duplicates("entrez")
    out = np.full(len(U), np.nan)
    m = U.entrez.isin(zz.entrez)
    out[m.to_numpy()] = zz.set_index("entrez").loc[U.entrez[m], "z"].to_numpy()
    return out, int(tw.STUDY.nunique()), sorted(tw.STUDY.unique().tolist())


# ---------------------------------------------------------------- model in gene space
VARIANTS = {"primary": ["n_all", "count_A", "none_H"],
            "unsigned": ["n_all", "count_U", "none_H"],
            "split": ["n_all", "count_P", "count_N", "none_H"]}


def weights(z, shet):
    zf = np.nan_to_num(z, nan=0.0)
    return {"n_all": np.ones(len(z)), "count_A": (zf > Z_CUT).astype(float) - (zf < -Z_CUT).astype(float),
            "count_U": (np.abs(zf) > Z_CUT).astype(float), "count_P": (zf > Z_CUT).astype(float), "count_N": (zf < -Z_CUT).astype(float),
            "none_H": (np.nan_to_num(shet, nan=-1) > SHET_CUT).astype(float)}


class RunData:
    """Sufficient statistics for one CNV set (one run). Only people with >= 1 gene hit enter F; Q'F comes from Q'G."""

    def __init__(self, Gd, Gu, Q, yr, yy, N, k):
        rows = np.flatnonzero((Gd.getnnz(axis=1) + Gu.getnnz(axis=1)) > 0)
        self.G = {"DEL": Gd[rows].tocsr(), "DUP": Gu[rows].tocsr()}
        self.QtG = {"DEL": np.asarray((Gd.T @ Q).T), "DUP": np.asarray((Gu.T @ Q).T)}     # k x M
        self.yr = yr[rows]
        self.yy = yy; self.N = N; self.k = k
        self.carriers = {"DEL": int((Gd.getnnz(axis=1) > 0).sum()), "DUP": int((Gu.getnnz(axis=1) > 0).sum()),
                         "any_gene_hit": int(len(rows))}

    def fit(self, W):
        """W: dict dosage -> M x p weights. Returns beta, se, t (DEL columns then DUP columns)."""
        F = np.hstack([np.asarray(self.G[d] @ W[d]) for d in ("DEL", "DUP")])
        QtF = np.hstack([self.QtG[d] @ W[d] for d in ("DEL", "DUP")])
        A = F.T @ F - QtF.T @ QtF
        Fy = F.T @ self.yr
        keep = np.diag(A) > 1e-8
        t = np.full(len(Fy), np.nan); bb = np.full(len(Fy), np.nan); se = np.full(len(Fy), np.nan)
        if keep.sum() == 0:
            return bb, se, t
        try:
            Ai = np.linalg.inv(A[np.ix_(keep, keep)])
        except np.linalg.LinAlgError:
            return bb, se, t
        b = Ai @ Fy[keep]
        s2 = (self.yy - b @ Fy[keep]) / (self.N - self.k - keep.sum())
        bb[keep] = b; se[keep] = np.sqrt(np.diag(Ai) * s2); t[keep] = b / se[keep]
        return bb, se, t


def readouts(t, cols):
    """t over DEL cols then DUP cols -> named readouts."""
    p = len(cols); r = {}
    for j, c in enumerate(cols):
        r[f"{c}_DEL"] = t[j]; r[f"{c}_DUP"] = t[p + j]
    for c in ["count_A", "count_U", "count_P", "count_N"]:
        if f"{c}_DEL" in r:
            r[f"{c}_mirror"] = (r[f"{c}_DUP"] - r[f"{c}_DEL"]) / np.sqrt(2)
    return r


_RUN_CTX = {}


def _process_run(rname):
    """Observed fit + all null replicates for one run (full / LOO / no_loci)."""
    X = _RUN_CTX
    c_ = X["cnv"][X["runs"][rname]]
    Gd = X["Gd0"] if rname == "full" else incidence(c_[c_.TYPE == 1], X["ex"], "overlap", X["N"], X["M"])
    Gu = X["Gu0"] if rname == "full" else incidence(c_[c_.TYPE == 3], X["tx"], "within", X["N"], X["M"])
    RD = RunData(Gd, Gu, X["Q"], X["yr"], X["yy"], X["N"], X["k"])
    o_rows, nst = [], {}
    for vname, cols in VARIANTS.items():
        W = X["Wmat"](X["W0"], cols)
        b, se, t = RD.fit({"DEL": W, "DUP": W})
        ro = readouts(t, cols)
        names = [f"{c}_DEL" for c in cols] + [f"{c}_DUP" for c in cols]
        for j, nm in enumerate(names):
            o_rows.append(dict(run=rname, variant=vname, term=nm, beta=b[j], se=se[j], t=t[j]))
        for kx in [kx for kx in ro if kx.endswith("mirror")]:
            o_rows.append(dict(run=rname, variant=vname, term=kx, beta=np.nan, se=np.nan, t=ro[kx]))
        NUL = np.empty((X["nrep"], len(ro)))
        for i, perm in enumerate(X["perms"]):
            Wp = X["Wmat"](X["W0"], cols, perm)
            _, _, tp = RD.fit({"DEL": Wp, "DUP": Wp})
            NUL[i] = list(readouts(tp, cols).values())
        nst[(rname, vname)] = (list(ro.keys()), NUL)
    print(f"run {rname} done", flush=True)
    return o_rows, nst, RD.carriers


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=1, help="parallel runs (fork); results identical to --workers 1")
    ap.add_argument("--pheno", required=True); ap.add_argument("--trait", required=True)
    ap.add_argument("--cnv", required=True); ap.add_argument("--twas", required=True)
    ap.add_argument("--panels", default="gtex_nonbrain"); ap.add_argument("--exons", required=True)
    ap.add_argument("--gene_cov", required=True); ap.add_argument("--loci", required=True, help="omar_loci_all.tsv")
    ap.add_argument("--loo_runs", required=True, help="loo2_runs.tsv (18 runs)")
    ap.add_argument("--nrep", type=int, default=2000); ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--numeric", default="age,age2,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10")
    ap.add_argument("--factors", default="sex,batch"); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    ph = pd.read_csv(a.pheno, sep="\t", dtype={"IID": str})
    num = [c for c in a.numeric.split(",") if c]; fac = [c for c in a.factors.split(",") if c]
    ph = ph[np.isfinite(ph[a.trait])].dropna(subset=num + fac).reset_index(drop=True)
    Z, _ = design_matrix(ph, num, fac)
    Q, _ = np.linalg.qr(Z)
    y = ph[a.trait].to_numpy(float); yr = y - Q @ (Q.T @ y); yy = float(yr @ yr)
    N, k = Z.shape
    cnv = read_plink_cnv(a.cnv, keep_iids=ph.IID)
    cnv = cnv[cnv.CHR <= 22].reset_index(drop=True)
    cnv["row"] = pd.Series(np.arange(N), index=ph.IID).loc[cnv.IID].to_numpy()
    U, ex, tx = load_genes(a.exons)
    M = len(U)
    gc = pd.read_csv(a.gene_cov, sep="\t", dtype={"entrez": str})
    gc.loc[gc.shet_imputed == 1, "shet"] = np.nan
    shet = U[["entrez"]].merge(gc[["entrez", "shet"]], on="entrez", how="left").shet.to_numpy()
    z, npan, panlist = twas_gene_z(a.twas, a.panels, U, gc)
    W0 = weights(z, shet)
    info = {"trait": a.trait, "N": int(N), "covariate_columns": int(k), "genes_in_circle": int(M), "panels": a.panels,
            "n_panels": npan, "panel_list": panlist, "genes_with_twas_z": int(np.isfinite(z).sum()),
            "genes_z_gt3": int((z > 3).sum()), "genes_z_lt_m3": int((z < -3).sum()), "genes_shet_gt_cut": int(W0["none_H"].sum()),
            "calls_DEL": int((cnv.TYPE == 1).sum()), "calls_DUP": int((cnv.TYPE == 3).sum())}

    # ---- runs: full, 18 LOO, all loci removed
    loci = pd.read_csv(a.loci, sep="\t")
    lab = assign_loci(cnv, loci)
    lr = pd.read_csv(a.loo_runs, sep="\t", keep_default_na=False)
    runs = {"full": np.ones(len(cnv), bool)}
    is_del = cnv.TYPE.to_numpy() == 1
    union = np.zeros(len(cnv), bool)
    for r in lr.itertuples():
        md = set(x for x in str(r.del_members).split(",") if x); mu = set(x for x in str(r.dup_members).split(",") if x)
        drop = (is_del & np.isin(lab, list(md))) | (~is_del & np.isin(lab, list(mu)))
        runs[f"loo_{r.run}"] = ~drop; union |= drop
    runs["no_loci"] = ~union
    info["calls_dropped_per_run"] = {k_: int((~v).sum()) for k_, v in runs.items()}

    # ---- null shifts (ConGA true_perm, burden-matched on the full run)
    Gd0 = incidence(cnv[(cnv.TYPE == 1)], ex, "overlap", N, M)
    Gu0 = incidence(cnv[(cnv.TYPE == 3)], tx, "within", N, M)
    cc = np.asarray(Gd0.sum(0)).ravel() + np.asarray(Gu0.sum(0)).ravel()     # carriers per gene, DEL + DUP
    cats = [W0["count_P"] > 0, W0["count_N"] > 0, W0["none_H"] > 0]
    obs_b = np.array([cc[m].mean() if m.any() else 0 for m in cats])
    perms = []; devs = []
    for rep in range(1, a.nrep + 1):
        rng = np.random.default_rng(a.seed * 100003 + rep); best = None
        for _ in range(200):
            kk = int(rng.integers(1000, M - 1000)); perm = (np.arange(M) - kk) % M      # gene g takes the annotation of perm[g]
            nb = np.array([cc[m[perm]].mean() if m.any() else 0 for m in cats])
            use = obs_b > 0
            dev = abs(np.median(nb[use] / obs_b[use]) - 1)
            if best is None or dev < best[1]:
                best = (perm, dev)
            if dev <= 0.1:
                break
        perms.append(best[0]); devs.append(best[1])
    info["null_burden_dev_median"] = float(np.median(devs)); info["null_burden_dev_max"] = float(np.max(devs))

    def Wmat(Wd, cols, perm=None):
        return np.column_stack([(Wd[c] if perm is None else Wd[c][perm]) for c in cols])

    obs_rows, null_store = [], {}
    global _RUN_CTX
    _RUN_CTX = dict(cnv=cnv, runs=runs, Gd0=Gd0, Gu0=Gu0, ex=ex, tx=tx, N=N, M=M, Q=Q, yr=yr, yy=yy, k=k, W0=W0, perms=perms,
                    nrep=a.nrep, Wmat=Wmat)
    names_runs = list(runs)
    if a.workers > 1:                                  # fork: workers share the inputs copy-on-write (Linux, incl. AoU VMs)
        import multiprocessing as mp
        with mp.get_context("fork").Pool(a.workers) as pool:
            results = pool.map(_process_run, names_runs)
    else:
        results = [_process_run(r) for r in names_runs]
    for rname, (o_rows, nst, carriers) in zip(names_runs, results):   # deterministic: same order and values either way
        obs_rows.extend(o_rows); null_store.update(nst)
        info.setdefault("carriers_per_run", {})[rname] = carriers
    print(f"[{time.time() - t0:7.1f}s] all {len(names_runs)} runs done (workers={a.workers})", flush=True)

    obs = pd.DataFrame(obs_rows)
    obs.to_csv(os.path.join(a.out, "observed_terms.tsv"), sep="\t", index=False, float_format="%.6g")
    # ---- empirical z / P, LOO hold, min-over-LOO with matched null
    out = []
    loo_names = [r for r in runs if r.startswith("loo_")]
    for vname in VARIANTS:
        keys, _ = null_store[("full", vname)]
        for j, term in enumerate(keys):
            def emp(rn):
                o = obs[(obs.run == rn) & (obs.variant == vname) & (obs.term == term)].t.iloc[0]
                nu = null_store[(rn, vname)][1][:, j]; nu = nu[np.isfinite(nu)]
                mu, sd = nu.mean(), nu.std(ddof=1)
                pe = (1 + np.sum(np.abs(nu - mu) >= abs(o - mu))) / (len(nu) + 1)
                return o, (o - mu) / sd, pe
            o, ze, pe = emp("full")
            lz = {rn: emp(rn)[1] for rn in loo_names}
            s = np.sign(ze)
            held = sum(1 for v in lz.values() if np.isfinite(v) and abs(v) > 1.96 and np.sign(v) == s)
            reversed_ = sum(1 for v in lz.values() if np.isfinite(v) and np.sign(v) == -s)
            weakest = min(lz, key=lambda rn: s * lz[rn] if np.isfinite(lz[rn]) else np.inf)
            # min over LOO of the signed analytic t, against the same min in each null replicate
            T_obs = np.array([obs[(obs.run == rn) & (obs.variant == vname) & (obs.term == term)].t.iloc[0] for rn in loo_names])
            T_nul = np.column_stack([null_store[(rn, vname)][1][:, j] for rn in loo_names])
            mo = np.nanmin(s * T_obs); mn = np.nanmin(s * T_nul, axis=1)
            p_min = (1 + np.sum(mn >= mo)) / (len(mn) + 1)
            o_nl, ze_nl, pe_nl = emp("no_loci")
            out.append(dict(trait=a.trait, panels=a.panels, variant=vname, readout=term, t_analytic=o, z_emp=ze, P_emp=pe,
                            loo_held=held, loo_n=len(loo_names), loo_reversed=reversed_, weakest=weakest.replace("loo_", ""),
                            weakest_z=lz[weakest], P_minLOO_matched=p_min, noloci_t=o_nl, noloci_z_emp=ze_nl, noloci_P_emp=pe_nl,
                            nrep=a.nrep))
    res = pd.DataFrame(out)
    res.to_csv(os.path.join(a.out, "axes_summary.tsv"), sep="\t", index=False, float_format="%.6g")
    loo_tab = []
    for vname in VARIANTS:
        keys, _ = null_store[("full", vname)]
        for j, term in enumerate(keys):
            for rn in runs:
                o = obs[(obs.run == rn) & (obs.variant == vname) & (obs.term == term)].t.iloc[0]
                nu = null_store[(rn, vname)][1][:, j]
                loo_tab.append(dict(variant=vname, readout=term, run=rn, t=o, z_emp=(o - np.nanmean(nu)) / np.nanstd(nu, ddof=1)))
    pd.DataFrame(loo_tab).to_csv(os.path.join(a.out, "per_run_z.tsv"), sep="\t", index=False, float_format="%.6g")
    info["minutes"] = (time.time() - t0) / 60
    json.dump(info, open(os.path.join(a.out, "run_info.json"), "w"), indent=1, default=str)
    print(res[res.variant == "primary"].to_string())


if __name__ == "__main__":
    main()
