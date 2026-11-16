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
LOO (--loo_mode): 'all' = matched null for every single-family removal; 'drivers' (Jon 2026-11-15, to save compute) = an observed-only pass
over all removals picks the major families (|dt| >= --major_dt on any readout); each major gets its own matched null, all minor families are
dropped together in one matched null ('minor_all'), plus the all-loci run ('no_loci').
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
    sh = np.nan_to_num(shet, nan=-1); hc = sh > SHET_CUT
    if name.startswith("pres:"):        # presence of any gene in the count measure (person-level 0/1 is applied in RunData.cols)
        return measure(name[5:], Z, shet, q95, hi)
    if name == "n_all":
        return np.ones(len(shet))
    if name == "none_H":
        return hc.astype(float)
    if name == "n95":
        return (sh >= q95).astype(float)
    if name == "n95HI":                 # union: top-5% s_het OR ClinGen HI 2-3
        return ((sh >= q95) | hi).astype(float)
    if name == "HI":
        return hi.astype(float)
    if name == "sh_mid":                # GeneBayes s_het middle tertile (cuts 0.0033 / 0.0296; constraint grid 2026-11-16)
        return ((sh >= 0.0033) & (sh < 0.0296)).astype(float)
    if name == "sh_high":               # top tertile
        return (sh >= 0.0296).astype(float)
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

    def __init__(self, Gd, Gu, Q, yr, yy, N, k, extra=None):
        """extra: optional {"DEL:L5": person x gene incidence from a call-size class, ...} (subsets of Gd / Gu)."""
        rows = np.flatnonzero((Gd.getnnz(axis=1) + Gu.getnnz(axis=1)) > 0)
        allG = {"DEL": Gd, "DUP": Gu, **(extra or {})}
        self.G = {key: g[rows].tocsr() for key, g in allG.items()}
        self.QtG = {key: np.asarray((g.T @ Q).T) for key, g in allG.items()}
        self.Qr = Q[rows]
        self.yr = yr[rows]; self.yy = yy; self.N = N; self.k = k

    def cols(self, dosage, terms, Wfun):
        F, QF = [], []
        for t in terms:
            key, base = (f"{dosage}:{t.split('~')[0]}", t.split("~")[1]) if "~" in t else (dosage, t)   # 'L5~<measure>' = calls hitting > 5 genes
            w = Wfun(base)
            f = np.asarray(self.G[key] @ w).ravel()
            if base.startswith("pres:"):
                f = (f > 0).astype(float); QF.append(self.Qr.T @ f)
            else:
                QF.append(self.QtG[key] @ w)
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
        self.rss = float(self.yy - b @ Fy[keep]); self.npar = int(self.k + keep.sum())   # for log-likelihood / AIC (observed fits)
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
    ap.add_argument("--loo_mode", default="all", choices=["all", "drivers", "all_minor"], help="all = matched null for every single-locus removal; drivers = only major loci individually + minor loci together")
    ap.add_argument("--major_dt", type=float, default=1.0, help="drivers: a locus family is major if dropping it moves any readout t by >= this")
    ap.add_argument("--max_major", type=int, default=8)
    ap.add_argument("--major_loci", default="", help="drivers: fixed comma list of major LOO families (e.g. from UKBB) instead of the observed selection, so cohorts share the same runs for meta-analysis")
    ap.add_argument("--perms_only", action="store_true", help="compute/save the shared shifts (--perms_file) and exit")
    ap.add_argument("--perms_file", default="", help="shared genome shifts: load the shift offsets from this .npy if it exists (cohorts then share null replicates); else compute and save them here")
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
    if a.perms_file and os.path.exists(a.perms_file):
        offs = np.load(a.perms_file)
        if len(offs) < a.nrep:
            raise SystemExit(f"{a.perms_file} holds {len(offs)} shifts < --nrep {a.nrep}")
        perms = [(np.arange(M) - int(kk)) % M for kk in offs[:a.nrep]]
        info["perms"] = f"loaded {a.nrep} shared shifts from {os.path.basename(a.perms_file)}"
    offs_new = []
    for rep in (range(1, a.nrep + 1) if not perms else []):
        rng = np.random.default_rng(a.seed * 100003 + rep); best = None
        for _ in range(200):
            kk = int(rng.integers(1000, M - 1000)); perm = (np.arange(M) - kk) % M
            nb = np.array([cc[m[perm]].mean() if m.any() else 0 for m in cats]); use = obs_b > 0
            dev = abs(np.median(nb[use] / obs_b[use]) - 1)
            if best is None or dev < best[1]:
                best = (perm, dev, kk)
            if dev <= 0.1:
                break
        perms.append(best[0]); offs_new.append(best[2])
    import hashlib
    _offs = np.asarray(offs[:a.nrep] if not offs_new else offs_new, dtype=np.int64)
    info["perms_sha1"] = hashlib.sha1(_offs.tobytes()).hexdigest()[:16]      # identical across cohorts that share shifts
    if offs_new and a.perms_file:
        np.save(a.perms_file, np.array(offs_new, dtype=np.int64)); info["perms"] = f"computed {len(offs_new)} shifts, saved to {os.path.basename(a.perms_file)}"
    if a.perms_only:
        json.dump(info, open(os.path.join(a.out, "perms_info.json"), "w"), indent=1); print("shifts saved:", a.perms_file, info["perms_sha1"]); return

    terms = sorted({x for v in design for d in ("DEL", "DUP") for x in v[d]})
    W0 = {b: measure(b, Z, shet, q95, hi) for b in {x.split("~")[-1] for x in terms}}
    # call-size classes ('L5~' = calls hitting > 5 genes of their dosage, 'S5~' = <= 5; any integer): per-call gene counts under the same CNV->gene rule
    classes = sorted({x.split("~")[0] for x in terms if "~" in x})
    if classes:
        cc_ = cnv.copy(); cc_["row"] = np.arange(len(cc_))
        ng = np.zeros(len(cnv))
        for typ, feats, mode in [(1, ex, "overlap"), (3, tx, "within")]:
            m = (cnv.TYPE == typ).to_numpy()
            ng[m] = incidence(cc_[m], feats, mode, len(cnv), M).getnnz(axis=1)[m]
        cnv["n_genes_call"] = ng
        info["calls_by_class"] = {c: {"DEL": int(((cnv.TYPE == 1) & ((ng > int(c[1:])) if c[0] == "L" else (ng <= int(c[1:])))).sum()),
                                      "DUP": int(((cnv.TYPE == 3) & ((ng > int(c[1:])) if c[0] == "L" else (ng <= int(c[1:])))).sum())} for c in classes}
    def build_rd(rname, sel):
        c_ = cnv[sel]
        Gd = Gd0 if rname == "full" else incidence(c_[c_.TYPE == 1], ex, "overlap", N, M)
        Gu = Gu0 if rname == "full" else incidence(c_[c_.TYPE == 3], tx, "within", N, M)
        extra = {}
        for c in classes:
            kk = int(c[1:]); cm = (c_.n_genes_call > kk) if c[0] == "L" else (c_.n_genes_call <= kk)
            extra[f"DEL:{c}"] = incidence(c_[(c_.TYPE == 1) & cm], ex, "overlap", N, M)
            extra[f"DUP:{c}"] = incidence(c_[(c_.TYPE == 3) & cm], tx, "within", N, M)
        return RunData(Gd, Gu, Q, yr, yy, N, k, extra)

    def observed(RD, rname):
        rows = []
        for v in design:
            b, se, t = RD.fit(v, lambda x: W0[x]); ro = readouts(t, v)
            rss, npar = getattr(RD, "rss", np.nan), getattr(RD, "npar", np.nan)
            ll = -0.5 * RD.N * (np.log(2 * np.pi * rss / RD.N) + 1)          # Gaussian ML log-likelihood of the OLS fit
            for j, nm in enumerate(readout_names(v)):
                rows.append(dict(run=rname, variant=v["name"], readout=nm, beta=b[j] if j < len(b) else np.nan,
                                 se=se[j] if j < len(se) else np.nan, t=ro[j], rss=rss, npar=npar, loglik=ll, AIC=2 * (npar + 1) - 2 * ll))
        return rows

    obs, nul = [], {}
    loo_all = [r for r in runs if r.startswith("loo_")]
    if a.loo_mode == "drivers":
        # Jon 2026-11-15: observed-only pass over every single-family removal; families that move any readout by |dt| >= major_dt are
        # 'major' and get their own matched-null run; all other ('minor') families are dropped together in one matched-null run.
        for rname in ["full"] + loo_all:
            obs += observed(build_rd(rname, runs[rname]), rname)
        O_ = pd.DataFrame(obs); O_ = O_[~O_.readout.str.startswith("n_all")]
        tf = O_[O_.run == "full"].set_index(["variant", "readout"]).t
        dt = {r: np.nanmax(np.abs(O_[O_.run == r].set_index(["variant", "readout"]).t - tf).to_numpy()) for r in loo_all}
        majors = [r for r in sorted(dt, key=lambda r: -dt[r]) if dt[r] >= a.major_dt][:a.max_major]
        if a.major_loci:
            majors = [f"loo_{x}" for x in a.major_loci.split(",") if f"loo_{x}" in runs]
            info["major_loci_source"] = "fixed (--major_loci)"
        minor_drop = np.zeros(len(cnv), bool)
        for r in loo_all:
            if r not in majors:
                minor_drop |= ~runs[r]
        runs["minor_all"] = ~minor_drop
        null_runs = ["full"] + majors + ["minor_all", "no_loci"]
        info["loo_mode"] = "drivers"; info["major_dt"] = a.major_dt; info["observed_dt"] = {r.replace("loo_", ""): round(float(v), 3) for r, v in dt.items()}
        info["major_loci"] = [r.replace("loo_", "") for r in majors]
        obs = [o for o in obs if o["run"] not in null_runs]      # null runs are refitted (observed + null) below
        loo = majors
    else:
        if a.loo_mode == "all_minor":       # Jon 2026-11-15: every single-unit removal AND the drivers-kept null (all non-driver units dropped)
            majors = [f"loo_{x}" for x in a.major_loci.split(",") if f"loo_{x}" in runs]
            if not majors:
                raise SystemExit("--loo_mode all_minor needs --major_loci")
            minor_drop = np.zeros(len(cnv), bool)
            for r in loo_all:
                if r not in majors:
                    minor_drop |= ~runs[r]
            runs["minor_all"] = ~minor_drop; info["major_loci"] = [r.replace("loo_", "") for r in majors]
        null_runs = list(runs); loo = loo_all
        info["loo_mode"] = a.loo_mode
    for rname in null_runs:
        RD = build_rd(rname, runs[rname])
        obs += observed(RD, rname)
        for v in design:
            names = readout_names(v)
            NU = np.empty((a.nrep, len(names)))
            for i, perm in enumerate(perms):
                _, _, tp = RD.fit(v, lambda x: W0[x][perm]); NU[i] = readouts(tp, v)
            nul[(rname, v["name"])] = NU
        print(f"[{time.time() - t0:7.1f}s] run {rname} done", flush=True)
    obs = pd.DataFrame(obs); obs.to_csv(os.path.join(a.out, "observed_terms.tsv"), sep="	", index=False, float_format="%.6g")

    out = []; per = []
    ot = obs.set_index(["run", "variant", "readout"]).t
    for v in design:
        for j, nm in enumerate(readout_names(v)):
            def emp(rn):
                o = ot[(rn, v["name"], nm)]
                if (rn, v["name"]) not in nul:
                    return o, np.nan, np.nan
                nu = nul[(rn, v["name"])][:, j]; nu = nu[np.isfinite(nu)]
                if not np.isfinite(o) or len(nu) < 10:
                    return o, np.nan, np.nan
                mu, sd = nu.mean(), nu.std(ddof=1)
                return o, (o - mu) / sd, (1 + np.sum(np.abs(nu - mu) >= abs(o - mu))) / (len(nu) + 1)
            o, ze, pe = emp("full"); s = np.sign(ze)
            lz = {rn: emp(rn)[1] for rn in loo}
            for rn in runs:
                if (rn, v["name"], nm) in ot.index:
                    per.append(dict(variant=v["name"], readout=nm, run=rn, t=ot[(rn, v["name"], nm)], z_emp=emp(rn)[1]))
            held = sum(1 for q in lz.values() if np.isfinite(q) and abs(q) > 1.96 and np.sign(q) == s)
            rev = sum(1 for q in lz.values() if np.isfinite(q) and np.sign(q) == -s)
            fin = {rn: q for rn, q in lz.items() if np.isfinite(q)}
            wk = min(fin, key=lambda rn: s * fin[rn]) if fin else ""
            if loo:
                T_obs = np.array([ot[(rn, v["name"], nm)] for rn in loo]); T_nul = np.column_stack([nul[(rn, v["name"])][:, j] for rn in loo])
                mo = np.nanmin(s * T_obs); mn = np.nanmin(s * T_nul, axis=1); pmin = (1 + np.sum(mn >= mo)) / (len(mn) + 1)
            else:
                pmin = np.nan
            o2, z2, p2 = emp("no_loci")
            row = dict(trait=a.trait, panels=a.panels, variant=v["name"], readout=nm, t_analytic=o, z_emp=ze, P_emp=pe,
                       loo_held=held, loo_n=len(loo), loo_reversed=rev, weakest=wk.replace("loo_", ""), weakest_z=fin.get(wk, np.nan),
                       P_minLOO_matched=pmin, noloci_z_emp=z2, noloci_P_emp=p2, nrep=a.nrep, loo_mode=info["loo_mode"])
            if info["loo_mode"] in ("drivers", "all_minor"):
                o3, z3, p3 = emp("minor_all")
                row.update(minor_all_z_emp=z3, minor_all_P_emp=p3, major_loci=";".join(info["major_loci"]))
            out.append(row)
    pd.DataFrame(out).to_csv(os.path.join(a.out, "grid_summary.tsv"), sep="\t", index=False, float_format="%.6g")
    pd.DataFrame(per).to_csv(os.path.join(a.out, "per_run_z.tsv"), sep="\t", index=False, float_format="%.6g")
    np.savez_compressed(os.path.join(a.out, "null_runs.npz"), **{f"{rn}__{vn}": arr for (rn, vn), arr in nul.items()})   # every matched-null run, for meta-analysis
    np.savez_compressed(os.path.join(a.out, "null_full.npz"), **{f"{v['name']}": nul[("full", v["name"])] for v in design},
                        **{f"names__{v['name']}": np.array(readout_names(v)) for v in design})
    info["minutes"] = (time.time() - t0) / 60
    json.dump(info, open(os.path.join(a.out, "run_info.json"), "w"), indent=1)
    print(pd.DataFrame(out)[["variant", "readout", "z_emp", "P_emp", "loo_held", "noloci_z_emp"]].to_string())


if __name__ == "__main__":
    main()
