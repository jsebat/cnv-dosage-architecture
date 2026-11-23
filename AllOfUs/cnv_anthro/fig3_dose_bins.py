"""Binned dose response for every stratum of the Fig. 3 joint model (Jon 2026-11-23: "make sure that we can make dose response plots like this
for ALL of the major strata in new joint model (Fig. 3) ... DUP z and DEL z for all CNVs, Single Gene, >2 Gene; z>+3, z<-3, weak or no TWAS, total
genes"). The display of the dosage paper's Fig. 2a-h (class_dose_bins.py), extended by CNV size.

People who carry deletions only, duplications only, or neither (carriers of gene-hitting CNVs of both dosages are left out, as Fig. 2); per gene
class and CNV-size stratum, the number of genes deleted / duplicated is cut into FIXED bins (identical in every cohort). Each bin: adjusted
difference in the trait (SD) from people with no gene of that class in that CNV-size stratum; SE; carriers; mean genes.
  gene classes   up_LC, up_HC (TWAS z > +3), down_LC, down_HC (z < -3), other_LC, other_HC (weak or no TWAS: -3 <= z <= 3 or no z), and
                 'all' (total genes hit); LC = s_het <= 0.029 or missing, HC = s_het > 0.029 (as conga_grid_qt.py cls:*)
  CNV size       N1_inf = all CNVs; N1_1 = single-gene CNVs; N2_inf = CNVs hitting >= 2 genes (genes of the call's own dosage rule)
  partitions     uni         UNIVARIATE (Jon 2026-11-23: "This comes before the joint model, so these are univariate tests of stratified
                             CNVs"): one OLS per class and CNV-size stratum, trait ~ covariates + that class's DEL and DUP bins, no other
                             class terms. 7 classes x 3 CNV sizes. Every class has its own curve, weak-or-no-TWAS LC included.
                 all         total genes, all CNVs (= uni / N1_inf / all; = class_dose_bins.py partition 'all', Fig. 2a-b)
                 twas_x_con  the six classes, all CNVs, in ONE model (= class_dose_bins.py partition 'twas_x_con', Fig. 2c-h); kept so the
                             Fig. 2 tables can be reproduced and compared with the univariate curves
  shape terms per class and size (within the same model): linear and quadratic term of the signed count, and DEL / DUP slopes.
  The univariate DEL z and DUP z of each class and stratum against the genome-shift null are NOT computed here: they are the univariate
  variants of conga_grid_qt.py (designs/fig3_univariate_2026-11-23.json, run_fig3.py fits uni_*).
Runs: full, and no_loci (calls assigned to the 34 recurrent intervals removed, as conga_grid_qt.py).
CNV -> gene: DEL any exon by >= 1 bp, DUP the whole transcript. TWAS z = scaled tissue combination of --panels.
null model: none here (analytic OLS SEs: the bin points and 95% CIs are descriptive); the z of each class is tested against the genome shift
of the gene annotation by the univariate conga_grid_qt.py fits (2,000 shared shifts, 9-run LOO, no_loci); total genes: analytic t.
Output (aggregate): fig3_dose_bins.tsv, fig3_dose_shape.tsv, run_info.json. Small-cell rule (--min_cell k): a bin with 1..k people is written
without n, mean_genes, beta, SE, P; when k > 0 the no_loci run carries no n at all (n differences between the two runs could reveal a small
count; the meta-analysis needs only beta and SE).
"""
import argparse
import json
import os
import numpy as np
import pandas as pd
from scipy import stats
from common import read_plink_cnv, design_matrix
from conga_axes_qt import load_genes, incidence, assign_loci, SHET_CUT
from conga_grid_qt import twas_z_both

BINS_WIDE = [(1, 2, "1-2"), (3, 5, "3-5"), (6, 10, "6-10"), (11, 20, "11-20"), (21, 10 ** 6, "21+")]   # as class_dose_bins.py
BINS_NARROW = [(1, 1, "1"), (2, 2, "2"), (3, 4, "3-4"), (5, 9, "5-9"), (10, 10 ** 6, "10+")]
CLS6 = [("up_LC", BINS_NARROW), ("up_HC", BINS_NARROW), ("down_LC", BINS_NARROW), ("down_HC", BINS_NARROW),
        ("other_LC", BINS_WIDE), ("other_HC", BINS_NARROW)]
SIZES = {"N1_inf": lambda ng: ng >= 1, "N1_1": lambda ng: ng == 1, "N2_inf": lambda ng: ng >= 2}


# (partition, [(size, class, bins), ...]): the terms of one OLS
MODELS = [("all", [("N1_inf", "all", BINS_WIDE)]), ("twas_x_con", [("N1_inf", c, B) for c, B in CLS6])] +          [("uni", [(sz, c, B)]) for sz in SIZES for c, B in CLS6 + [("all", BINS_WIDE)]]


def ols(yr, X, Q):
    Xr = X - Q @ (Q.T @ X); A = Xr.T @ Xr; Ai = np.linalg.pinv(A); b = Ai @ (Xr.T @ yr)
    df = len(yr) - Q.shape[1] - X.shape[1]; s2 = (yr @ yr - b @ (Xr.T @ yr)) / df; se = np.sqrt(np.diag(Ai) * s2)
    return b, se, 2 * stats.t.sf(np.abs(b / se), df)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pheno", required=True); ap.add_argument("--trait", required=True); ap.add_argument("--cnv", required=True)
    ap.add_argument("--twas", required=True); ap.add_argument("--panels", default="gtex_nonbrain")
    ap.add_argument("--exons", required=True); ap.add_argument("--gene_cov", required=True)
    ap.add_argument("--loci", required=True); ap.add_argument("--loo_runs", required=True)
    ap.add_argument("--numeric", default="age,age2,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"); ap.add_argument("--factors", default="sex,batch")
    ap.add_argument("--min_cell", type=int, default=20, help="bins with 1..k people are withheld"); ap.add_argument("--cohort", default="UKBB")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); k = a.min_cell
    ph = pd.read_csv(a.pheno, sep="\t", dtype={"IID": str})
    num = [c for c in a.numeric.split(",") if c]; fac = [c for c in a.factors.split(",") if c]
    ph = ph[np.isfinite(ph[a.trait])].dropna(subset=num + fac).reset_index(drop=True); N = len(ph)
    cnv = read_plink_cnv(a.cnv, keep_iids=ph.IID); cnv = cnv[cnv.CHR <= 22].reset_index(drop=True)
    cnv["row"] = pd.Series(np.arange(N), index=ph.IID).loc[cnv.IID].to_numpy()
    U, ex, tx = load_genes(a.exons); M = len(U)
    gc = pd.read_csv(a.gene_cov, sep="\t", dtype={"entrez": str}); gc.loc[gc.shet_imputed == 1, "shet"] = np.nan
    shet = U[["entrez"]].merge(gc[["entrez", "shet"]], on="entrez", how="left").shet.to_numpy()
    Z, npan = twas_z_both(a.twas, a.panels, U, gc); z = np.nan_to_num(Z["scaled"], nan=0.0)
    hc = np.nan_to_num(shet, nan=-1) > SHET_CUT
    up, dn = z > 3, z < -3; oth = ~(up | dn)
    CLS = {"all": np.ones(M, bool)}
    for s_, m_ in (("LC", ~hc), ("HC", hc)):
        CLS[f"up_{s_}"] = up & m_; CLS[f"down_{s_}"] = dn & m_; CLS[f"other_{s_}"] = oth & m_
    Zc, _ = design_matrix(ph, num, fac); y = ph[a.trait].to_numpy(float)
    loci = pd.read_csv(a.loci, sep="\t"); lab = assign_loci(cnv, loci); lr = pd.read_csv(a.loo_runs, sep="\t", keep_default_na=False)
    is_del = cnv.TYPE.to_numpy() == 1; union = np.zeros(len(cnv), bool)
    for r in lr.itertuples():
        md = set(x for x in str(r.del_members).split(",") if x); mu = set(x for x in str(r.dup_members).split(",") if x)
        union |= (is_del & np.isin(lab, list(md))) | (~is_del & np.isin(lab, list(mu)))
    # genes hit per call (its own dosage rule): the CNV-size class of a call, as conga_grid_qt.py
    ng = np.zeros(len(cnv)); cc_ = cnv.copy(); cc_["row"] = np.arange(len(cc_))
    for typ, feats, mode in [(1, ex, "overlap"), (3, tx, "within")]:
        m = (cnv.TYPE == typ).to_numpy()
        ng[m] = incidence(cc_[m], feats, mode, len(cnv), M).getnnz(axis=1)[m]
    info = dict(cohort=a.cohort, trait=a.trait, N=N, panels=a.panels, n_panels=npan, genes=int(M),
                genes_per_class={c: int(v.sum()) for c, v in CLS.items()}, min_cell=k,
                bins=dict(wide=[b[2] for b in BINS_WIDE], narrow=[b[2] for b in BINS_NARROW]), models=sorted({m[0] for m in MODELS}), sizes=list(SIZES))
    mask = lambda n: (1 <= n <= k)
    bins_out, shape_out = [], []
    for rn, sel in [("full", np.ones(len(cnv), bool)), ("no_loci", ~union)]:
        nD, nU = {}, {}
        for sz, f in SIZES.items():
            c_ = cnv[sel & f(ng)]
            Gd = incidence(c_[c_.TYPE == 1], ex, "overlap", N, M); Gu = incidence(c_[c_.TYPE == 3], tx, "within", N, M)
            for c, v in CLS.items():
                nD[(sz, c)] = np.asarray(Gd @ v.astype(float)).ravel(); nU[(sz, c)] = np.asarray(Gu @ v.astype(float)).ravel()
        keep = (nD[("N1_inf", "all")] == 0) | (nU[("N1_inf", "all")] == 0)          # single-dosage carriers and non-carriers
        if rn == "full" or k == 0:
            info[f"{rn}_N_single"] = int(keep.sum())
        Q, _ = np.linalg.qr(Zc[keep]); yr = y[keep] - Q @ (Q.T @ y[keep])
        for part, terms in MODELS:
            cols, meta = [], []
            if True:
                for sz, c, B in terms:
                    for side, cnt, sg in (("DEL", nD[(sz, c)][keep], -1), ("DUP", nU[(sz, c)][keep], 1)):
                        other = nU[(sz, c)][keep] if side == "DEL" else nD[(sz, c)][keep]
                        for i, (lo, hi, bl) in enumerate(B, start=1):
                            ind = ((cnt >= lo) & (cnt <= hi) & (other == 0)).astype(float)
                            cols.append(ind); meta.append(dict(size=sz, cls=c, side=side, bin=bl, order=sg * i, n=int(ind.sum()),
                                                               mean_genes=float(cnt[ind > 0].mean()) if ind.sum() else np.nan))
            X = np.column_stack(cols); live = X.sum(0) > 0
            b = np.full(X.shape[1], np.nan); se = b.copy(); p = b.copy()
            b[live], se[live], p[live] = ols(yr, X[:, live], Q)
            for j, m in enumerate(meta):
                ok = m["n"] > 0 and not mask(m["n"])
                nn = (m["n"] if ok else (f"<={k}" if mask(m["n"]) else 0)) if (rn == "full" or k == 0) else np.nan
                bins_out.append({**dict(cohort=a.cohort, trait=a.trait, run=rn, partition=part), **m, "n": nn,
                                 "mean_genes": m["mean_genes"] if ok else np.nan, "beta": b[j] if ok else np.nan, "SE": se[j] if ok else np.nan,
                                 "P": p[j] if ok else np.nan})
            if True:
                for sz, c, _ in terms:
                    n0 = int(((nD[(sz, c)][keep] == 0) & (nU[(sz, c)][keep] == 0)).sum())
                    nn = (n0 if not mask(n0) else f"<={k}") if (rn == "full" or k == 0) else np.nan
                    bins_out.append(dict(cohort=a.cohort, trait=a.trait, run=rn, partition=part, size=sz, cls=c, side="none", bin="0", order=0,
                                         n=nn, mean_genes=0.0, beta=0.0, SE=0.0, P=np.nan))
            for kind in ("signed", "slopes"):
                S, names = [], []
                if True:
                    for sz, c, B in terms:
                        sc = 10.0 if B is BINS_WIDE else 1.0; d_, u_ = nD[(sz, c)][keep], nU[(sz, c)][keep]
                        if kind == "signed":
                            s = (u_ - d_) / sc; S += [s, s ** 2]; names += [(sz, c, "linear (signed count)", sc), (sz, c, "quadratic (U)", sc)]
                        else:
                            S += [d_ / sc, u_ / sc]; names += [(sz, c, "DEL slope", sc), (sz, c, "DUP slope", sc)]
                S = np.column_stack(S); nz = (S != 0).sum(0)
                bb, ss, pp = ols(yr, S, Q)
                for (sz, c, term, sc), b_, s_, p_, n_ in zip(names, bb, ss, pp, nz):
                    ok = n_ > 0 and not mask(int(n_))
                    shape_out.append(dict(cohort=a.cohort, trait=a.trait, run=rn, partition=part, size=sz, cls=c, term=term, genes_per_unit=sc,
                                          beta=b_ if ok else np.nan, SE=s_ if ok else np.nan, P=p_ if ok else np.nan, N=int(keep.sum()) if (rn == "full" or k == 0) else np.nan))
        print(f"run {rn} done", flush=True)
    B = pd.DataFrame(bins_out); B.to_csv(os.path.join(a.out, "fig3_dose_bins.tsv"), sep="\t", index=False, float_format="%.6g")
    pd.DataFrame(shape_out).to_csv(os.path.join(a.out, "fig3_dose_shape.tsv"), sep="\t", index=False, float_format="%.6g")
    json.dump(info, open(os.path.join(a.out, "run_info.json"), "w"), indent=1)
    print(json.dumps(info)); print("bins:", len(B), "| withheld (1..k people):", int((B.n.astype(str) == f"<={k}").sum()))


if __name__ == "__main__":
    main()
