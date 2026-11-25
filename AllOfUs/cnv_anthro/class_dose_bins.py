"""Genome-wide CNV dose response by gene class (dosage paper Fig. 2, Jon 2026-11-22: "show the U shaped burden dose responses
Fig2a,b. Then show the same thing for CNV-TWAS and for constraint").

Same display as Fig. 2A-B (dup_burden.py signed burden): people who carry deletions only, duplications only, or neither (carriers of
both dosages are left out); x = genes deleted (left) / genes duplicated (right). Here the genes are counted by class, and all classes
of one partition enter ONE joint model, because the classes travel together (a large CNV hits TWAS genes and constrained genes alike):
a raw binned curve for one class would inherit the U of the classes it co-occurs with.
  partitions  all         all genes (reproduces Fig. 2A-B on the joint-model scale)
              twas        up = TWAS z > 3 (higher predicted expression, higher trait), down = z < -3, other = the rest (incl. no z)
              constraint  HC = GeneBayes s_het > 0.029, LC = the rest (incl. imputed / missing s_het)
CNV -> gene: DEL any exon by >= 1 bp, DUP the whole transcript (conga_axes_qt.py). TWAS z = scaled tissue combination of --panels.
Model per partition and run (OLS, covariates removed by QR):
  bins     trait ~ covariates + sum over classes of indicators of the class's signed-count bin (0 = reference)
           -> per class and bin: adjusted difference from people with no gene of that class (others held fixed), SE, n
  shape    trait ~ covariates + sum over classes of [s_c, s_c^2], s_c = genes duplicated - genes deleted of class c
           -> linear (mirror-like, monotone) and quadratic (U) term per class; also DEL and DUP slopes per class
Runs: full, and no_loci (every call assigned to one of the recurrent intervals of --loci / --loo_runs removed, as conga_grid_qt.py).
null model: none (analytic OLS SEs); the genome-shift tests of the same classes are the band tests (band_trend_summary.py).
Output: aggregate tables only; cells with fewer than --min_cell people are written without estimates (AoU rule: >= 21).
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

BINS_WIDE = [(1, 2, "1-2"), (3, 5, "3-5"), (6, 10, "6-10"), (11, 20, "11-20"), (21, 10 ** 6, "21+")]   # Fig. 2A-B bins
BINS_NARROW = [(1, 1, "1"), (2, 2, "2"), (3, 4, "3-4"), (5, 9, "5-9"), (10, 10 ** 6, "10+")]          # sparse classes (Fig. 2I-J bins)
PARTS = {"all": [("all", BINS_WIDE)],
         "twas": [("up", BINS_NARROW), ("down", BINS_NARROW), ("other", BINS_WIDE)],
         "constraint": [("HC", BINS_NARROW), ("LC", BINS_WIDE)],
         # Jon 2026-11-22: "the mirror effects in panels D and E [may be] obscured by constrained genes" -> TWAS classes within one
         # constraint stratum; the other stratum enters the same joint model as one class, so every gene is still in the model
         "twas_LC": [("up_LC", BINS_NARROW), ("down_LC", BINS_NARROW), ("other_LC", BINS_WIDE), ("HC", BINS_NARROW)],
         "twas_HC": [("up_HC", BINS_NARROW), ("down_HC", BINS_NARROW), ("other_HC", BINS_NARROW), ("LC", BINS_WIDE)],
         # Jon 2026-11-22: "partition the other way. Top row twas Z > +3 (LC vs HC), bottom row twas Z < -3 (LC vs HC)": all six
         # TWAS x constraint classes in ONE joint model, so the LC and HC lines of a panel come from the same fit
         "twas_x_con": [("up_LC", BINS_NARROW), ("up_HC", BINS_NARROW), ("down_LC", BINS_NARROW), ("down_HC", BINS_NARROW),
                        ("other_LC", BINS_WIDE), ("other_HC", BINS_NARROW)],
         # Jon 2026-11-22: third row, the unsigned count sum(|z| > 3) (two-axis report v3 G1-G2: height's mirror is carried by the unsigned
         # count): genes with |TWAS z| > 3 in either direction, by constraint, with the same 'other' classes in the same joint model
         "abs3_x_con": [("abs3_LC", BINS_NARROW), ("abs3_HC", BINS_NARROW), ("other_LC", BINS_WIDE), ("other_HC", BINS_NARROW)]}


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
    ap.add_argument("--min_cell", type=int, default=21); ap.add_argument("--cohort", default="UKBB"); ap.add_argument("--out", required=True)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
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
    CLS = {"all": np.ones(M, bool), "up": z > 3, "down": z < -3, "other": ~((z > 3) | (z < -3)), "HC": hc, "LC": ~hc}
    for s_, m_ in (("LC", ~hc), ("HC", hc)):
        CLS[f"up_{s_}"] = CLS["up"] & m_; CLS[f"down_{s_}"] = CLS["down"] & m_; CLS[f"other_{s_}"] = CLS["other"] & m_
        CLS[f"abs3_{s_}"] = (CLS["up"] | CLS["down"]) & m_
    Zc, _ = design_matrix(ph, num, fac); y = ph[a.trait].to_numpy(float)
    loci = pd.read_csv(a.loci, sep="\t"); lab = assign_loci(cnv, loci); lr = pd.read_csv(a.loo_runs, sep="\t", keep_default_na=False)
    is_del = cnv.TYPE.to_numpy() == 1; union = np.zeros(len(cnv), bool)
    for r in lr.itertuples():
        md = set(x for x in str(r.del_members).split(",") if x); mu = set(x for x in str(r.dup_members).split(",") if x)
        union |= (is_del & np.isin(lab, list(md))) | (~is_del & np.isin(lab, list(mu)))
    info = dict(cohort=a.cohort, trait=a.trait, N=N, panels=a.panels, n_panels=npan, genes=int(M),
                genes_per_class={k: int(v.sum()) for k, v in CLS.items()}, calls_removed_no_loci=int(union.sum()), min_cell=a.min_cell)
    bins_out, shape_out = [], []
    for rn, sel in [("full", np.ones(len(cnv), bool)), ("no_loci", ~union)]:
        c_ = cnv[sel]
        Gd = incidence(c_[c_.TYPE == 1], ex, "overlap", N, M); Gu = incidence(c_[c_.TYPE == 3], tx, "within", N, M)
        nD = {k: np.asarray(Gd @ v.astype(float)).ravel() for k, v in CLS.items()}
        nU = {k: np.asarray(Gu @ v.astype(float)).ravel() for k, v in CLS.items()}
        keep = (nD["all"] == 0) | (nU["all"] == 0)                     # single-dosage carriers and non-carriers
        info[f"{rn}_N_single"] = int(keep.sum()); info[f"{rn}_DEL_only"] = int((keep & (nD["all"] > 0)).sum())
        info[f"{rn}_DUP_only"] = int((keep & (nU["all"] > 0)).sum())
        Q, _ = np.linalg.qr(Zc[keep]); yr = y[keep] - Q @ (Q.T @ y[keep])
        for part, cls in PARTS.items():
            cols, meta = [], []
            for c, B in cls:
                for side, cnt, sg in (("DEL", nD[c][keep], -1), ("DUP", nU[c][keep], 1)):
                    other = nU[c][keep] if side == "DEL" else nD[c][keep]
                    for i, (lo, hi, bl) in enumerate(B, start=1):
                        ind = ((cnt >= lo) & (cnt <= hi) & (other == 0)).astype(float)
                        cols.append(ind); meta.append(dict(cls=c, side=side, bin=bl, order=sg * i, n=int(ind.sum()),
                                                           mean_genes=float(cnt[ind > 0].mean()) if ind.sum() else np.nan))
            X = np.column_stack(cols); live = X.sum(0) > 0
            b = np.full(X.shape[1], np.nan); se = b.copy(); p = b.copy()
            b[live], se[live], p[live] = ols(yr, X[:, live], Q)
            for j, m in enumerate(meta):
                ok = m["n"] >= a.min_cell
                bins_out.append({**dict(cohort=a.cohort, trait=a.trait, run=rn, partition=part), **m, "n": m["n"] if ok else f"<{a.min_cell}",
                                 "mean_genes": m["mean_genes"] if ok else np.nan,
                                 "beta": b[j] if ok else np.nan, "SE": se[j] if ok else np.nan, "P": p[j] if ok else np.nan})
            for c, _ in cls:
                n0 = int(((nD[c][keep] == 0) & (nU[c][keep] == 0)).sum())
                bins_out.append(dict(cohort=a.cohort, trait=a.trait, run=rn, partition=part, cls=c, side="none", bin="0", order=0,
                                     n=n0 if n0 >= a.min_cell else f"<{a.min_cell}", mean_genes=0.0, beta=0.0, SE=0.0, P=np.nan))
            # shape terms: per class [s, s^2] jointly; and DEL / DUP slopes jointly (genes per 10 for the wide classes)
            S, names = [], []
            for c, B in cls:
                sc = 10.0 if B is BINS_WIDE else 1.0
                s = (nU[c][keep] - nD[c][keep]) / sc; S += [s, s ** 2]; names += [(c, "linear (signed count)", sc), (c, "quadratic (U)", sc)]
            bb, ss, pp = ols(yr, np.column_stack(S), Q)
            for (c, term, sc), b_, s_, p_ in zip(names, bb, ss, pp):
                shape_out.append(dict(cohort=a.cohort, trait=a.trait, run=rn, partition=part, cls=c, term=term, genes_per_unit=sc, beta=b_, SE=s_, P=p_, N=int(keep.sum())))
            S, names = [], []
            for c, B in cls:
                sc = 10.0 if B is BINS_WIDE else 1.0
                S += [nD[c][keep] / sc, nU[c][keep] / sc]; names += [(c, "DEL slope", sc), (c, "DUP slope", sc)]
            bb, ss, pp = ols(yr, np.column_stack(S), Q)
            for (c, term, sc), b_, s_, p_ in zip(names, bb, ss, pp):
                shape_out.append(dict(cohort=a.cohort, trait=a.trait, run=rn, partition=part, cls=c, term=term, genes_per_unit=sc, beta=b_, SE=s_, P=p_, N=int(keep.sum())))
        print(f"run {rn} done", flush=True)
    pd.DataFrame(bins_out).to_csv(os.path.join(a.out, "class_dose_bins.tsv"), sep="\t", index=False, float_format="%.6g")
    pd.DataFrame(shape_out).to_csv(os.path.join(a.out, "class_dose_shape.tsv"), sep="\t", index=False, float_format="%.6g")
    json.dump(info, open(os.path.join(a.out, "run_info.json"), "w"), indent=1)
    pd.set_option("display.width", 220); print(json.dumps(info)); print(pd.DataFrame(shape_out).to_string())


if __name__ == "__main__":
    main()
