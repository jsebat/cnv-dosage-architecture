"""Gene-level TOTAL mirror and the part explained by TWAS, for a quantitative trait (dosage paper Fig. 3; Jon 2026-11-22: "As per the Shanta
paper. The yellow arrow is the total mirror coming from correlation of Del and Dup effects. Orange arrow is the portion explained by GWAS. You
need to add a total Del vs Dup to this like Shanta"). Height / BMI port of conga_gene_mirror.R + gene_mirror_summary.R (Engchuan et al.,
companion manuscript; Shanta et al. Fig. 3I), on the UK Biobank model of conga_grid_qt.py (same phenotype, covariates, CNV-to-gene rules:
deletion = any exon overlap, duplication = full-length transcript inside the call; same loci and LOO runs).

Per run, CNV-size class, dosage and gene g: efficient-score z of carrying a CNV of that dosage that hits g (indicator column of the person x gene
incidence), linear model, from the null fit  trait ~ covariates + total genes deleted + total genes duplicated  adjusted for
  model A: constrained-gene counts (s_het > 0.029) per dosage                                   -> zA[g, dos]
  model B: model A + signed TWAS count per dosage, n(z > 3) - n(z < -3) (own-trait TWAS)         -> zB[g, dos]
Statistics (pre-specified as in Shanta), over genes with >= 5 deletion AND >= 5 duplication carriers in the run (and in the gene stratum:
all genes, less constrained, constrained):
  total     = Pearson r(zA[, DEL], zA[, DUP])     (negative = mirror line DUP = -DEL)
  residual  = Pearson r(zB[, DEL], zB[, DUP])     (not explained by TWAS)
  explained = total - residual; TWAS share = explained / total
null model: DUP-vs-DEL GENE-CIRCLE SHIFT (as gene_mirror_summary.R). Autosomal genes in genome order; the DUP z vector is shifted k genes along the
  circle against the DEL z vector, k >= KMIN = (largest number of genes in any 10-Mb window) + 1, so no gene is re-paired with a gene of its own
  locus; r recomputed per replicate (10,000). KEEPS each dosage's own per-gene profile (clustering of neighbouring genes that share CNVs);
  REMOVES only which DEL effect is paired with which DUP effect. The same k set for every model, stratum, size class and run (paired, so
  'explained' has its own null). NOT the segment flip, NOT the genomic shift of annotations (neither moves r). LOO: each major-locus run and
  minor_all against the shift null of the same run (matched); held = same sign and |z| > 1.96; no_loci = all 34 recurrent intervals removed.
Outputs (aggregate only): gene_mirror_summary.tsv; gene_z.npz (per-gene zA/zB and carrier counts per run, for meta-analysis: cohorts combine the
per-gene z by sqrt(N)-weighted Stouffer and recompute r); run_info.json.
  python gene_mirror_qt.py --pheno ... --trait height_rint --cnv ... --twas ... --exons ... --gene_cov ... --loci ... --loo_runs ...
         --major_loci 16p11_2_BP4_BP5,... --out DIR
"""
import argparse
import json
import os
import time
import numpy as np
import pandas as pd
from scipy import sparse
from common import read_plink_cnv, design_matrix
from conga_axes_qt import load_genes, incidence, assign_loci, twas_gene_z, SHET_CUT, Z_CUT

MINCAR = 5


def gene_scores(G, Q, y):
    """Efficient-score z per gene column of G (N x M sparse 0/1) for y ~ span(Q) (Q orthonormal, N x p): z = x'r / sqrt(s2 (x'x - |Q'x|^2))."""
    r = y - Q @ (Q.T @ y); s2 = float(r @ r) / (len(y) - Q.shape[1])
    xr = np.asarray(G.T @ r).ravel(); xx = np.asarray(G.sum(0)).ravel()
    QtX = np.asarray((G.T @ Q))                     # M x p
    v = xx - (QtX ** 2).sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        z = xr / np.sqrt(s2 * v)
    z[(xx == 0) | ~(v > 1e-9)] = np.nan
    return z, xx


def augment(Q, F):
    """Orthonormal basis of [Q, F]."""
    F = F - Q @ (Q.T @ F)
    Q2, R = np.linalg.qr(F)
    keep = np.abs(np.diag(R)) > 1e-8 * max(1.0, np.abs(np.diag(R)).max())
    return np.hstack([Q, Q2[:, keep]])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pheno", required=True); ap.add_argument("--trait", required=True); ap.add_argument("--cnv", required=True)
    ap.add_argument("--twas", required=True); ap.add_argument("--panels", default="gtex_nonbrain")
    ap.add_argument("--exons", required=True); ap.add_argument("--gene_cov", required=True)
    ap.add_argument("--loci", required=True); ap.add_argument("--loo_runs", required=True)
    ap.add_argument("--major_loci", required=True, help="comma list of major LOO families (each its own matched run); all others dropped together (minor_all)")
    ap.add_argument("--numeric", default="age,age2,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10")
    ap.add_argument("--factors", default="sex,batch")
    ap.add_argument("--nrep", type=int, default=10000); ap.add_argument("--seed", type=int, default=20261122)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); t0 = time.time()

    ph = pd.read_csv(a.pheno, sep="\t", dtype={"IID": str})
    num = [c for c in a.numeric.split(",") if c]; fac = [c for c in a.factors.split(",") if c]
    ph = ph[np.isfinite(ph[a.trait])].dropna(subset=num + fac).reset_index(drop=True)
    Zc, _ = design_matrix(ph, num, fac); Q0, _ = np.linalg.qr(Zc)
    y = ph[a.trait].to_numpy(float); N = len(y)
    cnv = read_plink_cnv(a.cnv, keep_iids=ph.IID); cnv = cnv[cnv.CHR <= 22].reset_index(drop=True)
    cnv["row"] = pd.Series(np.arange(N), index=ph.IID).loc[cnv.IID].to_numpy()
    U, ex, tx = load_genes(a.exons); M = len(U)
    gc = pd.read_csv(a.gene_cov, sep="\t", dtype={"entrez": str}); gc.loc[gc.shet_imputed == 1, "shet"] = np.nan
    shet = U[["entrez"]].merge(gc[["entrez", "shet"]], on="entrez", how="left").shet.to_numpy()
    hc = np.nan_to_num(shet, nan=-1) > SHET_CUT
    tz, npan, _ = twas_gene_z(a.twas, a.panels, U, gc); tzf = np.nan_to_num(tz, nan=0.0)
    sgn = (tzf > Z_CUT).astype(float) - (tzf < -Z_CUT).astype(float)

    # gene circle and KMIN (genes in genome order = U order)
    ch = U.c.to_numpy(); pos = U.pos.to_numpy()
    win = np.zeros(M, int)
    for c in np.unique(ch):
        idx = np.flatnonzero(ch == c); p = pos[idx]
        win[idx] = np.searchsorted(p, p + 1e7, "left") - np.arange(len(p))
    KMIN = int(win.max()) + 1
    rng = np.random.default_rng(a.seed); KS = rng.integers(KMIN, M - KMIN + 1, a.nrep)

    # runs (as conga_grid_qt.py --loo_mode drivers --major_loci)
    loci = pd.read_csv(a.loci, sep="\t"); lab = assign_loci(cnv, loci)
    lr = pd.read_csv(a.loo_runs, sep="\t", keep_default_na=False)
    is_del = cnv.TYPE.to_numpy() == 1; drops = {}; union = np.zeros(len(cnv), bool)
    for r in lr.itertuples():
        md = set(x for x in str(r.del_members).split(",") if x); mu = set(x for x in str(r.dup_members).split(",") if x)
        d = (is_del & np.isin(lab, list(md))) | (~is_del & np.isin(lab, list(mu))); drops[r.run] = d; union |= d
    majors = [x for x in a.major_loci.split(",") if x in drops]
    minor = np.zeros(len(cnv), bool)
    for k_, d in drops.items():
        if k_ not in majors:
            minor |= d
    runs = {"full": np.ones(len(cnv), bool)}
    runs.update({f"loo_{m}": ~drops[m] for m in majors}); runs["minor_all"] = ~minor; runs["no_loci"] = ~union
    loo = [r for r in runs if r.startswith("loo_")] + ["minor_all"]

    # per-call gene count (size class), same CNV-to-gene rule
    ng = np.zeros(len(cnv)); cc_ = cnv.copy(); cc_["row"] = np.arange(len(cc_))
    for typ, feats, mode in [(1, ex, "overlap"), (3, tx, "within")]:
        m = (cnv.TYPE == typ).to_numpy()
        ng[m] = incidence(cc_[m], feats, mode, len(cnv), M).getnnz(axis=1)[m]
    SIZES = {"all sizes": np.ones(len(cnv), bool), "N1_1": ng == 1}

    ZA, ZB, NC = {}, {}, {}
    for rn, sel in runs.items():
        c_ = cnv[sel]
        Gd = incidence(c_[c_.TYPE == 1], ex, "overlap", N, M); Gu = incidence(c_[c_.TYPE == 3], tx, "within", N, M)
        # covariates use ALL CNVs of the run (every size), as the joint model; the tested indicator is per size class
        F_A = np.column_stack([np.asarray(Gd.sum(1)).ravel(), np.asarray(Gu.sum(1)).ravel(), Gd @ hc.astype(float), Gu @ hc.astype(float)])
        F_B = np.column_stack([F_A, Gd @ sgn, Gu @ sgn])
        QA = augment(Q0, F_A); QB = augment(Q0, F_B)
        for sz, sm in SIZES.items():
            cs = c_[sm[sel]]
            Hd = Gd if sz == "all sizes" else incidence(cs[cs.TYPE == 1], ex, "overlap", N, M)
            Hu = Gu if sz == "all sizes" else incidence(cs[cs.TYPE == 3], tx, "within", N, M)
            for dos, H in (("DEL", Hd), ("DUP", Hu)):
                ZA[(rn, sz, dos)], NC[(rn, sz, dos)] = gene_scores(H, QA, y)
                ZB[(rn, sz, dos)], _ = gene_scores(H, QB, y)
        print(f"[{time.time() - t0:7.1f}s] run {rn} done", flush=True)

    STRATA = {"all": np.ones(M, bool), "LC": ~hc, "HC": hc}
    idx0 = np.arange(M)

    def rr(zd, zu, cd, cu, keep, k=0):
        j = (idx0 + k) % M
        b = zu[j]; ok = keep & (cd >= MINCAR) & (cu[j] >= MINCAR) & np.isfinite(zd) & np.isfinite(b)
        return (np.corrcoef(zd[ok], b[ok])[0, 1], int(ok.sum())) if ok.sum() > 2 else (np.nan, int(ok.sum()))

    res = {}
    for rn in runs:
        for sz in SIZES:
            cd, cu = NC[(rn, sz, "DEL")], NC[(rn, sz, "DUP")]
            for st, keep in STRATA.items():
                o = {}; nul = {}
                for mdl, Zm in (("A", ZA), ("B", ZB)):
                    zd, zu = Zm[(rn, sz, "DEL")], Zm[(rn, sz, "DUP")]
                    o[mdl], ngn = rr(zd, zu, cd, cu, keep)
                    nul[mdl] = np.array([rr(zd, zu, cd, cu, keep, int(k))[0] for k in KS])
                res[(rn, sz, st)] = dict(n_genes=ngn, total=(o["A"], nul["A"]), residual=(o["B"], nul["B"]),
                                         explained=(o["A"] - o["B"], nul["A"] - nul["B"]))
        print(f"[{time.time() - t0:7.1f}s] null {rn} done", flush=True)

    def emp(o, nu):
        nu = nu[np.isfinite(nu)]
        if not np.isfinite(o) or len(nu) < 10:
            return np.nan, np.nan
        mu, sd = nu.mean(), nu.std(ddof=1)
        return (o - mu) / sd, (1 + np.sum(np.abs(nu - mu) >= abs(o - mu))) / (len(nu) + 1)

    rows = []
    for sz in SIZES:
        for st in STRATA:
            R = res[("full", sz, st)]
            for stat in ("total", "residual", "explained"):
                o, nu = R[stat]; z, p = emp(o, nu); s = np.sign(z)
                lz = {rn: emp(*res[(rn, sz, st)][stat])[0] for rn in loo}
                held = sum(1 for q in lz.values() if np.isfinite(q) and abs(q) > 1.96 and np.sign(q) == s)
                fin = {k_: q for k_, q in lz.items() if np.isfinite(q)}
                wk = min(fin, key=lambda k_: s * fin[k_]) if fin else ""
                zn, pn = emp(*res[("no_loci", sz, st)][stat])
                tot = R["total"][0]
                rows.append(dict(trait=a.trait, size=sz, stratum=st, stat=stat, obs=o, null_mean=np.nanmean(nu), z=z, P=p,
                                 n_genes=R["n_genes"], share=(o / tot if stat == "explained" and np.isfinite(tot) and tot != 0 else np.nan),
                                 loo_held=held, n_loo=len(loo), weakest=wk.replace("loo_", ""), z_weakest=fin.get(wk, np.nan),
                                 noloci_obs=res[("no_loci", sz, st)][stat][0], noloci_z=zn, noloci_P=pn, n_genes_noloci=res[("no_loci", sz, st)]["n_genes"]))
    out = pd.DataFrame(rows); out.to_csv(os.path.join(a.out, "gene_mirror_summary.tsv"), sep="\t", index=False, float_format="%.6g")
    np.savez_compressed(os.path.join(a.out, "gene_z.npz"), entrez=U.entrez.to_numpy().astype(str), hc=hc, twas_sign=sgn, KS=KS,
                        **{f"zA__{rn}__{sz}__{d}": ZA[(rn, sz, d)] for (rn, sz, d) in ZA},
                        **{f"zB__{rn}__{sz}__{d}": ZB[(rn, sz, d)] for (rn, sz, d) in ZB},
                        **{f"n__{rn}__{sz}__{d}": NC[(rn, sz, d)] for (rn, sz, d) in NC})
    json.dump(dict(trait=a.trait, N=int(N), genes=int(M), KMIN=KMIN, nrep=a.nrep, seed=a.seed, panels=a.panels, n_panels=npan,
                   major_loci=majors, runs=list(runs), mincar=MINCAR, minutes=(time.time() - t0) / 60),
              open(os.path.join(a.out, "run_info.json"), "w"), indent=1)
    pd.set_option("display.width", 220); print(out.drop(columns=["trait"]).to_string())


if __name__ == "__main__":
    main()
