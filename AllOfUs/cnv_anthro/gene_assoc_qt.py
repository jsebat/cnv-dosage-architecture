"""Gene-based CNV association for a quantitative trait (Jon 2026-11-23: "Should add a gene based test so Milad can output gene based sum stats
(using the standard gene/exon intersect rules for dels and dups)").

Per autosomal gene and dosage, one linear model:   trait ~ covariates + carrier
  carrier (0/1 INDICATOR, not a count) = the person has at least one CNV of that dosage hitting the gene, by the standard rules of this
  pipeline (conga_axes_qt.incidence): deletion = overlaps any exon of the gene by >= 1 bp; duplication = contains a full-length transcript.
  Covariates only: no burden or gene-class terms. This is NOT the per-gene z of gene_mirror_qt.py, which is a score z adjusted for the total
  and constrained gene counts (and, in its model B, the signed TWAS count) and has no beta or SE.
Exact OLS per gene from sufficient statistics (covariates removed by QR): beta in trait units (SD for a rank-normalised trait), SE, t, two-sided
P (t distribution), N, carriers, carrier frequency.
No filter on CNV frequency (Jon: "Don't filter now. Just record ... frequency"): every call of the cohort's CNV table is used.
Annotation columns: chr, start, end (exon span on the gene's home chromosome, hg38), s_het (GeneBayes; imputed values blank), HC (s_het > 0.029),
TWAS z of --twas (scaled tissue combination of --panels), fraction of exonic bases in segmental duplications (--segdup_bed), and with
--pos_freq (position_freq.tsv.gz of cnv_freq_filter.py) max_pos_freq = the highest same-dosage per-position carrier frequency of the whole
dataset inside the gene's exon span (Jon: "Just record probe level frequency"); max_pos_n is its carrier count (dropped at export).
Rows: every gene with >= --min_carriers carriers of that dosage (default 1: all counts stay inside the cohort so that the exporter can apply the
small-cell rule across strata, and so that a supplement at a lower floor can still be made from the full table; export_gene_assoc.py keeps
genes with >= 21 carriers, the floor of every cohort (Jon 2026-11-23), and above the small-cell limit).
null model: none (analytic OLS; a single-gene test has no gene annotation to shift).
Output (aggregate): gene_assoc.tsv.gz, run_info.json. The small-cell rule is applied at export (export_gene_assoc.py), not here.
  python gene_assoc_qt.py --pheno ... --trait height_rint --cnv ... --twas ... --exons ... --gene_cov ... [--segdup_bed ...] --out DIR
"""
import argparse
import json
import os
import numpy as np
import pandas as pd
from scipy import stats
from common import read_plink_cnv, design_matrix
from conga_axes_qt import load_genes, incidence, twas_gene_z, SHET_CUT
from conga_grid_qt import segdup_fraction


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pheno", required=True); ap.add_argument("--trait", required=True); ap.add_argument("--cnv", required=True)
    ap.add_argument("--twas", required=True); ap.add_argument("--panels", default="gtex_nonbrain")
    ap.add_argument("--exons", required=True); ap.add_argument("--gene_cov", required=True); ap.add_argument("--segdup_bed", default="")
    ap.add_argument("--pos_freq", default="", help="position_freq.tsv.gz of cnv_freq_filter.py (dataset-level per-position carrier frequency)")
    ap.add_argument("--numeric", default="age,age2,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"); ap.add_argument("--factors", default="sex,batch")
    ap.add_argument("--cnv_filter", default="", help="label of the call filter applied to --cnv (recorded in run_info.json)")
    ap.add_argument("--min_carriers", type=int, default=1); ap.add_argument("--cohort", default="UKBB"); ap.add_argument("--stratum", default="all")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    ph = pd.read_csv(a.pheno, sep="\t", dtype={"IID": str})
    num = [c for c in a.numeric.split(",") if c]; fac = [c for c in a.factors.split(",") if c]
    ph = ph[np.isfinite(ph[a.trait])].dropna(subset=num + fac).reset_index(drop=True); N = len(ph)
    Zc, _ = design_matrix(ph, num, fac); Q, _ = np.linalg.qr(Zc); p = Q.shape[1]
    y = ph[a.trait].to_numpy(float); r = y - Q @ (Q.T @ y); rr = float(r @ r)
    cnv = read_plink_cnv(a.cnv, keep_iids=ph.IID); cnv = cnv[cnv.CHR <= 22].reset_index(drop=True)
    cnv["row"] = pd.Series(np.arange(N), index=ph.IID).loc[cnv.IID].to_numpy()
    U, ex, tx = load_genes(a.exons); M = len(U)
    span = ex.merge(U[["gidx", "c"]], on=["gidx", "c"]).groupby("gidx").agg(start=("start", "min"), end=("end", "max")).reindex(np.arange(M))
    gc = pd.read_csv(a.gene_cov, sep="\t", dtype={"entrez": str}); gc.loc[gc.shet_imputed == 1, "shet"] = np.nan
    shet = U[["entrez"]].merge(gc[["entrez", "shet"]], on="entrez", how="left").shet.to_numpy()
    tz, npan, _ = twas_gene_z(a.twas, a.panels, U, gc)
    sd = segdup_fraction(ex, M, a.segdup_bed) if a.segdup_bed else np.full(M, np.nan)
    ann = pd.DataFrame(dict(entrez=U.entrez.astype(str), symbol=U.symbol, chr=U.c.to_numpy(), start=span.start.to_numpy(), end=span.end.to_numpy(),
                            shet=shet, HC=(np.nan_to_num(shet, nan=-1) > SHET_CUT).astype(int), twas_z=tz, segdup_frac=np.round(sd, 4)))
    out = []; df = N - p - 1
    PF = pd.read_csv(a.pos_freq, sep="\t") if a.pos_freq else None

    def max_pos(dos):
        """Per gene: (max carriers, max frequency) of the dataset's same-dosage position frequency over the gene's exon span."""
        mn = np.full(M, np.nan); mf = np.full(M, np.nan)
        if PF is None:
            return mn, mf
        for c, f in PF[PF.TYPE == dos].groupby("CHR"):
            f = f.sort_values("START"); s = f.START.to_numpy(); e = f.END.to_numpy(); n = f.n_carriers.to_numpy(float); fr = f.freq.to_numpy(float)
            for g in np.flatnonzero((ann.chr == c).to_numpy()):
                lo = np.searchsorted(e, ann.start.iat[g], "left"); hi = np.searchsorted(s, ann.end.iat[g], "right")
                if hi > lo:
                    j = lo + int(np.argmax(n[lo:hi])); mn[g] = n[j]; mf[g] = fr[j]
        return mn, mf
    for dos, typ, feats, mode in (("DEL", 1, ex, "overlap"), ("DUP", 3, tx, "within")):
        G = incidence(cnv[cnv.TYPE == typ], feats, mode, N, M)
        n = np.asarray(G.sum(0)).ravel(); xr = np.asarray(G.T @ r).ravel(); QtX = np.asarray(G.T @ Q)
        v = n - (QtX ** 2).sum(1)                                   # residual sum of squares of the indicator on the covariates
        with np.errstate(invalid="ignore", divide="ignore"):
            beta = xr / v; se = np.sqrt((rr - xr ** 2 / v) / df / v); t = beta / se
        ok = (n >= a.min_carriers) & (v > 1e-9) & np.isfinite(t)
        d = ann[ok].copy(); d.insert(2, "dosage", dos)
        d["N"] = N; d["n_carriers"] = n[ok].astype(int); d["freq"] = n[ok] / N
        d["beta"] = beta[ok]; d["SE"] = se[ok]; d["t"] = t[ok]; d["P"] = 2 * stats.t.sf(np.abs(t[ok]), df)
        mn, mf = max_pos(dos); d["max_pos_freq"] = mf[ok]; d["max_pos_n"] = mn[ok]
        out.append(d)
    R = pd.concat(out, ignore_index=True)
    R.insert(0, "trait", a.trait); R.insert(0, "stratum", a.stratum); R.insert(0, "cohort", a.cohort)
    R.to_csv(os.path.join(a.out, "gene_assoc.tsv.gz"), sep="\t", index=False, float_format="%.6g")
    info = dict(cohort=a.cohort, stratum=a.stratum, trait=a.trait, N=int(N), covariates=num + fac, panels=a.panels, n_panels=npan, genes=int(M),
                min_carriers=a.min_carriers, rows=int(len(R)), rows_DEL=int((R.dosage == "DEL").sum()), rows_DUP=int((R.dosage == "DUP").sum()),
                test="OLS, trait ~ covariates + carrier indicator; DEL = any exon overlap >= 1 bp, DUP = full-length transcript inside the call",
                segdup_bed=os.path.basename(a.segdup_bed) if a.segdup_bed else None, cnv_table=os.path.basename(a.cnv),
                cnv_filter=a.cnv_filter or "none (every call of the cohort's CNV table)")
    json.dump(info, open(os.path.join(a.out, "run_info.json"), "w"), indent=1)
    lam = np.median(R.t ** 2) / stats.chi2.ppf(.5, 1) if len(R) else np.nan
    print(json.dumps(info)); print(f"lambda (median t^2): {lam:.3f}; smallest P: {R.P.min() if len(R) else float('nan'):.3g}")


if __name__ == "__main__":
    main()
