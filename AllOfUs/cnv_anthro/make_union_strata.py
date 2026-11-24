"""Phenotype subsets that exclude carriers of a CNV hitting a single-hit ('union') gene (Jon 2026-11-15).

Union gene = top-5% GeneBayes s_het or ClinGen haploinsufficiency score 2-3 (same as dup_burden.py / additivity_constraint.py).
A deletion hits a gene if it overlaps any exon; a duplication if it spans a whole transcript (conga_axes_qt.incidence).
Writes, next to --out (Expanse only; person-level):
  pheno_noUnionDUP.tsv   people with no duplication of a union gene
  pheno_noUnionAny.tsv   people with no deletion or duplication of a union gene
and union_strata_info.json (counts only).
"""
import argparse, json, os
import numpy as np, pandas as pd
from common import read_plink_cnv, rcg
from conga_axes_qt import load_genes, incidence


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pheno", required=True); ap.add_argument("--cnv", required=True); ap.add_argument("--exons", required=True)
    ap.add_argument("--gene_cov", required=True); ap.add_argument("--clingen", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--asym_loci", default=None, help="TSV NAME CHR START END of the asymmetric intervals (S2/S3 strata)")
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    ph = pd.read_csv(a.pheno, sep="\t", dtype={"IID": str}); N = len(ph)
    cnv = read_plink_cnv(a.cnv, keep_iids=ph.IID); cnv = cnv[cnv.CHR <= 22].reset_index(drop=True)
    cnv["row"] = pd.Series(np.arange(N), index=ph.IID).loc[cnv.IID].to_numpy()
    U, ex, tx = load_genes(a.exons); M = len(U)
    gc = pd.read_csv(a.gene_cov, sep="\t", dtype={"entrez": str}); gc.loc[gc.shet_imputed == 1, "shet"] = np.nan
    shet = U[["entrez"]].merge(gc[["entrez", "shet"]], on="entrez", how="left").shet.to_numpy()
    q95 = float(np.nanquantile(gc.shet.to_numpy(float), 0.95))
    cg = pd.read_csv(a.clingen, sep="\t", skiprows=5, dtype=str, quoting=3); cg.columns = [c.strip().lstrip("#") for c in cg.columns]
    hi = U.symbol.isin(set(cg.loc[cg["Haploinsufficiency Score"].isin(["2", "3"]), cg.columns[0]])).to_numpy()
    w = ((np.nan_to_num(shet, nan=-1) >= q95) | hi).astype(float)
    Gd = incidence(cnv[cnv.TYPE == 1], ex, "overlap", N, M); Gu = incidence(cnv[cnv.TYPE == 3], tx, "within", N, M)
    hDEL = np.asarray(Gd @ w).ravel() > 0; hDUP = np.asarray(Gu @ w).ravel() > 0
    ph[~hDUP].to_csv(os.path.join(a.out, "pheno_noUnionDUP.tsv"), sep="\t", index=False)
    ph[~(hDEL | hDUP)].to_csv(os.path.join(a.out, "pheno_noUnionAny.tsv"), sep="\t", index=False)
    info = dict(N=int(N), union_genes=int(w.sum()), q95_shet=q95, union_DEL_carriers=int(hDEL.sum()), union_DUP_carriers=int(hDUP.sum()),
                union_any_carriers=int((hDEL | hDUP).sum()), N_noUnionDUP=int((~hDUP).sum()), N_noUnionAny=int((~(hDEL | hDUP)).sum()))
    # S0 (Jon 2026-11-15): all people, union presence per arm as fixed person-level covariates (observed annotation, same in every null replicate)
    p0 = ph.copy(); p0["unionDEL"] = hDEL.astype(float); p0["unionDUP"] = hDUP.astype(float)
    p0.to_csv(os.path.join(a.out, "pheno_withUnionCov.tsv"), sep="\t", index=False)
    if a.asym_loci:   # S2/S3: carriers of Molly's asymmetric intervals (either dosage; Molly's recurrent CNV genotyping), as in dup_burden.py
        asym = set()
        asym = rcg.carriers(rcg.genotype_calls(cnv), pd.read_csv(a.asym_loci, sep="\t").NAME, "both")   # Molly's genotyping (Jon 2026-11-24)
        am = ph.IID.isin(asym).to_numpy()
        ph[~am].to_csv(os.path.join(a.out, "pheno_noAsym13.tsv"), sep="\t", index=False)
        ph[~am & ~hDUP].to_csv(os.path.join(a.out, "pheno_noUnionDUP_noAsym13.tsv"), sep="\t", index=False)
        info.update(asym13_carriers=int(am.sum()), N_noAsym13=int((~am).sum()), N_noUnionDUP_noAsym13=int((~am & ~hDUP).sum()),
                    asym13_and_unionDUP=int((am & hDUP).sum()))
    json.dump(info, open(os.path.join(a.out, "union_strata_info.json"), "w"), indent=1); print(json.dumps(info))


if __name__ == "__main__":
    main()
