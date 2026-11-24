"""Per-locus DEL and DUP effects on quantitative traits (aggregate output only).

For each locus (NAME CHR START END [RULE]) and each dosage, carriers are people with a call of that type matching the locus by RULE:
  recip50  - recurrent CNV loci: carriers from Molly's recurrent CNV genotyping (recurrentCNVgenotyping/aou, rcg.carriers;
             reciprocal / one-way sheet rules and segment patterns). No overlap rule here (Jon 2026-11-24).
  cover50  - the call covers >= 50% of the locus interval (probe-group loci from the CNV GWAS; zero-length loci = call covers the point)
The effect is the exact OLS coefficient of carrier status in  trait ~ covariates + carrier  (same model as cnv_gwas_qt.py / alleles.py).
Usage: python locus_effects.py --pheno P --cnv C --loci L.tsv --out O.tsv [--traits height_rint,bmi_rint] [--numeric ...] [--factors ...]
"""
import argparse
import json
import numpy as np
import pandas as pd
from common import read_plink_cnv, design_matrix, Residualiser, rcg


def carriers(cnv, chrom, s, e, typ, rule, name=None, G=None):
    c = cnv[(cnv.CHR == chrom) & (cnv.TYPE == typ) & (cnv.BP2 >= s) & (cnv.BP1 <= e)]
    ov = (np.minimum(c.BP2, e) - np.maximum(c.BP1, s)).clip(lower=0)
    L = max(e - s, 0)
    if rule == "recip50":
        m = c.IID.isin(rcg.carriers(G, name, "DEL" if typ == 1 else "DUP")).to_numpy()
    elif rule == "cover50":
        m = (ov >= 0.5 * L) if L > 0 else (c.BP1 <= s) & (c.BP2 >= e)
    else:
        raise ValueError(rule)
    c = c[m]
    return c.IID.unique(), (c.BP2 - c.BP1).median() if len(c) else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pheno", required=True); ap.add_argument("--cnv", required=True); ap.add_argument("--loci", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--default_rule", default="recip50")
    ap.add_argument("--traits", default="height_rint,bmi_rint")
    ap.add_argument("--numeric", default="age,age2,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10")
    ap.add_argument("--factors", default="sex,batch")
    ap.add_argument("--min_carriers", type=int, default=12)
    a = ap.parse_args()
    ph = pd.read_csv(a.pheno, sep="\t", dtype={"IID": str})
    cnv = read_plink_cnv(a.cnv, keep_iids=ph.IID)
    loci = pd.read_csv(a.loci, sep="\t")
    G = rcg.genotype_calls(cnv)
    if "RULE" not in loci:
        loci["RULE"] = a.default_rule
    num = [x for x in a.numeric.split(",") if x]; fac = [x for x in a.factors.split(",") if x]
    rows = []
    for tr in a.traits.split(","):
        P = ph[np.isfinite(ph[tr])].reset_index(drop=True)
        pos = pd.Series(np.arange(len(P)), index=P.IID)
        Z, _ = design_matrix(P, num, fac); R = Residualiser(Z)
        yr = R.resid(P[tr].to_numpy(float)); yy = yr @ yr
        for l in loci.itertuples():
            for typ, lab in [(1, "DEL"), (3, "DUP")]:
                ids, med = carriers(cnv, int(l.CHR), int(l.START), int(l.END), typ, l.RULE, l.NAME, G)
                idx = pos.reindex(ids).dropna().astype(int).to_numpy()
                b, se, t, p = R.test_carriers(idx, yr, yy) if len(idx) >= a.min_carriers else (np.nan,) * 4
                rows.append(dict(trait=tr.replace("_rint", ""), NAME=l.NAME, CHR=l.CHR, START=l.START, END=l.END, RULE=l.RULE, CNV=lab,
                                 n_carriers=len(idx), median_call_kb=round(med / 1e3, 1) if np.isfinite(med) else np.nan,
                                 BETA_SD=b, SE=se, T=t, P=p, N=len(P)))
    out = pd.DataFrame(rows)
    out.to_csv(a.out, sep="\t", index=False, float_format="%.6g")
    print(json.dumps(dict(loci=len(loci), traits=a.traits, model=f"OLS trait ~ carrier + {'+'.join(num + fac)}", min_carriers=a.min_carriers)))


if __name__ == "__main__":
    main()
