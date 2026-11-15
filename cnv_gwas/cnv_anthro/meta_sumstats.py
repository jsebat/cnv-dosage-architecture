"""Inverse-variance meta-analysis of exported CNV summary statistics across cohorts (e.g. UKBB array + AoU srWGS).

Tables (matched by key; rows with suppressed statistics are skipped for that cohort):
  probes   taskA_gwas/<trait>/sumstats_{del,dup}.tsv.gz   key CHR, BP (+ TYPE from the file name). AoU tests the UKBB probe positions too
           (gwas.ref_probes), so the UKBB grid of probes is the meta grid; use --ref_label to restrict to the reference cohort's probes.
  alleles  taskA_alleles/<trait>/allele_assoc.tsv          key allele
  loci     locus_effects/locus_effects_<trait>.tsv         key NAME, CNV
Output per table: fixed-effect beta/SE/P, DerSimonian-Laird random-effects beta/SE/P, Cochran Q, I2, k cohorts, total carriers where available.

  python meta_sumstats.py --cohort UKBB=<export>/all --cohort AoU=<export>/all --trait height --out meta/height_all [--ref_label UKBB]
"""
import argparse
import os
import numpy as np
import pandas as pd
from scipy import stats


def ivw(b, s):
    b = np.asarray(b, float); s = np.asarray(s, float); ok = np.isfinite(b) & np.isfinite(s) & (s > 0)
    b, s = b[ok], s[ok]; k = len(b)
    if k == 0:
        return dict(k=0)
    w = 1 / s ** 2; bf = np.sum(w * b) / w.sum(); sf = np.sqrt(1 / w.sum())
    Q = float(np.sum(w * (b - bf) ** 2)) if k > 1 else 0.0
    tau2 = max(0.0, (Q - (k - 1)) / (w.sum() - np.sum(w ** 2) / w.sum())) if k > 1 else 0.0
    wr = 1 / (s ** 2 + tau2); br = np.sum(wr * b) / wr.sum(); sr = np.sqrt(1 / wr.sum())
    return dict(k=k, BETA_FE=bf, SE_FE=sf, P_FE=2 * stats.norm.sf(abs(bf / sf)), BETA_RE=br, SE_RE=sr, P_RE=2 * stats.norm.sf(abs(br / sr)),
                Q=Q, P_Q=stats.chi2.sf(Q, k - 1) if k > 1 else np.nan, I2=max(0.0, (Q - (k - 1)) / Q) if Q > 0 else 0.0, tau2=tau2)


def load(cohorts, rel, key, bcol, scol, ncol):
    rows = []
    for lab, root in cohorts:
        p = os.path.join(root, rel)
        if not os.path.exists(p):
            continue
        d = pd.read_csv(p, sep="\t", dtype=str)
        d["cohort"] = lab
        d["b"] = pd.to_numeric(d[bcol], errors="coerce"); d["s"] = pd.to_numeric(d[scol], errors="coerce")
        d["n"] = pd.to_numeric(d[ncol], errors="coerce") if ncol in d.columns else np.nan
        rows.append(d[key + ["cohort", "b", "s", "n"]])
    return pd.concat(rows) if rows else None


def meta_table(D, key):
    out = []
    for kv, g in D.groupby(key, sort=False):
        r = ivw(g.b, g.s)
        if r["k"] == 0:
            continue
        kv = kv if isinstance(kv, tuple) else (kv,)
        out.append({**dict(zip(key, kv)), **r, "cohorts": ";".join(g.cohort[np.isfinite(g.b) & np.isfinite(g.s)]),
                    "carriers_total": g.n.sum(min_count=1)})
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", action="append", required=True, help="LABEL=export stratum dir (repeat)")
    ap.add_argument("--trait", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--ref_label", default=None, help="keep only probe positions present in this cohort")
    a = ap.parse_args()
    C = [tuple(x.split("=", 1)) for x in a.cohort]; os.makedirs(a.out, exist_ok=True); t = a.trait
    for typ in ["del", "dup"]:
        D = load(C, f"taskA_gwas/{t}/sumstats_{typ}.tsv.gz", ["CHR", "BP"], "BETA", "SE", "NCNVs")
        if D is None:
            continue
        if a.ref_label:
            ref = set(map(tuple, D[D.cohort == a.ref_label][["CHR", "BP"]].values))
            D = D[[tuple(x) in ref for x in D[["CHR", "BP"]].values]]
        M = meta_table(D, ["CHR", "BP"]); M.insert(2, "TYPE", typ.upper())
        M.to_csv(os.path.join(a.out, f"meta_probes_{typ}.tsv.gz"), sep="\t", index=False, float_format="%.6g")
    D = load(C, f"taskA_alleles/{t}/allele_assoc.tsv", ["allele"], "BETA_SD", "SE", "n_carriers")
    if D is not None:
        meta_table(D, ["allele"]).to_csv(os.path.join(a.out, "meta_alleles.tsv"), sep="\t", index=False, float_format="%.6g")
    D = load(C, f"locus_effects/locus_effects_{t}.tsv", ["NAME", "CNV"], "BETA_SD", "SE", "n_carriers")
    if D is not None:
        meta_table(D, ["NAME", "CNV"]).to_csv(os.path.join(a.out, "meta_loci.tsv"), sep="\t", index=False, float_format="%.6g")
    print("wrote", a.out)


if __name__ == "__main__":
    main()
