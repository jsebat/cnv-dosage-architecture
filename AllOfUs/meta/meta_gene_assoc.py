"""Meta-analysis of the gene-based CNV association tables (gene_assoc_qt.py -> export_gene_assoc.py) across cohorts, at UCSD.

  python meta_gene_assoc.py --cohort UKBB=<export>/all/gene_assoc --cohort AoU_array=<export>/all/gene_assoc --out <new dir> [--traits height,bmi]

Per trait, gene and dosage: fixed-effect inverse-variance meta of beta (trait SD per carrier), SE, z, P; Cochran Q and its P, I2; cohorts
contributing (a gene a cohort did not export, e.g. under its small-cell rule, is left out for that cohort); carriers and N summed over the
contributing cohorts, carrier frequency = their ratio; per-cohort beta / SE / carriers kept; annotation columns (chr, start, end, s_het, HC,
TWAS z, segdup fraction) from the first cohort that has the gene. A one-cohort run returns that cohort's estimates.
Main meta: UK Biobank + All of Us arrays; All of Us WGS is reported separately (run it as its own --cohort set).
null model: none (analytic SEs, as gene_assoc_qt.py).
Output: <out>/<trait>/gene_assoc_meta.tsv.gz, meta_info.json.
"""
import argparse
import json
import os
import numpy as np
import pandas as pd
from scipy import stats

ANN = ["symbol", "chr", "start", "end", "shet", "HC", "twas_z", "segdup_frac"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", action="append", required=True, help="LABEL=<export>/<stratum>/gene_assoc (repeat)")
    ap.add_argument("--traits", default="height,bmi"); ap.add_argument("--out", required=True)
    a = ap.parse_args(); C = dict(x.split("=", 1) for x in a.cohort); labs = list(C)
    if os.path.exists(a.out):
        raise SystemExit(f"{a.out} exists; choose a new output folder")
    info = {}
    for t in a.traits.split(","):
        T = {}
        for l, d in C.items():
            f = os.path.join(d, t, "gene_assoc.tsv.gz")
            if os.path.exists(f):
                T[l] = pd.read_csv(f, sep="\t", dtype={"entrez": str}).drop_duplicates(["entrez", "dosage"]).set_index(["entrez", "dosage"])
        if not T:
            continue
        idx = T[list(T)[0]].index
        for l in list(T)[1:]:
            idx = idx.union(T[l].index)
        B = np.column_stack([T[l].beta.reindex(idx).to_numpy(float) for l in T]); S = np.column_stack([T[l].SE.reindex(idx).to_numpy(float) for l in T])
        ok = np.isfinite(B) & np.isfinite(S) & (S > 0); W = np.where(ok, 1 / np.where(ok, S, 1) ** 2, 0.0); sw = W.sum(1)
        with np.errstate(invalid="ignore", divide="ignore"):
            b = (W * np.where(ok, B, 0)).sum(1) / sw; se = np.sqrt(1 / sw); kc = ok.sum(1)
            Qs = (W * (np.where(ok, B, 0) - b[:, None]) ** 2).sum(1); I2 = np.where(Qs > 0, np.maximum(0, (Qs - (kc - 1)) / Qs), 0.0)
        if len(T) == 1:                 # one cohort: its own numbers, not a recomputation
            b, se = B[:, 0], S[:, 0]
        out = pd.DataFrame(index=idx).reset_index()
        for c in ANN:
            v = pd.Series(index=idx, dtype=object)
            for l in T:
                if c in T[l]:
                    v = v.where(v.notna(), T[l][c].reindex(idx))
            out[c] = v.to_numpy()
        nc = np.column_stack([T[l].n_carriers.reindex(idx).to_numpy(float) for l in T]); NN = np.column_stack([T[l].N.reindex(idx).to_numpy(float) for l in T])
        out["n_carriers"] = np.where(ok, nc, 0).sum(1).astype(int); out["N"] = np.where(ok, NN, 0).sum(1).astype(int); out["freq"] = out.n_carriers / out.N
        z = b / se
        out["beta"] = b; out["SE"] = se; out["z"] = z; out["P"] = 2 * stats.norm.sf(np.abs(z)) if len(T) > 1 else np.column_stack([T[l].P.reindex(idx) for l in T])[:, 0]
        out["k_cohorts"] = kc; out["cohorts"] = [";".join(l for l, o in zip(T, r) if o) for r in ok]
        out["Q"] = np.where(kc > 1, Qs, np.nan); out["P_Q"] = np.where(kc > 1, stats.chi2.sf(Qs, np.maximum(kc - 1, 1)), np.nan); out["I2"] = np.where(kc > 1, I2, np.nan)
        for j, l in enumerate(T):
            out[f"beta_{l}"] = B[:, j]; out[f"SE_{l}"] = S[:, j]; out[f"n_carriers_{l}"] = nc[:, j]
        out.insert(0, "trait", t); out = out.sort_values(["chr", "start", "dosage"])
        d = os.path.join(a.out, t); os.makedirs(d)
        out.to_csv(os.path.join(d, "gene_assoc_meta.tsv.gz"), sep="\t", index=False, float_format="%.6g")
        info[t] = dict(rows=int(len(out)), rows_2plus_cohorts=int((kc > 1).sum()), cohorts=list(T), rows_by_cohort={l: int(len(T[l])) for l in T})
        print(t, info[t])
    json.dump(dict(cohorts=labs, inputs=C, method="fixed-effect inverse variance per gene x dosage", traits=info), open(os.path.join(a.out, "meta_info.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
