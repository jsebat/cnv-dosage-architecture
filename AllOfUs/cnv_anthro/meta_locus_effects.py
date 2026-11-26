"""Fixed-effect inverse-variance meta-analysis of per-locus CNV effects (locus_effects.py exports) across cohorts (dosage paper Fig. 2A,
Jon 2026-11-22: "Use a meta-analysis of just UK biobank and AoU for the effect size estimates"). Per trait x locus x dosage: pooled beta
(SD units), SE, P, Cochran Q P, number of cohorts contributing (cells suppressed or with no carriers in a cohort are left out for that
cohort and flagged), carriers summed. null model: none (analytic SEs).
  python meta_locus_effects.py --cohort UKBB=<dir with locus_effects_<trait>.tsv> --cohort AoU=<dir> --traits height,bmi --out <dir>
"""
import argparse
import os
import numpy as np
import pandas as pd
from scipy import stats


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--cohort", action="append", required=True)
    ap.add_argument("--traits", default="height,bmi"); ap.add_argument("--out", required=True)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    C = [tuple(x.split("=", 1)) for x in a.cohort]
    for t in a.traits.split(","):
        X = []
        for lab, d in C:
            x = pd.read_csv(os.path.join(d, f"locus_effects_{t}.tsv"), sep="\t", dtype={"n_carriers": str}); x["cohort"] = lab; X.append(x)
        X = pd.concat(X)
        for c in ("BETA_SD", "SE"):
            X[c] = pd.to_numeric(X[c], errors="coerce")
        X["usable"] = np.isfinite(X.BETA_SD) & np.isfinite(X.SE) & (X.SE > 0)
        rows = []
        for (nm, cnv), g in X.groupby(["NAME", "CNV"]):
            u = g[g.usable]; r = dict(trait=t, NAME=nm, CNV=cnv)
            for _, z in g.iterrows():
                r[f"{z.cohort}_beta"] = z.BETA_SD; r[f"{z.cohort}_SE"] = z.SE; r[f"{z.cohort}_n"] = z.n_carriers
            if len(u):
                w = 1 / u.SE ** 2; b = (w * u.BETA_SD).sum() / w.sum(); se = w.sum() ** -0.5; Q = (w * (u.BETA_SD - b) ** 2).sum()
                r.update(beta=b, SE=se, P=2 * stats.norm.sf(abs(b / se)), Q_P=stats.chi2.sf(Q, len(u) - 1) if len(u) > 1 else np.nan,
                         k_cohorts=len(u), cohorts_used="+".join(u.cohort), n_carriers=int(pd.to_numeric(u.n_carriers, errors="coerce").sum()))
            else:
                r.update(beta=np.nan, SE=np.nan, P=np.nan, Q_P=np.nan, k_cohorts=0, cohorts_used="", n_carriers=0)
            rows.append(r)
        pd.DataFrame(rows).to_csv(os.path.join(a.out, f"meta_locus_effects_{t}.tsv"), sep="\t", index=False, float_format="%.6g")
        print(t, len(rows), "locus x dosage cells")


if __name__ == "__main__":
    main()
