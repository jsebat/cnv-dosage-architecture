"""Check that the all-CNV rows of fig3_dose_bins.py reproduce the dosage paper's Fig. 2a-h tables (class_dose_bins.py, partitions 'all' and
'twas_x_con'): same bins, carriers, mean genes, beta and SE, in the full and no_loci runs. Exit 1 on any difference beyond --tol.
null model: none (comparison of analytic OLS tables).
  python check_bins_vs_fig2.py --new <fig3 dose_bins dir with <trait>/fig3_dose_bins.tsv> --old <class_bins dir with <trait>/class_dose_bins.tsv>
"""
import argparse, os, sys
import numpy as np, pandas as pd

ap = argparse.ArgumentParser(); ap.add_argument("--new", required=True); ap.add_argument("--old", required=True); ap.add_argument("--tol", type=float, default=1e-5)
a = ap.parse_args(); bad = 0; key = ["run", "partition", "cls", "side", "bin"]
for t in ("height", "bmi"):
    x = pd.read_csv(os.path.join(a.new, t, "fig3_dose_bins.tsv"), sep="\t"); y = pd.read_csv(os.path.join(a.old, t, "class_dose_bins.tsv"), sep="\t")
    x = x[(x["size"] == "N1_inf") & x.partition.isin(["all", "twas_x_con"])]; y = y[y.partition.isin(["all", "twas_x_con"])]
    m = x.merge(y, on=key, how="outer", suffixes=("_new", "_old"), indicator=True)
    print(f"{t}: rows new {len(x)}, old {len(y)}, matched {(m._merge == 'both').sum()}, unmatched {(m._merge != 'both').sum()}"); bad += int((m._merge != "both").sum())
    for c in ("beta", "SE", "mean_genes", "n"):
        u = pd.to_numeric(m[f"{c}_new"], errors="coerce").to_numpy(float); v = pd.to_numeric(m[f"{c}_old"], errors="coerce").to_numpy(float)
        both = np.isfinite(u) & np.isfinite(v); d = np.abs(u[both] - v[both]) / np.maximum(1, np.abs(v[both]))
        only = int((np.isfinite(u) != np.isfinite(v)).sum())      # expected: n of the no_loci run (not written by the new script), bins with exactly 20 people
        print(f"   {c:10s} compared {int(both.sum()):3d}  max rel diff {d.max() if len(d) else 0:.2e}  present in one only {only}")
        bad += int((d > a.tol).sum())
        if c in ("beta", "SE") and only:
            w = m[np.isfinite(u) != np.isfinite(v)]; print(w[key + [f"n_new", f"n_old"]].to_string(index=False))
print("OK: the all-CNV bins equal the Fig. 2 tables" if not bad else f"DIFFERENCES: {bad}")
sys.exit(1 if bad else 0)
