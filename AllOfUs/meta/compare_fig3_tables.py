"""Compare two sets of Fig. 3 tables (e.g. the UK Biobank run folders vs a UKBB-only meta_fig3.py output, or two cohorts' code paths).

  python compare_fig3_tables.py --a <meta dir | size2t_mirror=<dir>,all2t_mirror=<dir>,totdose=<dir>,gene_mirror=<dir>> --b <same form> [--tol 1e-4]

Compares every numeric column of grid_summary.tsv, mirror4_cross_summary.tsv, dose_summary.tsv and gene_mirror_summary.tsv row by row (rows
matched on their key columns), and the calls the figure makes from them (P < 0.05; held in every LOO run). Exit 1 if a table is missing on
one side, a value differs by more than --tol (relative, or absolute for |x| < 1) or a figure call differs.
null model: none (a comparison of finished tables).
"""
import argparse
import glob
import os
import sys
import numpy as np
import pandas as pd

SUB = ["size2t_mirror", "all2t_mirror", "totdose", "gene_mirror", "dose_bins", "uni_N1_inf", "uni_N1_1", "uni_N2_inf"]
KEYS = {"grid_summary.tsv": ["variant", "readout"], "mirror4_cross_summary.tsv": ["trait", "variant", "size", "stratum", "stat"],
        "dose_summary.tsv": ["trait", "variant", "size", "stratum", "stat"], "gene_mirror_summary.tsv": ["size", "stratum", "stat"],
        "fig3_dose_bins.tsv": ["trait", "run", "partition", "size", "cls", "side", "bin"], "fig3_dose_shape.tsv": ["trait", "run", "partition", "size", "cls", "term"]}
SKIP = {"cohorts", "weights", "cohort"}


def roots(x):
    return dict(p.split("=", 1) for p in x.split(",")) if "=" in x else {s: os.path.join(x, s) for s in SUB}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--a", required=True); ap.add_argument("--b", required=True); ap.add_argument("--tol", type=float, default=1e-4)
    a = ap.parse_args(); A, B = roots(a.a), roots(a.b); bad = 0; rows = []
    for s in SUB:
        if s not in A or s not in B:
            continue
        fa = sorted(os.path.relpath(f, A[s]) for k in KEYS for f in glob.glob(os.path.join(A[s], "**", k), recursive=True) if "code_snapshot" not in f)
        for rel in fa:
            pb = os.path.join(B[s], rel)
            if not os.path.exists(pb):
                print("MISSING in b:", s, rel); bad += 1; continue
            if os.path.getsize(os.path.join(A[s], rel)) < 5:
                continue
            x = pd.read_csv(os.path.join(A[s], rel), sep="\t"); y = pd.read_csv(pb, sep="\t"); key = KEYS[os.path.basename(rel)]
            m = x.merge(y, on=key, how="outer", suffixes=("_a", "_b"), indicator=True)
            extra = int((m._merge == "right_only").sum()); m = m[m._merge != "right_only"]   # rows only in b (a newer summary script): reported, not an error
            nomatch = int((m._merge != "both").sum()); worst = 0.0; wcol = ""; calls = 0
            # total-genes terms do not depend on the gene annotation, so their shift null is degenerate (sd ~ 0) and its z is numerical noise:
            # only the analytic t is compared for them (the figure reads nothing else)
            nall = m.readout.str.rsplit("|", n=1).str[0].str.endswith("n_all").to_numpy() if "readout" in m.columns else np.zeros(len(m), bool)
            for c in [c for c in x.columns if c not in key and c in y.columns and c not in SKIP]:
                if nall.all() and c != "t_analytic":
                    continue
                u, v = m[f"{c}_a"], m[f"{c}_b"]
                if pd.api.types.is_numeric_dtype(u) and pd.api.types.is_numeric_dtype(v):
                    u, v = u.to_numpy(float), v.to_numpy(float)
                    if c != "t_analytic":
                        u, v = u[~nall], v[~nall]
                    both = np.isfinite(u) & np.isfinite(v)
                    nomatch += int((np.isfinite(u) != np.isfinite(v)).sum())
                    d = np.abs(u[both] - v[both]) / np.maximum(1.0, np.abs(u[both]))
                    if len(d) and d.max() > worst:
                        worst, wcol = float(d.max()), c
                    if c in ("P", "P_emp", "minor_all_P_emp", "noloci_P", "noloci_P_emp"):
                        calls += int(((u[both] < .05) != (v[both] < .05)).sum())
                    if c in ("loo_held", "loo_n", "n_loo", "n_genes"):
                        calls += int((u[both] != v[both]).sum())
                else:
                    nomatch += int((u.fillna("").astype(str)[~nall] != v.fillna("").astype(str)[~nall]).sum())
            ok = nomatch == 0 and calls == 0 and worst <= a.tol; bad += (not ok)
            rows.append(dict(set=s, table=rel.replace(os.sep, "/"), rows=len(x), extra_rows_b=extra, unmatched=nomatch, max_rel_diff=worst, column=wcol, call_changes=calls, ok=ok))
    R = pd.DataFrame(rows); pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 80)
    print(R.to_string(index=False)); print(f"\n{len(R)} tables compared; {int((~R.ok).sum()) if len(R) else 0} differ beyond tol {a.tol}; missing {bad - int((~R.ok).sum()) if len(R) else bad}")
    sys.exit(1 if bad or not len(R) else 0)


if __name__ == "__main__":
    main()
