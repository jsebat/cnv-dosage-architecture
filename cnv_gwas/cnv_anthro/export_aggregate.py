"""Copy aggregate results to an export folder, enforcing a small-cell rule (All of Us: no count 1..20 may leave the
Workbench, nor anything that lets one be recomputed). Applied to every cohort so UKBB and AoU exports are identical in form.

Rules: any count column (NCNVs, n_carriers, cases/carriers fields, json counts) with 1 <= x <= min_cell is replaced by
'<=min_cell'; rows of per-probe / per-allele tables whose carrier count is in that range also lose BETA/SE/T/P
(an effect on <= 20 people is not exported). Person-level files (*genotypes*, pheno, .cnv) are never copied.
"""
import argparse
import json
import os
import shutil
import numpy as np
import pandas as pd

COUNT_COLS = ["NCNVs", "n_carriers", "n", "cases", "ctrl", "carriers_port", "carriers_omar", "n_probes"]
STAT_COLS = ["BETA", "SE", "T", "P", "BETA_SD", "raw_mean_diff_unadjusted", "FREQ"]
NEVER = ("genotype", "pheno", ".cnv", "allele_genotypes")


def scrub_json(o, k):
    if isinstance(o, dict):
        return {kk: scrub_json(v, k) for kk, v in o.items()}
    if isinstance(o, list):
        return [scrub_json(v, k) for v in o]
    if isinstance(o, (int, np.integer)) and not isinstance(o, bool) and 1 <= o <= k:
        return f"<={k}"
    return o


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True); ap.add_argument("--export", required=True)
    ap.add_argument("--min_cell", type=int, default=20)
    a = ap.parse_args()
    k = a.min_cell; n_sup = 0
    for root, _, files in os.walk(a.results):
        for f in files:
            if any(x in f for x in NEVER):
                continue
            src = os.path.join(root, f); rel = os.path.relpath(src, a.results); dst = os.path.join(a.export, rel)
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            if f.endswith((".tsv", ".tsv.gz")):
                d = pd.read_csv(src, sep="\t")
                cc = [c for c in COUNT_COLS if c in d.columns]
                small = np.zeros(len(d), bool)
                for c in cc:
                    v = pd.to_numeric(d[c], errors="coerce")
                    small |= v.between(1, k).to_numpy()
                for c in [c for c in STAT_COLS if c in d.columns]:
                    d[c] = d[c].astype(object); d.loc[small, c] = "suppressed"
                for c in cc:
                    v = pd.to_numeric(d[c], errors="coerce"); d[c] = d[c].astype(object); d.loc[v.between(1, k), c] = f"<={k}"
                n_sup += int(small.sum())
                d.to_csv(dst, sep="\t", index=False)
            elif f.endswith(".json"):
                json.dump(scrub_json(json.load(open(src)), k), open(dst, "w"), indent=1)
            elif f.endswith((".png", ".pdf", ".md", ".txt")):
                shutil.copy(src, dst)
    print(f"exported to {a.export}; rows with small cells suppressed: {n_sup}")


if __name__ == "__main__":
    main()
