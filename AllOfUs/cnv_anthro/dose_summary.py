"""Dose-aligned contrasts of the size2t joint model (Jon 2026-11-23: 'add a "Dose aligned" term that is also derived from the 2x2: where dup/del
up/down all go in the same direction ... make an alternate version of the figure with that joint model'). Contrasts of the SAME fit (a term that is
a linear combination of existing terms cannot be added to the model; it is collinear), each per CNV-trait x TWAS-trait pair and CNV size:
  D4      = (t_DEL up + t_DUP up + t_DEL down + t_DUP down) / 2, per LC / HC   (orange / brown: the 2x2 tiles all pushing the trait the same way;
                                                                               companion of the mirror M4 = (-1, +1, +1, -1) / 2)
  doseOth = (t_DEL + t_DUP) / sqrt2 of 'weak or no TWAS' HC                     (green: the constraint U-shape beyond an ordinary gene)
  doseTot = (t_DEL + t_DUP) / sqrt2 of the total-genes term                   (yellow: per-gene U-shape of an ordinary gene; Jon: 'What do you
                                                                               think is best?' -> total genes, same joint model as every arrow)
Positive = both dosages raise the trait.
null model: D4 and doseOth: the genome shifts of the gene annotation of the fit (null_runs.npz, 2,000 shared shifts), statistic recomputed in every
replicate of every run; z = (obs - mean) / sd, two-sided empirical P; LOO = each major-locus run vs its own matched null, plus minor_all.
doseTot: total genes hit does not depend on the annotation, so no shift null exists: analytic, z = (t_DEL + t_DUP)/sqrt2 (DEL and DUP terms of
different carriers, taken as independent), two-sided normal P; LOO = analytic P < 0.05 with the same sign in each run.
  python dose_summary.py --grid <run>/grid --variant size2t --out <run>/dose_summary.tsv
"""
import argparse, os
import numpy as np, pandas as pd
from scipy.stats import norm


def readout_names(var):
    n = [f"{x}|DEL" for x in var["DEL"]] + [f"{x}|DUP" for x in var["DUP"]]
    return n + [f"{x}|mirror" for x in var["DEL"] if x in var["DUP"] and x != "n_all"]


def nz_cols(d):
    z = np.load(os.path.join(d, "null_runs.npz")); return z[f"full__{os.path.basename(d)}"].shape[1]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--grid", required=True); ap.add_argument("--variant", default="size2t"); ap.add_argument("--out", required=True)
    ap.add_argument("--design", required=True)
    a = ap.parse_args(); rows = []
    import json
    design = {v["name"]: v for v in json.load(open(a.design))}
    for pair in sorted(os.listdir(a.grid)):
        d = os.path.join(a.grid, pair, a.variant)
        if not os.path.exists(os.path.join(d, "null_runs.npz")):
            continue
        ob = pd.read_csv(os.path.join(d, "observed_terms.tsv"), sep="\t")
        names = readout_names(design[a.variant])                 # column order of the null matrices (as conga_grid_qt.py / mirror4_summary.py)
        assert nz_cols(d) == len(names), (pair, nz_cols(d), len(names))
        obs = ob.set_index(["run", "readout"]).t
        nz = np.load(os.path.join(d, "null_runs.npz"))
        runs = sorted({k.split("__")[0] for k in nz.files})
        loo = [r for r in runs if r.startswith("loo_")] + ["minor_all"]
        for size in ("N1_1", "N2_inf", "N1_inf"):   # N1_inf: totdose run (total genes hit, all CNV sizes; Jon 2026-11-23 yellow total)
            specs = [(con, "D4", [f"{size}~cls:{c}:{con}|{x}" for c in ("up", "down") for x in ("DEL", "DUP")], np.array([1., 1, 1, 1]) / 2, True)
                     for con in ("LC", "HC")]
            # consistency check: the mirror recomputed here must equal mirror4_summary.py's M4 (same columns, weights (-1, +1, +1, -1)/2)
            specs += [(con, "M4check", [f"{size}~cls:{c}:{con}|{x}" for c in ("up", "down") for x in ("DEL", "DUP")], np.array([-1., 1, 1, -1]) / 2, True)
                      for con in ("LC", "HC")]
            specs += [("HC", "doseOth", [f"{size}~cls:other:HC|DEL", f"{size}~cls:other:HC|DUP"], np.ones(2) / np.sqrt(2), True),
                      ("total", "doseTot", [f"{size}~n_all|DEL", f"{size}~n_all|DUP"], np.ones(2) / np.sqrt(2), False),
                      ("total", "mirTot", [f"{size}~n_all|DEL", f"{size}~n_all|DUP"], np.array([-1., 1]) / np.sqrt(2), False)]
            for con, stat, cols, w, shift in specs:
                if not all((("full", c) in obs.index) for c in cols):
                    continue
                res = {}
                for rn in ["full"] + loo + ["no_loci"]:
                    if not all(((rn, c) in obs.index) for c in cols):
                        continue
                    o = float(np.array([obs[(rn, c)] for c in cols]) @ w)
                    if shift:
                        key = f"{rn}__{a.variant}"
                        if key not in nz.files:
                            continue
                        nm = nz[key]; j = [names.index(c) for c in cols]; nu = nm[:, j] @ w; nu = nu[np.isfinite(nu)]
                        mu, sd = nu.mean(), nu.std(ddof=1)
                        res[rn] = ((o - mu) / sd, (1 + np.sum(np.abs(nu - mu) >= abs(o - mu))) / (len(nu) + 1), o)
                    else:
                        res[rn] = (o, 2 * norm.sf(abs(o)), o)
                if "full" not in res:
                    continue
                z, p, o = res["full"]; lr = [r for r in loo if r in res and (not shift or f"{r}__{a.variant}" in nz.files)]
                held = sum(1 for r in lr if res[r][1] < .05 and np.sign(res[r][0]) == np.sign(z))
                rows.append(dict(trait=pair, variant=a.variant, size=size, stratum=con, stat=stat, obs=o, z=z, P=p, loo_held=held, loo_n=len(lr),
                                 noloci_z=res.get("no_loci", (np.nan,) * 3)[0], noloci_P=res.get("no_loci", (np.nan,) * 3)[1],
                                 test="genome-shift null" if shift else "analytic"))
    out = pd.DataFrame(rows); out.to_csv(a.out, sep="\t", index=False); print(out.to_string(index=False))


if __name__ == "__main__":
    main()
