"""Four-tile mirror signature (Jon 2026-11-22: "Is there some way to integrate differential signals of dup and del in these 4 tiles that would
give a stronger MIRROR signature?"). For one constraint stratum (LC or HC) and one CNV-size class of a conga_grid_qt.py run, the four terms
  t_DEL(up), t_DUP(up), t_DEL(down), t_DUP(down)       (up = TWAS z > 3, down = z < -3; t = analytic t of each joint-model term)
carry a mirror when deletions of 'up' genes and duplications of 'down' genes lower the trait and the other two raise it, i.e. along the
pre-specified direction s = (-1, +1, +1, -1). Two statistics, both fixed before looking at results:
  M4   = s . t / 2                                (equal weights; = sum of the four tiles' z in the mirror direction / 2)
  Mopt = w . t / sqrt(w' S w), w = S^-1 s         (most powerful linear test of direction s; S = covariance of the four t across the null
                                                   replicates of the full run, so the weights never see the observed data)
For comparison, the current readouts: single-class mirrors (t_DUP - t_DEL)/sqrt2 and, from variant size4s, the signed-count mirror.
null model: the same genome shifts of the gene annotation (null_runs.npz), statistic recomputed in every replicate of every run; z = (obs - mean)
/ sd, two-sided empirical P; LOO = each major-locus run vs its own matched null; no_loci = the 34 recurrent intervals removed.
  python mirror4_summary.py --grid <class_models root>/grid --out mirror4_summary.tsv
"""
import argparse, os
import numpy as np, pandas as pd

S = np.array([-1.0, 1.0, 1.0, -1.0])


def readout_names(var):
    n = [f"{x}|DEL" for x in var["DEL"]] + [f"{x}|DUP" for x in var["DUP"]]
    return n + [f"{x}|mirror" for x in var["DEL"] if x in var["DUP"] and x != "n_all"]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--grid", required=True); ap.add_argument("--design", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    import json
    design = {v["name"]: v for v in json.load(open(a.design))}
    rows = []
    for trait in sorted(os.listdir(a.grid)):
        for variant, prefixes in (("size4", ["N1_1~", "N2_5~", "N6_20~", "N21_inf~"]), ("size2", ["N1_1~", "N2_inf~"]), ("size2j", ["N1_1~", "N2_inf~"]), ("size3c", ["N1_1~", "N2_inf~"]), ("size2t", ["N1_1~", "N2_inf~"]), ("all2t", ["N1_inf~"]), ("size3a", ["N1_1~", "N2_inf~"]), ("joint6", [""])):
            d = os.path.join(a.grid, trait, variant)
            if not os.path.exists(os.path.join(d, "null_runs.npz")):
                continue
            names = readout_names(design[variant]); nz = np.load(os.path.join(d, "null_runs.npz"))
            obs = pd.read_csv(os.path.join(d, "observed_terms.tsv"), sep="\t").set_index(["run", "readout"]).t
            runs = sorted({k.split("__")[0] for k in nz.files})
            for pre in prefixes:
                for con in ("LC", "HC", "all"):
                    cols = [f"{pre}cls:up:{con}|DEL", f"{pre}cls:up:{con}|DUP", f"{pre}cls:down:{con}|DEL", f"{pre}cls:down:{con}|DUP"]
                    if not all(c in names for c in cols):
                        continue
                    j = [names.index(c) for c in cols]
                    Sig = np.cov(nz[f"full__{variant}"][:, j].T); w = np.linalg.solve(Sig, S)
                    res = {}
                    for rn in runs:
                        T0 = np.array([obs[(rn, c)] for c in cols]); TN = nz[f"{rn}__{variant}"][:, j]
                        # armDEL / armDUP (Jon 2026-11-22, Shanta Fig. 3F style dosage grid): signed TWAS term of each dosage,
                        # (t_up - t_down)/sqrt2, raw sign (positive: genes whose higher expression raises the trait push it up); M4 = (armDUP - armDEL)/sqrt2
                        for nm, f in (("M4", lambda x: x @ S / 2), ("Mopt", lambda x: x @ w / np.sqrt(w @ Sig @ w)),
                                      ("armDEL", lambda x: x @ np.array([1.0, 0, -1.0, 0]) / np.sqrt(2)), ("armDUP", lambda x: x @ np.array([0, 1.0, 0, -1.0]) / np.sqrt(2))):
                            o = f(T0); nu = f(TN); nu = nu[np.isfinite(nu)]; mu, sd = nu.mean(), nu.std(ddof=1)
                            res[(nm, rn)] = ((o - mu) / sd, (1 + np.sum(np.abs(nu - mu) >= abs(o - mu))) / (len(nu) + 1), o)
                    loo = [r for r in runs if r.startswith("loo_")] + ["minor_all"]
                    for nm in ("M4", "Mopt", "armDEL", "armDUP"):
                        z, p, o = res[(nm, "full")]; held = sum(1 for r in loo if (nm, r) in res and res[(nm, r)][1] < .05 and np.sign(res[(nm, r)][0]) == np.sign(z))
                        rows.append(dict(trait=trait, variant=variant, size=pre.rstrip("~") or "all sizes", stratum=con, stat=nm, obs=o, z=z, P=p,
                                         loo_held=held, loo_n=sum(1 for r in loo if (nm, r) in res),
                                         noloci_z=res.get((nm, "no_loci"), (np.nan,))[0], noloci_P=res.get((nm, "no_loci"), (np.nan, np.nan))[1],
                                         weights=" ".join(f"{x:+.2f}" for x in (w / np.abs(w).sum() * 4)) if nm == "Mopt" else {"M4": "-1 +1 +1 -1", "armDEL": "+1 0 -1 0", "armDUP": "0 +1 0 -1"}[nm]))
    R = pd.DataFrame(rows); R.to_csv(a.out, sep="\t", index=False, float_format="%.4g")
    pd.set_option("display.width", 250); print(R.to_string(index=False, float_format=lambda x: f"{x:.3g}"))


if __name__ == "__main__":
    main()
