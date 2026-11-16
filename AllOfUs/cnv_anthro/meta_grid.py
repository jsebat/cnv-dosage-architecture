"""Meta-analysis of the genome-wide grid / joint models (conga_grid_qt.py) across cohorts, with a matched empirical null.

Requirement: every cohort ran the same design with the SAME genome shifts (--perms_file from the reference cohort) and the
same LOO runs (--loo_mode drivers --major_loci <fixed list>, or --loo_mode all with the same loci files). Then replicate i
is the same annotation shift in every cohort and the meta statistic can be formed replicate by replicate:
    Z_meta = sum_c w_c t_c / sqrt(sum_c w_c^2),  w_c = sqrt(N_c)   (observed and each null replicate, each matched run)
Empirical z = (Z_obs - mean null) / sd null; two-sided empirical P; LOO held over the shared major-locus runs; no_loci and
minor_all as in the single-cohort summaries.

  python meta_grid.py --cohort UKBB=<export>/all/grid/<design>/<trait> --cohort AoU=... --out meta_grid_<design>_<trait>.tsv
"""
import argparse
import json
import os
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", action="append", required=True, help="LABEL=<grid dir holding one sub-folder per variant> (repeat)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    C = [tuple(x.split("=", 1)) for x in a.cohort]
    variants = sorted(set.intersection(*[set(os.listdir(d)) for _, d in C]))
    out = []
    for v in variants:
        obs, nul, N, perms, pfile = {}, {}, {}, set(), set()
        for lab, d in C:
            vd = os.path.join(d, v)
            if not os.path.exists(os.path.join(vd, "null_runs.npz")):
                print(f"WARNING: {vd} has no null_runs.npz; variant {v} skipped (rerun that cohort with the current conga_grid_qt.py)", flush=True)
                break
            info = json.load(open(os.path.join(vd, "run_info.json"))); N[lab] = float(info["N"]); perms.add((info.get("nrep"), info.get("perms_sha1")))
            pi = str(info.get("perms", ""))      # "loaded <n> shared shifts from <file>": replicate i = row i of that file in every cohort
            pfile.add(pi.split(" from ", 1)[1] if pi.startswith("loaded ") and " from " in pi else None)
            obs[lab] = pd.read_csv(os.path.join(vd, "observed_terms.tsv"), sep="\t").set_index(["run", "readout"]).t
            nul[lab] = np.load(os.path.join(vd, "null_runs.npz"))
        else:
            vname = obs[C[0][0]].index.get_level_values(0)  # noqa
            keys = [kk for kk in nul[C[0][0]].files]
            runs = sorted({kk.split("__", 1)[0] for kk in keys})
            runs = [r for r in runs if all(any(f.startswith(r + "__") for f in nul[lab].files) for lab, _ in C)]
            rd_names = [r for r in obs[C[0][0]].loc["full"].index]
            W = np.array([np.sqrt(N[lab]) for lab, _ in C])
            def meta_run(r):
                fk = [f for f in nul[C[0][0]].files if f.startswith(r + "__")][0]
                X = [np.asarray(nul[lab][fk], float) for lab, _ in C]
                nr = min(x.shape[0] for x in X)                                             # cohorts with different nrep: the common
                T = np.stack([x[:nr] for x in X])          # cohorts x reps x readouts      # first nr shifts (same file prefix) are matched
                o = np.array([[obs[lab].get((r, nm), np.nan) for nm in rd_names] for lab, _ in C])
                zo = np.nansum(W[:, None] * o, 0) / np.sqrt(np.sum(W ** 2))
                zn = np.nansum(W[:, None, None] * T, 0) / np.sqrt(np.sum(W ** 2))
                return zo, zn
            res = {r: meta_run(r) for r in runs}
            loo = [r for r in runs if r.startswith("loo_")]
            for j, nm in enumerate(rd_names):
                def emp(r):
                    zo, zn = res[r]; nu = zn[:, j]; nu = nu[np.isfinite(nu)]
                    if len(nu) < 10 or not np.isfinite(zo[j]):
                        return np.nan, np.nan, np.nan
                    mu, sd = nu.mean(), nu.std(ddof=1)
                    return zo[j], (zo[j] - mu) / sd, (1 + np.sum(np.abs(nu - mu) >= abs(zo[j] - mu))) / (len(nu) + 1)
                zo, ze, pe = emp("full"); s = np.sign(ze)
                lz = [emp(r)[1] for r in loo]
                out.append(dict(variant=v, readout=nm, cohorts=";".join(l for l, _ in C), Z_meta=zo, z_emp=ze, P_emp=pe,
                                loo_held=sum(1 for q in lz if np.isfinite(q) and abs(q) > 1.96 and np.sign(q) == s), loo_n=len(loo),
                                noloci_z_emp=emp("no_loci")[1] if "no_loci" in res else np.nan,
                                minor_all_z_emp=emp("minor_all")[1] if "minor_all" in res else np.nan,
                                **{f"t_{lab}": obs[lab].get(("full", nm), np.nan) for lab, _ in C},
                                shared_shifts=len(perms) == 1 or (len(pfile) == 1 and None not in pfile)))
    pd.DataFrame(out).to_csv(a.out, sep="\t", index=False, float_format="%.6g")
    print("wrote", a.out, len(out), "rows")


if __name__ == "__main__":
    main()
