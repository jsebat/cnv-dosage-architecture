# Meta-analysis of the person-level ConGA context runs (conga_ctx_qt.py, run_pipeline.py step 'ctx') across cohorts, replicate by
# replicate under the shared shifts, written in the engine's own folder layout so conga_contrasts.py (cnv_gwas/meta_aou_2026-11-24)
# draws the mirror / dose-aligned readouts, FDR and locus-robust marks from UK Biobank + All of Us unchanged.
# null model: genome shift of the gene annotation; every cohort used the same shared shift file, so replicate i is the same shift
#   everywhere and the meta null is formed replicate by replicate (matched null); removal runs matched in the same way
"""As meta_grid.py does for the grid models:
    Z_meta = sum_c w_c z_c / sqrt(sum_c w_c^2),  w_c = sqrt(N_c)
for the observed score z of every set and for every null replicate of every run (full and each locus-removal run). Cohorts must share
the shift file, the 10-Mb minimum and the gene annotation (info.json: shifts_sha1, KMIN_10Mb, M); with different nrep the common first
replicates are used. A set is combined only where every cohort tested it. z_emp, P_emp and FDR are then computed exactly as in
conga_ctx_qt.py, so the output folder can be passed to conga_contrasts.py as --ctx_run and --loo_run.

Input: one folder per cohort holding <mode>_<trait>_<dosage>/ sub-folders (the 'ctx' folder of an export:
<export>/all/ctx for All of Us; the UKBB runs conga_ctx_*/ and conga_supp_*/, or their export). Output: the same sub-folders, meta.

  python meta_ctx.py --cohort UKBB=<dir> --cohort AoU=<dir> --out <new dir>
  python ../../cnv_gwas/meta_aou_2026-11-24/conga_contrasts.py --stage B --ctx_run <out> --loo_run <out> --groups <signal groups dir> \
      --out <figure dir> --cohort_label "UK Biobank + All of Us"
"""
import argparse
import json
import os
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control as fdc


def load1(d, mode):
    D = pd.read_csv(os.path.join(d, f"conga_ctx_{mode}.tsv.gz"), sep="\t")
    nul = {r: np.load(os.path.join(d, f"null_z_{r}.npy")) for r in D.run.unique() if os.path.exists(os.path.join(d, f"null_z_{r}.npy"))}
    return D, nul, json.load(open(os.path.join(d, "info.json")))


def load(dirs, sub, mode):
    """one cohort may come in several folders (UK Biobank: context run + signal-group run); per run, rows and null columns are
    concatenated in the same order. The folders must share N, shifts and annotation."""
    parts = [load1(os.path.join(d, sub), mode) for d in dirs if os.path.exists(os.path.join(d, sub, "info.json"))]
    parts = [p for p in parts if "note" not in p[2]]
    if not parts:
        return None
    if len({(p[2]["N"], p[2]["shifts_sha1"], p[2]["KMIN_10Mb"], p[2]["M"]) for p in parts}) > 1:
        raise SystemExit(f"{dirs} {sub}: folders of one cohort differ in N, shifts or annotation")
    runs = [r for r in parts[0][1] if all(r in p[1] for p in parts)]
    D = pd.concat([p[0][p[0].run == r] for r in runs for p in parts], ignore_index=True)
    nr = {r: min(p[1][r].shape[0] for p in parts) for r in runs}
    nul = {r: np.hstack([p[1][r][:nr[r]] for p in parts]) for r in runs}
    for r in runs:      # duplicated set ids across folders: keep the first
        m = D.run == r; dup = D[m].set_id.duplicated().to_numpy()
        if dup.any():
            nul[r] = nul[r][:, ~dup]; D = D.drop(D[m].index[dup])
    return D.reset_index(drop=True), nul, parts[0][2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", action="append", required=True, help="LABEL=<folder with <mode>_<trait>_<dosage>/ runs> (repeat)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    C = [(lab, d.split(",")) for lab, d in (x.split("=", 1) for x in a.cohort)]
    if os.path.exists(a.out):
        raise SystemExit(f"{a.out} exists; choose a new output folder")
    subs = sorted(set.intersection(*[{s for d in ds if os.path.isdir(d) for s in os.listdir(d)
                                       if os.path.exists(os.path.join(d, s, "info.json")) and s.split("_", 1)[0] in ("full", "loo")}
                                      for _, ds in C]))
    log = []
    for sub in subs:
        mode = sub.split("_", 1)[0]
        parts = {lab: load(ds, sub, mode) for lab, ds in C}
        if any(p is None for p in parts.values()):
            print("skip (no sets in a cohort):", sub); continue
        keyset = {(p[2]["shifts_sha1"], p[2]["KMIN_10Mb"], p[2]["M"]) for p in parts.values()}
        shared = len(keyset) == 1
        if not shared:
            print(f"WARNING {sub}: cohorts did not use the same shifts / annotation {keyset}; z_emp is NOT replicate-matched", flush=True)
        N = {lab: float(p[2]["N"]) for lab, p in parts.items()}; W = np.array([np.sqrt(N[lab]) for lab, _ in C])
        od = os.path.join(a.out, sub); os.makedirs(od)
        rows = []
        runs = [r for r in parts[C[0][0]][0].run.unique() if all(r in p[1] for p in parts.values())]
        for r in runs:
            tabs = {lab: p[0][p[0].run == r].reset_index(drop=True) for lab, p in parts.items()}
            sets = sorted(set.intersection(*[set(t.set_id) for t in tabs.values()]))
            obs = np.stack([tabs[lab].set_index("set_id").z.reindex(sets).to_numpy(float) for lab, _ in C])
            col = {lab: pd.Series(np.arange(len(tabs[lab])), index=tabs[lab].set_id).reindex(sets).to_numpy() for lab, _ in C}
            nr = min(parts[lab][1][r].shape[0] for lab, _ in C)
            T = np.stack([parts[lab][1][r][:nr][:, col[lab]].astype(float) for lab, _ in C])      # cohorts x reps x sets
            ok_o = np.isfinite(obs); zo = np.nansum(W[:, None] * obs, 0) / np.sqrt(np.sum((W[:, None] * ok_o) ** 2, 0))
            ok_n = np.isfinite(T); zn = np.nansum(W[:, None, None] * T, 0) / np.sqrt(np.sum((W[:, None, None] * ok_n) ** 2, 0))
            zo[~ok_o.all(0)] = np.nan                                     # a set is combined only where every cohort has its z
            zn = zn.astype(np.float32); np.save(os.path.join(od, f"null_z_{r}.npy"), zn)
            ge = (np.abs(zn) >= np.abs(zo)[None, :] - 1e-12).sum(0); ok = np.isfinite(zn).sum(0)
            mu = np.nanmean(zn, 0); sd = np.nanstd(zn, 0)
            t0 = tabs[C[0][0]].set_index("set_id").reindex(sets)
            d = pd.DataFrame(dict(trait=t0.trait.to_numpy(), dosage=t0.dosage.to_numpy(), run=r, set_id=sets, n_genes=t0.n_genes.to_numpy(),
                                  z=zo, z_emp=(zo - mu) / np.where(sd > 0, sd, np.nan), P_emp=(1 + ge) / (ok + 1), n_rep=ok,
                                  null_mean=mu, null_sd=sd, n_ge=ge, family=t0.family.to_numpy() if "family" in t0 else "all",
                                  shared_shifts=shared, **{f"z_{lab}": obs[i] for i, (lab, _) in enumerate(C)}))
            rows.append(d)
        D = pd.concat(rows, ignore_index=True); D["FDR_BH"] = np.nan
        for _, idx in D.groupby(["run", "family"]).groups.items():
            p = D.loc[idx, "P_emp"].to_numpy(); m = np.isfinite(p)
            if m.any():
                D.loc[np.asarray(idx)[m], "FDR_BH"] = fdc(p[m], method="bh")
        D.to_csv(os.path.join(od, f"conga_ctx_{mode}.tsv.gz"), sep="\t", index=False, float_format="%.5g")
        i0 = parts[C[0][0]][2]
        json.dump(dict(trait=i0["trait"], dosage=i0["dosage"], N=sum(N.values()), N_by_cohort=N, weights="sqrt(N)", M=i0["M"],
                       KMIN_10Mb=i0["KMIN_10Mb"], shifts_sha1=i0["shifts_sha1"], shared_shifts=shared, nrep=int(nr), runs=runs,
                       cohorts=[lab for lab, _ in C]), open(os.path.join(od, "info.json"), "w"), indent=1)
        log.append(f"{sub}\truns={len(runs)}\tsets={D.set_id.nunique()}\tshared_shifts={shared}\tnrep={nr}")
        print(log[-1], flush=True)
    open(os.path.join(a.out, "META_LOG.tsv"), "w").write("\n".join(log) + "\n")


if __name__ == "__main__":
    main()
