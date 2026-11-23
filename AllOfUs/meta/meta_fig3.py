"""Meta-analysis of the Fig. 3 joint model (Sacks 2026 dosage paper) across cohorts, at UCSD, from the aggregate exports of run_fig3.py.

  python meta_fig3.py --cohort UKBB=<UKBB export>/all/fig3 [--cohort AoU=<AoU export>/all/fig3 ...] --out <new meta dir> [--code ../cnv_anthro]

Output: the folders and tables figures/figs_dosage_story.py reads for Fig. 3, so the figure is redrawn from the meta-analysis with
GW2_FIG3_META=<out> (one cohort alone reproduces that cohort's own tables):
  <out>/size2t_mirror/  grid/<cnv>_x_<twas>/size2t/{observed_terms.tsv, null_runs.npz, grid_summary.tsv, per_run_z.tsv, run_info.json},
                        mirror4_cross_summary.tsv (M4, armDEL, armDUP), dose_summary.tsv (D4, doseOth, doseTot, mirTot)
  <out>/all2t_mirror/   the same for the model pooled over CNV size        <out>/totdose/   total genes hit only
  <out>/gene_mirror/<trait>/gene_mirror_summary.tsv, gene_z.npz            <out>/meta_info.json, meta_checks.tsv
How cohorts are combined
  joint-model terms   per run (full, each major-locus removal, minor_all, no_loci) and term: Z = sum_c w_c t_c / sqrt(sum_c w_c^2), w_c = sqrt(N_c),
                      over the cohorts with a value (a term withheld by a cohort's small-cell rule is left out for that cohort); the same
                      formula in every shift replicate. The meta tables are then written in the single-cohort format and the single-cohort
                      summary scripts (mirror4_summary.py, dose_summary.py) run on them unchanged, so M4, D4, arm and total contrasts, their
                      shift-null z, the 9 LOO runs and no_loci are all recomputed from the meta t.
  gene mirror         per gene, run, CNV size and dosage: the same sqrt(N)-weighted Stouffer of the per-gene score z, over the cohorts where the
                      gene is usable (carriers >= 5, and above the cohort's small-cell limit); a gene enters r(zDEL, zDUP) when both dosages
                      have a usable cohort. total / residual / explained, LOO and no_loci then as gene_mirror_qt.py.
Checks (meta_checks.tsv; any failure stops the run): same design and spec, same readouts, same shift list (perms_sha1, nrep) and same major
loci in every cohort; same gene order and gene-circle shifts for the gene mirror.
null model: joint terms: genome shift of the gene annotation, the SAME 2,000 UK Biobank shifts in every cohort, so replicate i is one shift
everywhere and the meta statistic is formed replicate by replicate. Gene mirror: DUP-vs-DEL gene-circle shift (10,000, same offsets in every
cohort), applied to the meta per-gene z. Total-genes terms: analytic (Stouffer of the cohorts' t).
"""
import argparse
import json
import os
import subprocess
import sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
MINCAR = 5


def stouffer(X, w):
    """X: cohorts x ... (NaN = cohort has no value); sqrt-N weights over the cohorts with a value."""
    if X.shape[0] == 1:                 # one cohort: its own statistic, bit for bit
        return np.array(X[0], dtype=float)
    ok = np.isfinite(X); ww = w.reshape((-1,) + (1,) * (X.ndim - 1)) * ok
    den = np.sqrt((ww ** 2).sum(0))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, (np.where(ok, X, 0.0) * ww).sum(0) / den, np.nan)


def grid_summary(trait, panels, vn, names, ot, nul, runs, loo, nrep, majors):
    """The summary block of conga_grid_qt.py (same columns), from observed t (Series indexed run, readout) and null arrays per run."""
    out, per = [], []
    for j, nm in enumerate(names):
        def emp(rn):
            o = ot.get((rn, nm), np.nan)
            if rn not in nul:
                return o, np.nan, np.nan
            nu = nul[rn][:, j]; nu = nu[np.isfinite(nu)]
            if not np.isfinite(o) or len(nu) < 10:
                return o, np.nan, np.nan
            mu, sd = nu.mean(), nu.std(ddof=1)
            return o, (o - mu) / sd, (1 + np.sum(np.abs(nu - mu) >= abs(o - mu))) / (len(nu) + 1)
        o, ze, pe = emp("full"); s = np.sign(ze)
        lz = {rn: emp(rn)[1] for rn in loo}
        for rn in runs:
            if (rn, nm) in ot.index:
                per.append(dict(variant=vn, readout=nm, run=rn, t=ot[(rn, nm)], z_emp=emp(rn)[1]))
        held = sum(1 for q in lz.values() if np.isfinite(q) and abs(q) > 1.96 and np.sign(q) == s)
        rev = sum(1 for q in lz.values() if np.isfinite(q) and np.sign(q) == -s)
        fin = {rn: q for rn, q in lz.items() if np.isfinite(q)}
        wk = min(fin, key=lambda rn: s * fin[rn]) if fin else ""
        if loo:
            T_obs = np.array([ot.get((rn, nm), np.nan) for rn in loo]); T_nul = np.column_stack([nul[rn][:, j] for rn in loo])
            with np.errstate(all="ignore"):
                mo = np.nanmin(s * T_obs) if np.isfinite(T_obs).any() else np.nan
                mn = np.nanmin(s * T_nul, axis=1) if np.isfinite(T_nul).any() else np.full(len(T_nul), np.nan)
            pmin = (1 + np.sum(mn >= mo)) / (len(mn) + 1)
        else:
            pmin = np.nan
        o2, z2, p2 = emp("no_loci"); o3, z3, p3 = emp("minor_all")
        out.append(dict(trait=trait, panels=panels, variant=vn, readout=nm, t_analytic=o, z_emp=ze, P_emp=pe, loo_held=held, loo_n=len(loo),
                        loo_reversed=rev, weakest=wk.replace("loo_", ""), weakest_z=fin.get(wk, np.nan), P_minLOO_matched=pmin,
                        noloci_z_emp=z2, noloci_P_emp=p2, nrep=nrep, loo_mode="drivers", minor_all_z_emp=z3, minor_all_P_emp=p3,
                        major_loci=";".join(majors)))
    return pd.DataFrame(out), pd.DataFrame(per)


def meta_fit(C, rel, out, vn, checks):
    """C: {label: fig3 export root}. rel: <fit>/grid/<pair>/<variant>. Writes the meta fit folder in the single-cohort format."""
    info, obs, nul = {}, {}, {}
    for lab, root in C.items():
        d = os.path.join(root, rel)
        info[lab] = json.load(open(os.path.join(d, "run_info.json")))
        o = pd.read_csv(os.path.join(d, "observed_terms.tsv"), sep="\t"); o["t"] = pd.to_numeric(o.t, errors="coerce")
        obs[lab] = o.drop_duplicates(["run", "readout"]).set_index(["run", "readout"])
        nul[lab] = np.load(os.path.join(d, "null_runs.npz"))
    labs = list(C); i0 = info[labs[0]]
    names = i0.get("readouts") or list(obs[labs[0]].loc["full"].index)
    same = lambda f: len({json.dumps(f(info[l]), sort_keys=True) for l in labs}) == 1
    ck = dict(fit=rel, cohorts=";".join(labs), same_shifts=same(lambda i: [i.get("perms_sha1"), i.get("nrep")]),
              same_major_loci=same(lambda i: i.get("major_loci")), same_panels=same(lambda i: i.get("panels")),
              same_readouts=all(list(obs[l].loc["full"].index) == list(obs[labs[0]].loc["full"].index) for l in labs))
    checks.append(ck)
    if not all(v for k_, v in ck.items() if k_.startswith("same_")):
        raise SystemExit(f"cohorts are not comparable for {rel}: {ck}")
    w = np.array([np.sqrt(float(info[l]["N"])) for l in labs])
    runs = [f.split("__", 1)[0] for f in nul[labs[0]].files]
    runs = [r for r in runs if all(f"{r}__{vn}" in nul[l].files for l in labs)]
    nrep = min(nul[l][f"full__{vn}"].shape[0] for l in labs)
    mnul = {r: stouffer(np.stack([np.asarray(nul[l][f"{r}__{vn}"], float)[:nrep] for l in labs]), w) for r in runs}
    rows = []
    all_runs = list(dict.fromkeys(obs[labs[0]].index.get_level_values(0)))
    for r in all_runs:                  # every observed run (the null runs, plus the observed-only single-locus removals)
        if not all(r in obs[l].index.get_level_values(0) for l in labs):
            continue
        T = np.array([[obs[l].t.get((r, nm), np.nan) for nm in names] for l in labs])
        z = stouffer(T, w)
        for j, nm in enumerate(names):
            row = dict(run=r, variant=vn, readout=nm, t=z[j], n_cohorts=int(np.isfinite(T[:, j]).sum()))
            row.update({f"t_{l}": T[i, j] for i, l in enumerate(labs)}); rows.append(row)
    O = pd.DataFrame(rows); d = os.path.join(out, rel); os.makedirs(d, exist_ok=True)
    O.to_csv(os.path.join(d, "observed_terms.tsv"), sep="\t", index=False, float_format="%.10g")
    np.savez_compressed(os.path.join(d, "null_runs.npz"), **{f"{r}__{vn}": a for r, a in mnul.items()})
    majors = list(i0["major_loci"]); loo = [f"loo_{m}" for m in majors if f"loo_{m}" in mnul]
    gs, per = grid_summary(i0["trait"], i0.get("panels"), vn, names, O.set_index(["run", "readout"]).t, mnul, all_runs, loo, nrep, majors)
    gs.to_csv(os.path.join(d, "grid_summary.tsv"), sep="\t", index=False, float_format="%.6g")
    per.to_csv(os.path.join(d, "per_run_z.tsv"), sep="\t", index=False, float_format="%.6g")
    json.dump(dict(trait=i0["trait"], N=int(sum(float(info[l]["N"]) for l in labs)), N_by_cohort={l: info[l]["N"] for l in labs}, cohorts=labs,
                   weights="sqrt(N)", panels=i0.get("panels"), nrep=int(nrep), perms_sha1=i0.get("perms_sha1"), loo_mode="drivers",
                   major_loci=majors, variants=[vn], readouts=names), open(os.path.join(d, "run_info.json"), "w"), indent=1)


def gene_mirror(C, trait, out, checks):
    Z = {l: np.load(os.path.join(r, "gene_mirror", trait, "gene_z.npz")) for l, r in C.items()}
    I = {l: json.load(open(os.path.join(r, "gene_mirror", trait, "run_info.json"))) for l, r in C.items()}
    labs = list(C); z0 = Z[labs[0]]
    ck = dict(fit=f"gene_mirror/{trait}", cohorts=";".join(labs),
              same_genes=all(np.array_equal(Z[l]["entrez"], z0["entrez"]) and np.array_equal(Z[l]["hc"], z0["hc"]) for l in labs),
              same_shifts=all(np.array_equal(Z[l]["KS"], z0["KS"]) for l in labs),
              same_major_loci=len({json.dumps(I[l]["major_loci"]) for l in labs}) == 1,
              same_twas_sign=all(np.array_equal(Z[l]["twas_sign"], z0["twas_sign"]) for l in labs))
    checks.append(ck)
    if not all(v for k_, v in ck.items() if k_.startswith("same_")):
        raise SystemExit(f"cohorts are not comparable for the gene mirror ({trait}): {ck}")
    w = np.array([np.sqrt(float(I[l]["N"])) for l in labs]); hc = z0["hc"].astype(bool); KS = z0["KS"]; M = len(hc)
    runs = I[labs[0]]["runs"]; loo = [r for r in runs if r.startswith("loo_")] + ["minor_all"]
    SIZES = ["all sizes", "N1_1"]; STRATA = {"all": np.ones(M, bool), "LC": ~hc, "HC": hc}

    def meta_z(mdl, rn, sz, dos):
        X = []
        for l in labs:
            z = np.array(Z[l][f"z{mdl}__{rn}__{sz}__{dos}"], float); n = np.asarray(Z[l][f"n__{rn}__{sz}__{dos}"], float)
            z[~(n >= max(MINCAR, int(I[l].get("mincar", MINCAR))))] = np.nan; X.append(z)
        return stouffer(np.stack(X), w)

    idx0 = np.arange(M)

    def rr(zd, zu, keep, k=0):
        b = zu[(idx0 + k) % M]; ok = keep & np.isfinite(zd) & np.isfinite(b)
        return (np.corrcoef(zd[ok], b[ok])[0, 1], int(ok.sum())) if ok.sum() > 2 else (np.nan, int(ok.sum()))

    def emp(o, nu):
        nu = nu[np.isfinite(nu)]
        if not np.isfinite(o) or len(nu) < 10:
            return np.nan, np.nan
        mu, sd = nu.mean(), nu.std(ddof=1)
        return (o - mu) / sd, (1 + np.sum(np.abs(nu - mu) >= abs(o - mu))) / (len(nu) + 1)

    res = {}; save = {}
    for rn in runs:
        for sz in SIZES:
            zz = {(m, d): meta_z(m, rn, sz, d) for m in "AB" for d in ("DEL", "DUP")}
            for (m, d), v in zz.items():
                save[f"z{m}__{rn}__{sz}__{d}"] = v
                if m == "A":
                    save[f"n__{rn}__{sz}__{d}"] = np.where(np.isfinite(v), MINCAR, 0).astype(np.int16)   # usable flag, not a count
            for st, keep in STRATA.items():
                o = {}; nu = {}; ng = 0
                for m in "AB":
                    o[m], ng = rr(zz[(m, "DEL")], zz[(m, "DUP")], keep)
                    nu[m] = np.array([rr(zz[(m, "DEL")], zz[(m, "DUP")], keep, int(k))[0] for k in KS])
                res[(rn, sz, st)] = dict(n_genes=ng, total=(o["A"], nu["A"]), residual=(o["B"], nu["B"]), explained=(o["A"] - o["B"], nu["A"] - nu["B"]))
        print(f"gene mirror {trait}: run {rn} done", flush=True)
    rows = []
    for sz in SIZES:
        for st in STRATA:
            R = res[("full", sz, st)]
            for stat in ("total", "residual", "explained"):
                o, nu = R[stat]; z, p = emp(o, nu); s = np.sign(z)
                lz = {rn: emp(*res[(rn, sz, st)][stat])[0] for rn in loo}
                held = sum(1 for q in lz.values() if np.isfinite(q) and abs(q) > 1.96 and np.sign(q) == s)
                fin = {k_: q for k_, q in lz.items() if np.isfinite(q)}; wk = min(fin, key=lambda k_: s * fin[k_]) if fin else ""
                zn, pn = emp(*res[("no_loci", sz, st)][stat]); tot = R["total"][0]
                rows.append(dict(trait=I[labs[0]]["trait"], size=sz, stratum=st, stat=stat, obs=o, null_mean=np.nanmean(nu), z=z, P=p, n_genes=R["n_genes"],
                                 share=(o / tot if stat == "explained" and np.isfinite(tot) and tot != 0 else np.nan), loo_held=held, n_loo=len(loo),
                                 weakest=wk.replace("loo_", ""), z_weakest=fin.get(wk, np.nan), noloci_obs=res[("no_loci", sz, st)][stat][0],
                                 noloci_z=zn, noloci_P=pn, n_genes_noloci=res[("no_loci", sz, st)]["n_genes"], cohorts=";".join(labs)))
    d = os.path.join(out, "gene_mirror", trait); os.makedirs(d, exist_ok=True)
    pd.DataFrame(rows).to_csv(os.path.join(d, "gene_mirror_summary.tsv"), sep="\t", index=False, float_format="%.6g")
    np.savez_compressed(os.path.join(d, "gene_z.npz"), entrez=z0["entrez"], hc=z0["hc"], twas_sign=z0["twas_sign"], KS=KS, **save)
    json.dump(dict(trait=I[labs[0]]["trait"], N=int(sum(float(I[l]["N"]) for l in labs)), N_by_cohort={l: I[l]["N"] for l in labs}, cohorts=labs,
                   weights="sqrt(N)", genes=int(M), nrep=int(len(KS)), major_loci=I[labs[0]]["major_loci"], runs=runs, mincar=MINCAR,
                   usable_carriers_by_cohort={l: I[l].get("usable_carriers", I[l].get("mincar")) for l in labs}),
              open(os.path.join(d, "run_info.json"), "w"), indent=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", action="append", required=True, help="LABEL=<export>/<stratum>/fig3 (repeat; the first is the reference)")
    ap.add_argument("--out", required=True); ap.add_argument("--code", default=os.path.normpath(os.path.join(HERE, "..", "cnv_anthro")))
    ap.add_argument("--parts", default="joint,gene_mirror")
    a = ap.parse_args(); C = dict(x.split("=", 1) for x in a.cohort); labs = list(C); parts = a.parts.split(",")
    if os.path.exists(a.out):
        raise SystemExit(f"{a.out} exists; choose a new output folder")
    os.makedirs(a.out)
    E = {l: json.load(open(os.path.join(r, "FIG3_EXPORT_INFO.json"))) for l, r in C.items()}
    if len({e["design_sha1"] for e in E.values()}) != 1 or len({e["spec"] for e in E.values()}) != 1:
        raise SystemExit(f"cohorts ran different designs or specs: { {l: (e['spec'], e['design_sha1']) for l, e in E.items()} }")
    ref = C[labs[0]]; spec = json.load(open(os.path.join(ref, "designs", E[labs[0]]["spec"])))
    design = os.path.join(ref, "designs", os.path.basename(spec["design"])); checks = []
    if "joint" in parts:
        for f in spec["fits"]:
            g = os.path.join(ref, f["name"], "grid")
            pairs = [p for p in sorted(os.listdir(g)) if all(os.path.isdir(os.path.join(r, f["name"], "grid", p, f["variant"])) for r in C.values())] if os.path.isdir(g) else []
            for p in pairs:
                meta_fit(C, os.path.join(f["name"], "grid", p, f["variant"]), a.out, f["variant"], checks)
            if not pairs:
                print("skip (not in every cohort):", f["name"]); continue
            d = os.path.join(a.out, f["name"])
            for scr, o, extra in (("mirror4_summary.py", "mirror4_cross_summary.tsv", []), ("dose_summary.py", "dose_summary.tsv", ["--variant", f["variant"]])):
                with open(os.path.join(d, o.replace(".tsv", ".log")), "w") as fh:
                    subprocess.run([sys.executable, os.path.join(a.code, scr), "--grid", os.path.join(d, "grid"), "--design", design,
                                    "--out", os.path.join(d, o)] + extra, check=True, stdout=fh, stderr=subprocess.STDOUT)
            print("meta", f["name"], pairs, flush=True)
    if "gene_mirror" in parts:
        for t in ("height", "bmi"):
            if all(os.path.exists(os.path.join(r, "gene_mirror", t, "gene_z.npz")) for r in C.values()):
                gene_mirror(C, t, a.out, checks)
    pd.DataFrame(checks).to_csv(os.path.join(a.out, "meta_checks.tsv"), sep="\t", index=False)
    name = {"UKBB": "UK Biobank", "AoU": "All of Us", "AoU_WGS": "All of Us WGS", "AoU_array": "All of Us arrays"}
    json.dump(dict(cohorts=labs, label=" + ".join(name.get(l, l) for l in labs) + (" (meta-analysis)" if len(labs) > 1 else ""),
                   exports=C, export_info=E, parts=parts), open(os.path.join(a.out, "meta_info.json"), "w"), indent=1)
    print("wrote", a.out, "| cohorts:", labs)


if __name__ == "__main__":
    main()
