"""Fig. 3 joint model of the Sacks 2026 dosage paper, for one cohort: fits, in-cohort summaries and the meta-ready export.

One specification for every cohort (designs/fig3_spec_2026-11-23.json): UK Biobank, All of Us WGS, All of Us arrays. Uses the cohort's
pipeline config (the same JSON as run_pipeline.py: pheno_table, cnv_table, traits with their covariates and TWAS files, grid resources) plus
its "fig3" section: {"spec": ..., "perms_file": "<shared_shifts>/{trait}_....npy", "jobs": 4, "strata": ["all"], "export_dir": optional}.

  python run_fig3.py --config configs/<cohort>.json --steps fit,summary,export [--jobs 4] [--stratum all]

Steps
  filter   call set of the spec ("cnv_filter", Jon 2026-11-23): the PGC freeze filters, which every cohort's CNV table already carries
           (>= 10 kb; arrays >= 10 probes; calls carried by > 10% of samples at >= 50% reciprocal overlap removed). This step drops nothing:
           it records the per-position carrier frequency of the dataset, per dosage (cnv_freq_filter.py) in
           <results_dir>/<stratum>/fig3/cnv_filter/{position_freq.tsv.gz, run_info.json}. Run automatically by 'fit'. Only if the spec set
           "position_max_freq" would calls mostly above that frequency be dropped (cnv_rare.cnv) and every fit read the filtered table;
           the export records the setting and the meta-analysis refuses to mix cohorts with different settings.
  fit      conga_grid_qt.py for size2t (4 CNV-trait x TWAS-trait pairs), all2t and totdose (matched pairs), and the UNIVARIATE tests of each
           class in each CNV-size stratum (fits uni_N1_inf, uni_N1_1, uni_N2_inf: 7 one-class models each; these come before the joint model,
           Jon 2026-11-23); 2,000 shared shifts, LOO = drivers with the 8 fixed major loci + minor_all + no_loci; gene_mirror_qt.py per trait; fig3_dose_bins.py per trait (binned dose response of
           every stratum: all / single-gene / >= 2-gene CNVs x TWAS-by-constraint classes and total genes; analytic OLS). Finished fits are
           skipped (restartable).
           -> <results_dir>/<stratum>/fig3/<fit>/grid/<cnv>_x_<twas>/<variant>/, .../fig3/gene_mirror/<trait>/ and .../fig3/dose_bins/<trait>/
  summary  mirror4_summary.py (M4, armDEL, armDUP) and dose_summary.py (D4, doseOth, doseTot, mirTot) per fit, for QC inside the cohort.
  export   aggregate-only, meta-ready copy -> <export_dir>/<stratum>/fig3/ (default export_dir = the config's export_dir + "_fig3"):
           per fit observed_terms.tsv (t, beta, se of every term in every run), null_runs.npz (the t of every term in every shift replicate
           of every run, float64) and run_info.json; per trait gene_z.npz (per-gene score z of the gene mirror, with the shift offsets)
           and the dose-response bin tables (fig3_dose_bins.tsv, fig3_dose_shape.tsv: beta, SE, carriers and mean genes per fixed bin);
           the design, the spec and EXPORT_MANIFEST.txt. The summary tables are NOT exported: meta_fig3.py recomputes them.
           Small-cell rule (config min_cell = k; 0 = off): a term carried by 1..k people in a run loses beta / se / t in that run, its null
           column and its mirror readout; a gene hit in 1..k people (per run, CNV size and dosage) loses its z; gene carrier counts are
           exported only as 0 (not usable) or k + 1 (usable), never as counts; counts 1..k in run_info.json are written as "<=k"; a
           dose-response bin with 1..max(k, 20) people has no estimate, and the no_loci bins carry no carrier counts (fig3_dose_bins.py).
null model: genome shift of the gene annotation, 2,000 shifts shared with UK Biobank (the shift file's sha1 is checked against the spec before
and after every fit); gene mirror: DUP-vs-DEL gene-circle shift, 10,000 shifts, same seed in every cohort (gene_mirror_qt.py). Total-genes
terms: analytic t (no annotation, so no shift null).
"""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
INFO_KEEP = ["trait", "N", "panels", "n_panels", "q95_shet", "clingen_HI_genes_in_circle", "genes_twas", "variants", "nrep", "seed", "perms",
             "perms_sha1", "loo_mode", "major_loci", "major_loci_source", "calls_by_class"]
GM_KEEP = ["trait", "N", "genes", "KMIN", "nrep", "seed", "panels", "n_panels", "major_loci", "runs", "mincar"]


def expand(o):
    if isinstance(o, str):
        return os.path.expandvars(o)
    if isinstance(o, dict):
        return {k: expand(v) for k, v in o.items()}
    if isinstance(o, list):
        return [expand(v) for v in o]
    return o


def sha1_shifts(path, nrep):
    offs = np.load(os.path.join(HERE, path) if not os.path.isabs(path) else path)
    if len(offs) < nrep:
        raise SystemExit(f"{path} holds {len(offs)} shifts < {nrep}")
    return hashlib.sha1(np.asarray(offs[:nrep], dtype=np.int64).tobytes()).hexdigest()[:16]


def variants_of(f):
    return f.get("variants") or [f["variant"]]


def dir_of(f):
    return f.get("dir") or f["variant"]


def pairs_of(fit, traits):
    return [(c, t) for c in traits for t in traits] if fit["pairs"] == "all" else [(c, c) for c in traits]


def run(cmd, log, threads):
    env = dict(os.environ, OPENBLAS_NUM_THREADS=str(threads), OMP_NUM_THREADS=str(threads), MKL_NUM_THREADS=str(threads), PYTHONUNBUFFERED="1")
    print("+", " ".join(cmd), ">", log, flush=True)
    with open(log, "w") as fh:
        subprocess.run([sys.executable, *cmd], check=True, cwd=HERE, env=env, stdout=fh, stderr=subprocess.STDOUT)


def scrub(o, k):
    """Counts 1..k -> '<=k' (the All of Us egress rule), recursively."""
    if isinstance(o, dict):
        return {kk: scrub(v, k) for kk, v in o.items()}
    if isinstance(o, list):
        return [scrub(v, k) for v in o]
    if isinstance(o, (int, np.integer)) and not isinstance(o, bool) and 1 <= o <= k:
        return f"<={k}"
    return o


def readout_names(var):
    n = [f"{x}|DEL" for x in var["DEL"]] + [f"{x}|DUP" for x in var["DUP"]]
    return n + [f"{x}|mirror" for x in var["DEL"] if x in var["DUP"] and x != "n_all"]


def export_fit(src, dst, vars_, k):
    """One conga_grid_qt.py output folder (one or more variants) -> its export. Returns the number of (run, term) cells withheld."""
    os.makedirs(dst, exist_ok=True)
    names = {v["name"]: readout_names(v) for v in vars_}
    ob = pd.read_csv(os.path.join(src, "observed_terms.tsv"), sep="\t")
    nz = np.load(os.path.join(src, "null_runs.npz")); nul = {f: np.array(nz[f], dtype=np.float64) for f in nz.files}
    small = set()
    if k > 0:
        if "n_carriers_term" not in ob.columns:
            raise SystemExit(f"{src}: observed_terms.tsv has no n_carriers_term column; refit with the current conga_grid_qt.py (min_cell {k})")
        n = pd.to_numeric(ob.n_carriers_term, errors="coerce")
        for r in ob[n.between(1, k)].itertuples():
            small.add((r.run, r.readout)); small.add((r.run, r.readout.rsplit("|", 1)[0] + "|mirror"))
        hit = np.array([(r, q) in small for r, q in zip(ob.run, ob.readout)])
        for c in ["beta", "se", "t"]:
            ob.loc[hit, c] = np.nan
        ob["n_carriers_term"] = ob.n_carriers_term.astype(object)
        ob.loc[n.between(1, k).to_numpy(), "n_carriers_term"] = f"<={k}"
        for (rn, ro) in small:
            for vn, nm in names.items():
                key = f"{rn}__{vn}"
                if key in nul and ro in nm:
                    nul[key][:, nm.index(ro)] = np.nan
    keep = [c for c in ["run", "variant", "readout", "beta", "se", "t", "n_carriers_term"] if c in ob.columns]
    ob[keep].to_csv(os.path.join(dst, "observed_terms.tsv"), sep="\t", index=False, float_format="%.10g")
    np.savez_compressed(os.path.join(dst, "null_runs.npz"), **nul)
    info = json.load(open(os.path.join(src, "run_info.json")))
    info = {kk: info[kk] for kk in INFO_KEEP if kk in info}; info["readouts"] = names; info["min_cell"] = k
    if k > 0 and "calls_by_class" in info:
        info["calls_by_class"] = scrub(info["calls_by_class"], k)
    json.dump(info, open(os.path.join(dst, "run_info.json"), "w"), indent=1)
    return len({s for s in small if not s[1].endswith("|mirror")})


def export_gene_mirror(src, dst, k):
    os.makedirs(dst, exist_ok=True)
    z = np.load(os.path.join(src, "gene_z.npz")); out = {}; nmask = 0
    for f in z.files:
        if f.startswith("n__"):
            n = np.asarray(z[f], float)
            out[f] = np.where(n > k, k + 1, 0).astype(np.int16) if k > 0 else n     # never the counts themselves under a small-cell rule
        elif f.startswith(("zA__", "zB__")):
            v = np.array(z[f], dtype=np.float64)
            if k > 0:
                n = np.asarray(z["n__" + f.split("__", 1)[1]], float); m = n <= k; nmask += int((m & np.isfinite(v)).sum()); v[m] = np.nan
            out[f] = v
        else:
            out[f] = z[f]
    np.savez_compressed(os.path.join(dst, "gene_z.npz"), **out)
    info = json.load(open(os.path.join(src, "run_info.json"))); info = {kk: info[kk] for kk in GM_KEEP if kk in info}
    info["min_cell"] = k; info["usable_carriers"] = max(int(info.get("mincar", 5)), k + 1)
    json.dump(info, open(os.path.join(dst, "run_info.json"), "w"), indent=1)
    return nmask


def audit(root, k):
    """No count 1..k and no person-level file in the export."""
    bad = []
    for d, _, fs in os.walk(root):
        for f in fs:
            p = os.path.join(d, f)
            if any(x in f for x in ("pheno", "genotype", ".cnv")):
                bad.append(f"person-level file name: {p}")
            if f == "observed_terms.tsv" and k > 0:
                n = pd.to_numeric(pd.read_csv(p, sep="\t").get("n_carriers_term"), errors="coerce")
                if n.between(1, k).any():
                    bad.append(f"count 1..{k} in {p}")
            if f == "fig3_dose_bins.tsv" and k > 0:
                d = pd.read_csv(p, sep="	"); n = pd.to_numeric(d.n, errors="coerce")
                if n.between(1, k).any() or (d.run.ne("full") & n.notna() & d.side.ne("none")).any():
                    bad.append(f"small or no_loci bin counts in {p}")
            if f == "position_freq.tsv.gz" and k > 0 and pd.read_csv(p, sep="\t").n_carriers.between(1, k).any():
                bad.append(f"position counts 1..{k} in {p}")
            if f == "gene_z.npz" and k > 0:
                z = np.load(p)
                for a in [a for a in z.files if a.startswith("n__")]:
                    if not np.isin(z[a], [0, k + 1]).all():
                        bad.append(f"gene carrier counts in {p}:{a}")
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True); ap.add_argument("--steps", default="fit,summary")
    ap.add_argument("--stratum", default="all"); ap.add_argument("--jobs", type=int, default=None)
    ap.add_argument("--only_fit", default="", help="comma list of fit names (size2t_mirror, all2t_mirror, totdose, gene_mirror, dose_bins)")
    ap.add_argument("--min_cell", type=int, default=None, help="override the config's min_cell for the export (0 = no suppression)")
    a = ap.parse_args()
    C = expand(json.load(open(a.config))); F = C.get("fig3")
    if not F:
        raise SystemExit('the config has no "fig3" section')
    S = json.load(open(os.path.join(HERE, F["spec"]))); G = C["grid"]; steps = a.steps.split(",")
    dfile = lambda f: f.get("design", S["design"])
    DES = {p_: {v["name"]: v for v in json.load(open(os.path.join(HERE, p_)))} for p_ in sorted({dfile(f) for f in S["fits"]})}
    traits = [t for t in ("height", "bmi") if t in C["traits"]]
    strat = a.stratum if C.get("strata") else ""
    root = os.path.join(C["results_dir"], strat, "fig3")
    pheno = C["pheno_table"] if (not strat or not C["strata"].get(strat)) else os.path.join(C["results_dir"], strat, "stratum_pheno.tsv")
    only = set(x for x in a.only_fit.split(",") if x)
    fits = [f for f in S["fits"] if not only or f["name"] in only]
    pf = {t: F["perms_file"].replace("{trait}", t) for t in traits}
    cpu = (len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count()) or 1
    jobs = a.jobs or F.get("jobs", 4); thr = max(1, cpu // jobs)
    CF = S.get("cnv_filter") or {}; pmax = CF.get("position_max_freq")       # None: PGC freeze filters only (Jon 2026-11-23)
    FL = pmax is not None                                                    # True: the optional position-frequency filter is on
    fdir = os.path.join(root, "cnv_filter"); cnv_use = os.path.join(fdir, "cnv_rare.cnv") if FL else C["cnv_table"]
    if "filter" in steps or ("fit" in steps and (FL or not only)):
        if os.path.exists(os.path.join(fdir, "run_info.json")):
            print("skip (done):", fdir, flush=True)
        else:
            import glob
            old = glob.glob(os.path.join(root, "*", "grid", "*", "*", "run_info.json")) + glob.glob(os.path.join(root, "gene_mirror", "*", "run_info.json")) \
                + glob.glob(os.path.join(root, "dose_bins", "*", "run_info.json"))
            if old and FL:               # fits made on the unfiltered table would be skipped as 'done' and mixed with filtered ones
                raise SystemExit(f"{root} holds {len(old)} fits made WITHOUT the rare-CNV filter. Move that folder away "
                                 f"(e.g. mv {root} {root}_unfiltered) and rerun: every Fig. 3 fit must use the filtered calls.")
            os.makedirs(os.path.join(root, "logs"), exist_ok=True)
            run(["cnv_freq_filter.py", "--cnv", C["cnv_table"], "--pheno", pheno, "--min_frac", str(CF.get("min_frac", 0.5)), "--out_dir", fdir]
                + (["--max_freq", str(pmax), "--out_cnv", os.path.join(fdir, "cnv_rare.tmp")] if FL else []), os.path.join(root, "logs", "cnv_filter.log"), cpu)
            if FL:
                os.replace(os.path.join(fdir, "cnv_rare.tmp"), cnv_use)
        fi = json.load(open(os.path.join(fdir, "run_info.json")))
        print(f"call set: {CF.get('rule', 'as the CNV table')}; position filter {'off' if not FL else pmax}; calls {fi['calls_in']} -> {fi['calls_out']}; "
              f"highest per-position carrier frequency {fi['max_position_freq']}", flush=True)

    if "fit" in steps:
        for t in traits:                # the cohort must load the UK Biobank shift list of the CNV trait
            h = sha1_shifts(pf[t], S["nrep"])
            if h != S["perms_sha1"][t]:
                raise SystemExit(f"shift list {pf[t]}: sha1 {h} != {S['perms_sha1'][t]} (UK Biobank). Stop and tell Jon.")
            print(f"shift list {t}: {os.path.basename(pf[t])} sha1 {h} OK", flush=True)
        os.makedirs(os.path.join(root, "logs"), exist_ok=True); calls = []
        for f in fits:
            for c, t in pairs_of(f, traits):
                T = C["traits"][c]; out = os.path.join(root, f["name"], "grid", f"{c}_x_{t}", dir_of(f))
                if os.path.exists(os.path.join(out, "run_info.json")):
                    print("skip (done):", out, flush=True); continue
                calls.append((["conga_grid_qt.py", "--pheno", pheno, "--trait", T["column"], "--cnv", cnv_use, "--twas", C["traits"][t]["twas"],
                               "--panels", S["panels"], "--exons", G["exons"], "--gene_cov", G["gene_cov"], "--clingen", G["clingen"],
                               "--loci", G["loci"], "--loo_runs", G["loo_runs"], "--design", dfile(f), "--only", ",".join(variants_of(f)),
                               "--nrep", str(S["nrep"]), "--seed", str(S["seed"]), "--loo_mode", S["loo_mode"], "--major_loci", S["major_loci"],
                               "--perms_file", pf[c], "--numeric", ",".join(T["numeric"]), "--factors", ",".join(T["factors"]), "--out", out],
                              os.path.join(root, "logs", f"fit_{f['name']}_{c}_x_{t}.log")))
        if not only or "gene_mirror" in only:
            gm = S["gene_mirror"]
            for t in traits:
                T = C["traits"][t]; out = os.path.join(root, "gene_mirror", t)
                if os.path.exists(os.path.join(out, "run_info.json")):
                    print("skip (done):", out, flush=True); continue
                calls.append((["gene_mirror_qt.py", "--pheno", pheno, "--trait", T["column"], "--cnv", cnv_use, "--twas", T["twas"],
                               "--panels", S["panels"], "--exons", G["exons"], "--gene_cov", G["gene_cov"], "--loci", G["loci"],
                               "--loo_runs", G["loo_runs"], "--major_loci", S["major_loci"], "--nrep", str(gm["nrep"]), "--seed", str(gm["seed"]),
                               "--numeric", ",".join(T["numeric"]), "--factors", ",".join(T["factors"]), "--out", out],
                              os.path.join(root, "logs", f"gene_mirror_{t}.log")))
        if S.get("dose_bins") and (not only or "dose_bins" in only):
            kb = max(int(C.get("min_cell", 20)), int(S["dose_bins"].get("min_cell_floor", 20)))
            for t in traits:
                T = C["traits"][t]; out = os.path.join(root, "dose_bins", t)
                if os.path.exists(os.path.join(out, "run_info.json")):
                    print("skip (done):", out, flush=True); continue
                calls.append((["fig3_dose_bins.py", "--pheno", pheno, "--trait", T["column"], "--cnv", cnv_use, "--twas", T["twas"],
                               "--panels", S["panels"], "--exons", G["exons"], "--gene_cov", G["gene_cov"], "--loci", G["loci"],
                               "--loo_runs", G["loo_runs"], "--numeric", ",".join(T["numeric"]), "--factors", ",".join(T["factors"]),
                               "--min_cell", str(kb), "--cohort", str(C.get("cohort", "cohort")), "--out", out],
                              os.path.join(root, "logs", f"dose_bins_{t}.log")))
        with ThreadPoolExecutor(max_workers=jobs) as ex:
            list(ex.map(lambda c: run(c[0], c[1], thr), calls))
        for f in fits:                  # every fit used the shared shifts and the fixed major loci
            for c, t in pairs_of(f, traits):
                i = json.load(open(os.path.join(root, f["name"], "grid", f"{c}_x_{t}", dir_of(f), "run_info.json")))
                assert i["perms_sha1"] == S["perms_sha1"][c] and i["nrep"] == S["nrep"], (f["name"], c, t, i["perms_sha1"])
                assert ",".join(i["major_loci"]) == S["major_loci"], (f["name"], c, t, i["major_loci"])
        print("fits done:", root, flush=True)

    if "summary" in steps:
        for f in [f for f in fits if f.get("summaries", True)]:     # the univariate fits need none: their grid_summary.tsv is the summary
            d = os.path.join(root, f["name"])
            for scr, o, extra in (("mirror4_summary.py", "mirror4_cross_summary.tsv", []), ("dose_summary.py", "dose_summary.tsv", ["--variant", f["variant"]])):
                run([scr, "--grid", os.path.join(d, "grid"), "--design", dfile(f), "--out", os.path.join(d, o)] + extra,
                    os.path.join(d, o.replace(".tsv", ".log")), cpu)
        print("summaries written (z and P only; for QC, not exported):", root, flush=True)

    if "export" in steps:
        k = C.get("min_cell", 20) if a.min_cell is None else a.min_cell
        E = os.path.join(F.get("export_dir", C["export_dir"].rstrip("/") + "_fig3"), strat, "fig3")
        if os.path.exists(E):
            raise SystemExit(f"{E} exists; move it away first (an export is never overwritten)")
        nterm = ngene = 0
        for f in fits:
            for c, t in pairs_of(f, traits):
                rel = os.path.join(f["name"], "grid", f"{c}_x_{t}", dir_of(f))
                nterm += export_fit(os.path.join(root, rel), os.path.join(E, rel), [DES[dfile(f)][v] for v in variants_of(f)], k)
        for t in traits:
            if os.path.exists(os.path.join(root, "gene_mirror", t, "gene_z.npz")):
                ngene += export_gene_mirror(os.path.join(root, "gene_mirror", t), os.path.join(E, "gene_mirror", t), k)
            db = os.path.join(root, "dose_bins", t)             # already small-cell masked when written (fig3_dose_bins.py --min_cell)
            if os.path.exists(os.path.join(db, "run_info.json")):
                kb = json.load(open(os.path.join(db, "run_info.json")))["min_cell"]
                if kb < k:
                    raise SystemExit(f"{db} was written with min_cell {kb} < {k}: delete it and rerun the fit step")
                os.makedirs(os.path.join(E, "dose_bins", t))
                for fn in ("fig3_dose_bins.tsv", "fig3_dose_shape.tsv", "run_info.json"):
                    shutil.copy(os.path.join(db, fn), os.path.join(E, "dose_bins", t, fn))
        if os.path.exists(os.path.join(fdir, "run_info.json")):      # call-set record: counts scrubbed, position frequencies above the small-cell limit only
            os.makedirs(os.path.join(E, "cnv_filter"))
            fi = json.load(open(os.path.join(fdir, "run_info.json")))
            json.dump(scrub(fi, k) if k > 0 else fi, open(os.path.join(E, "cnv_filter", "run_info.json"), "w"), indent=1)
            pf = pd.read_csv(os.path.join(fdir, "position_freq.tsv.gz"), sep="\t")
            pf[pf.n_carriers > k].to_csv(os.path.join(E, "cnv_filter", "position_freq.tsv.gz"), sep="\t", index=False, float_format="%.6g")
        os.makedirs(os.path.join(E, "designs"), exist_ok=True)
        for p in list(DES) + [F["spec"]]:
            shutil.copy(os.path.join(HERE, p), os.path.join(E, "designs", os.path.basename(p)))
        json.dump(dict(cohort=C.get("cohort"), stratum=strat or "all", min_cell=k, spec=os.path.basename(F["spec"]),
                       cnv_filter=dict(rule=CF.get("rule", "PGC freeze filters"), position_max_freq=pmax),
                       design_sha1=hashlib.sha1(b"".join(json.dumps(DES[p_], sort_keys=True).encode() for p_ in DES)).hexdigest()[:16],
                       spec_sha1=hashlib.sha1(json.dumps({k_: v for k_, v in S.items() if k_ not in ("note", "cnv_filter")}, sort_keys=True).encode()).hexdigest()[:16],
                       traits=traits, fits=[f["name"] for f in fits]), open(os.path.join(E, "FIG3_EXPORT_INFO.json"), "w"), indent=1)
        bad = audit(E, k); tot = 0; lines = []
        for d, _, fs in os.walk(E):
            for fn in sorted(fs):
                p = os.path.join(d, fn); sz = os.path.getsize(p); tot += sz; lines.append(f"{sz}\t{os.path.relpath(p, E)}")
        with open(os.path.join(E, "EXPORT_MANIFEST.txt"), "w") as fh:
            fh.write(f"Fig. 3 joint model export; small-cell rule: counts 1..{k} withheld; (run, term) cells withheld {nterm}; gene z withheld {ngene}\n"
                     f"audit: {'OK' if not bad else 'FAILED'}\ntotal bytes {tot} ({tot / 1e6:.1f} MB) in {len(lines)} files\nbytes\tfile\n" + "\n".join(lines) + "\n")
        print(f"exported to {E}: {tot / 1e6:.1f} MB, {len(lines)} files; term cells withheld {nterm}; gene z withheld {ngene}")
        print("\n".join("AUDIT VIOLATION " + b for b in bad) or f"AUDIT OK: no count 1..{k}, no person-level file")
        if bad:
            raise SystemExit("export audit failed: do not download this folder")


if __name__ == "__main__":
    main()
