"""Single entry point. One code base, one JSON config per cohort (configs/ukbb_expanse.json, configs/aou_workbench.json).

  python run_pipeline.py --config configs/<cohort>.json --steps prep,gwas,plots,alleles,loci,burden,grid,additivity,scores,joint,export
        [--only_trait height] [--only_panels gtex_nonbrain] [--strata all,eur,female] [--jobs 8]

Steps
  prep     cohort adapters (UKBB: prep_ukbb.py; AoU: aou_sv_to_cnv.py + aou_pheno_extract.py)
  gwas     probe-level CNV GWAS (cnv_gwas_qt.py); with gwas.ref_probes also at the reference cohort's probe positions (meta-analysis)
  plots    Miami plot per trait
  alleles  Table S9 allele tests (alleles.py)
  burden   genome-wide deletion/duplication burden: signed-burden bins, U-shape test, removal of monogenic-like duplication carriers
  loci     per-locus DEL/DUP effects at the growth loci (locus_effects.py; Sacks recurrent intervals + UKBB-discovered loci)
  axes     ConGA consensus2 axes (conga_axes_qt.py)
  grid     genome-wide univariate grid, call-size split and joint models (conga_grid_qt.py), shared shifts for meta-analysis
  additivity  constrained-gene dose response, sub-additivity, asymmetric vs mirror duplications, single hit vs count (additivity_constraint.py)
  scores   cross-fitted mirror/constraint axis scores and their additivity across CNVs (score_axes.py)
  joint    joint models behind the joint-model heatmap: all CNVs + 4 carrier removals (make_union_strata.py + conga_grid_qt.py, UKBB shifts)
  export   aggregate-only copy with the small-cell rule, including cross-stratum (secondary) suppression (export_aggregate.py)

Strata (config "strata"): each stratum is a row filter on the phenotype table, e.g. {"eur": {"ancestry_pred": ["eur"]}, "female": {"sex": [2]}}.
Within a stratum the traits are re-ranked (RINT within sex in that stratum) and every step writes to <results_dir>/<stratum>/.
The cluster (SLURM on Expanse) or the Workbench (Jupyter cell / dsub task) only decides *where* this runs. Paths may use ${ENV} variables.
"""
import argparse
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from common import rint_within  # noqa: E402


def expand(o):
    if isinstance(o, str):
        return os.path.expandvars(o)
    if isinstance(o, dict):
        return {k: expand(v) for k, v in o.items()}
    if isinstance(o, list):
        return [expand(v) for v in o]
    return o


def run(args, threads=None):
    cmd = [sys.executable, *args]
    env = dict(os.environ)
    if threads:
        env.update(OPENBLAS_NUM_THREADS=str(threads), OMP_NUM_THREADS=str(threads), MKL_NUM_THREADS=str(threads))
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, cwd=HERE, env=env)


def make_stratum(pheno, filt, raw_of, out_path):
    """Filter the phenotype table and re-rank the traits within sex inside the stratum."""
    P = pd.read_csv(pheno, sep="\t", dtype={"IID": str})
    for col, vals in (filt or {}).items():
        P = P[P[col].astype(str).isin([str(v) for v in vals])]
    for rcol, raw in raw_of.items():
        if raw in P.columns:
            P[rcol] = rint_within(P[raw].to_numpy(float), P["sex"].to_numpy())
    tmp = f"{out_path}.{os.getpid()}.tmp"                 # atomic: parallel tasks of one stratum write identical content
    P.to_csv(tmp, sep="\t", index=False); os.replace(tmp, out_path)
    return len(P)


def safe(n):
    return n.replace(":", "_").replace("|", "__").replace("~", "-")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True); ap.add_argument("--steps", default="prep,gwas,plots,alleles,loci,grid,export")
    ap.add_argument("--only_trait", default=None); ap.add_argument("--only_panels", default=None)
    ap.add_argument("--strata", default=None, help="comma list of strata to run (default: all in the config)")
    ap.add_argument("--jobs", type=int, default=None, help="parallel grid variants (default config grid.jobs or 1)")
    ap.add_argument("--threads", type=int, default=None, help="total cores this run may use (default: all visible cores)")
    a = ap.parse_args()
    C = expand(json.load(open(a.config)))
    root = C["results_dir"]; os.makedirs(root, exist_ok=True)
    steps = a.steps.split(",")
    traits = {k: v for k, v in C["traits"].items() if a.only_trait in (None, k)}
    if "prep" in steps:
        for cmd in C["prep_commands"]:
            run(cmd)
    strata = C.get("strata", {"all": {}})
    if a.strata:
        strata = {k: v for k, v in strata.items() if k in a.strata.split(",")}
    raw_of = {T["column"]: T.get("raw", T["column"].replace("_rint", "")) for T in C["traits"].values()}
    cpu = a.threads or (len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count()) or 1
    failed = []
    for sname, filt in strata.items():
        out = os.path.join(root, sname) if len(C.get("strata", {})) > 0 else root
        os.makedirs(out, exist_ok=True)
        pheno = os.path.join(out, "stratum_pheno.tsv") if filt else C["pheno_table"]
        if filt:
            n = make_stratum(C["pheno_table"], filt, raw_of, pheno)
            print(f"== stratum {sname}: N = {n}", flush=True)
            if n < C.get("min_stratum_n", 1000):
                print(f"== stratum {sname} skipped: N < min_stratum_n ({C.get('min_stratum_n', 1000)})", flush=True)
                continue
        cnv = C["cnv_table"]
        try:
            for t, T in traits.items():
                num = ",".join(T["numeric"]); fac = ",".join(f for f in T["factors"] if not (f == "sex" and "sex" in (filt or {})))
                if "gwas" in steps:
                    g = C["gwas"]
                    args = ["cnv_gwas_qt.py", "--pheno", pheno, "--trait", T["column"], "--cnv", cnv, "--numeric", num, "--factors", fac,
                            "--chroms", g["chroms"], "--out", f"{out}/taskA_gwas/{t}", "--label", f"{C['cohort']} {sname} {t}",
                            "--min_ncnv", str(g["min_ncnv"]), "--freq_max", str(g["freq_max"])]
                    for opt in ["bp_dir", "known_loci", "segdups", "filtered_cnv_beds", "ref_probes"]:
                        if g.get(opt):
                            args += [f"--{opt}", g[opt].replace("{trait}", t)]
                    run(args)
                if "plots" in steps:
                    run(["plots.py", "--gwas_dir", f"{out}/taskA_gwas/{t}", "--title", f"{C['cohort']} {sname} {t} (RINT within sex)",
                         "--out", f"{out}/taskA_gwas/{t}/miami_{t}.png", "--min_ncnv_plot", str(C["gwas"]["min_ncnv"])])
                if "alleles" in steps:
                    run(["alleles.py", "--pheno", pheno, "--cnv", cnv, "--resources", C["alleles"]["resources"], "--traits", T["column"],
                         "--numeric", num, "--factors", fac, "--out", f"{out}/taskA_alleles/{t}"]
                        + (["--omar_matrix", C["alleles"]["omar_matrix"]] if C["alleles"].get("omar_matrix") and not filt else []))
                if "loci" in steps and C.get("loci"):
                    os.makedirs(f"{out}/locus_effects", exist_ok=True)
                    run(["locus_effects.py", "--pheno", pheno, "--cnv", cnv, "--loci", C["loci"]["file"], "--traits", T["column"],
                         "--numeric", num, "--factors", fac, "--min_carriers", str(C["loci"].get("min_carriers", 2)),
                         "--out", f"{out}/locus_effects/locus_effects_{t}.tsv"])
                if "burden" in steps and C.get("burden"):        # genome-wide signed burden, U-shape and monogenic-hit removal (dup_burden.py)
                    Bd = C["burden"]
                    run(["dup_burden.py", "--pheno", pheno, "--trait", T["column"], "--cnv", cnv, "--exons", Bd["exons"], "--gene_cov", Bd["gene_cov"],
                         "--clingen", Bd["clingen"], "--asym_loci", Bd["asym_loci"], "--numeric", num, "--factors", fac, "--out", f"{out}/burden/{t}",
                         "--min_n", str(Bd.get("min_n", 50)), "--single_def", Bd.get("single_def", "all_calls"),
                         "--rare_freq", str(Bd.get("rare_freq", 0.01))])
                if "axes" in steps:
                    X = C["axes"]
                    for pan in X["panels"]:
                        if a.only_panels and pan not in a.only_panels.split(","):
                            continue
                        run(["conga_axes_qt.py", "--pheno", pheno, "--trait", T["column"], "--cnv", cnv, "--twas", T["twas"], "--panels", pan,
                             "--exons", X["exons"], "--gene_cov", X["gene_cov"], "--loci", X["loci"], "--loo_runs", X["loo_runs"],
                             "--nrep", str(X["nrep"]), "--seed", str(X["seed"]), "--workers", str(X.get("workers", 1)), "--numeric", num,
                             "--factors", fac, "--out", f"{out}/taskB_axes/{t}_{pan}"])
                if "grid" in steps and C.get("grid") and sname in C["grid"].get("strata", list(strata)):
                    # Per design (2026-11-16): "panels", "loo_mode", "major_loci", "perms_file" and "outname" may override the grid defaults,
                    # so each design reproduces the UKBB run it is meta-analysed with. perms_file may use {trait} and {panels}.
                    G = C["grid"]; jobs = a.jobs or G.get("jobs", 1); thr = max(1, cpu // jobs)
                    calls = []; firsts = {}
                    for D in G["designs"]:
                        if t not in D.get("traits", [t]):
                            continue
                        pan = D.get("panels", G.get("panels", "gtex_nonbrain"))
                        if isinstance(pan, dict):
                            pan = pan.get(t, "gtex_nonbrain")
                        pf = D.get("perms_file", G.get("perms_file", "")).replace("{trait}", t).replace("{panels}", pan)
                        mode = D.get("loo_mode", G.get("loo_mode", "all")); maj = D.get("major_loci", G.get("major_loci", ""))
                        oname = D.get("outname", os.path.splitext(os.path.basename(D["file"]))[0])
                        vs = [v["name"] for v in json.load(open(os.path.join(HERE, D["file"])))]
                        vs = [v for v in vs if not D.get("only") or v in D["only"]]
                        for v in vs:
                            args = ["conga_grid_qt.py", "--pheno", pheno, "--trait", T["column"], "--cnv", cnv, "--twas", T["twas"],
                                    "--panels", pan, "--exons", G["exons"], "--gene_cov", G["gene_cov"],
                                    "--clingen", G["clingen"], "--loci", G["loci"], "--loo_runs", G["loo_runs"], "--design", D["file"], "--only", v,
                                    "--nrep", str(D.get("nrep", G.get("nrep", 1000))), "--seed", str(G.get("seed", 2026)),
                                    "--loo_mode", mode, "--numeric", num, "--factors", fac,
                                    "--out", f"{out}/grid/{oname}/{t}/{safe(v)}"]
                            if maj and mode != "all":
                                args += ["--major_loci", maj]
                            if pf:
                                args += ["--perms_file", pf]
                                if not os.path.exists(pf):
                                    if G.get("require_shared_shifts"):
                                        raise SystemExit(f"shared shift file missing: {pf} (copy the UKBB shift files into resources/shared_shifts; "
                                                         "a cohort's own shifts cannot be meta-analysed replicate by replicate)")
                                    if pf not in firsts:
                                        firsts[pf] = args; continue
                            calls.append(args)
                    # without require_shared_shifts: the first call per missing shift file creates it; the rest then load it
                    for pf, args in firsts.items():
                        run(args, threads=cpu)
                    with ThreadPoolExecutor(max_workers=jobs) as ex:
                        list(ex.map(lambda c: run(c, threads=thr), calls))
                F_ = C.get("figures", {})
                if "additivity" in steps and F_ and sname in F_.get("strata", ["all"]):
                    # Sacks 2026 Fig. 2I-J (constrained-gene dose response), Fig. 3H-I (asymmetric vs mirror duplications, sub-additivity),
                    # model figure (single hit vs constraint count). Same runs as UKBB additivity_2026-11-15, additivity_only{asym,mirror}_2026-11-15
                    # and additivity_hit3_2026-11-15 (the 'all' run = top5_hi with no exclusion = UKBB additivity_2026-11-15).
                    base = ["additivity_constraint.py", "--pheno", pheno, "--trait", T["column"], "--cnv", cnv, "--twas", T["twas"],
                            "--exons", F_["exons"], "--gene_cov", F_["gene_cov"], "--clingen", F_["clingen"], "--numeric", num, "--factors", fac]
                    if F_.get("phaplo"):
                        base += ["--phaplo", F_["phaplo"]]
                    ADD = [("all", [])]
                    # "only" runs: reference group = figures.only_ref (nocall = arrays/UKBB; no_rare_gene = WGS, Jon 2026-11-17), named
                    # onlyasym_<dos> as before; figures.only_ref_sensitivity adds the same runs with other references, suffixed _<ref>
                    refs = [(F_.get("only_ref", "nocall"), "")] + [(r, f"_{r}") for r in F_.get("only_ref_sensitivity", [])]
                    RX = ["--rare_freq", str(F_.get("rare_freq", 0.01))] + (["--segdup_bed", F_["segdups"]] if F_.get("segdups") else []) \
                         + (["--ref_exclude_loci", F_["sacks_loci"]] if F_.get("sacks_loci") else [])       # Jon 2026-11-18: no 34-interval carrier in any reference
                    for grp, lf in [("onlyasym", F_["asym_loci"]), ("onlymirror", F_["mirror_loci"])]:
                        for dos in ["both", "DUP"]:
                            for ref, sfx in refs:
                                ADD.append((f"{grp}_{dos}{sfx}", ["--exclude_loci", lf, "--exclude_dosage", dos, "--loci_mode", "only",
                                                                  "--only_ref", ref] + RX))
                    for h in (["top5_hi", "phaplo", "union"] if F_.get("phaplo") else ["top5_hi"]):
                        for ex_ in ["all", "noSacks34"]:
                            if h == "top5_hi" and ex_ == "all":
                                continue                                    # = the 'all' run
                            X_ = ["--hit_def", h] + (["--exclude_loci", F_["sacks_loci"], "--exclude_dosage", "both"] if ex_ == "noSacks34" else [])
                            ADD.append((f"hit_{h}_{ex_}", X_))
                    jobs = a.jobs or F_.get("jobs", 4); thr = max(1, cpu // jobs)
                    with ThreadPoolExecutor(max_workers=jobs) as ex:
                        list(ex.map(lambda r: run(base + r[1] + ["--out", f"{out}/additivity/{t}/{r[0]}"], threads=thr), ADD))
                if "scores" in steps and F_ and sname in F_.get("strata", ["all"]):
                    # Sacks 2026 Fig. 2L: 5-fold cross-fitted mirror/constraint scores and their additivity (UKBB score_axes_2026-11-15)
                    S_ = F_["scores"]
                    for v in S_["variants"].get(t, []):
                        run(["score_axes.py", "--pheno", pheno, "--trait", T["column"], "--cnv", cnv, "--twas", T["twas"],
                             "--panels", S_.get("panels", {}).get(t, "gtex_nonbrain"), "--exons", F_["exons"], "--gene_cov", F_["gene_cov"],
                             "--clingen", F_["clingen"], "--design", S_["design"], "--variant", v, "--folds", str(S_.get("folds", 5)),
                             "--seed", str(S_.get("seed", 2026)), "--numeric", num, "--factors", fac, "--out", f"{out}/scores/{t}/{v}"], threads=cpu)
                if "joint" in steps and C.get("joint") and sname in C["joint"].get("strata", ["all"]):
                    # Joint models of the gene-dosage axes behind the joint-model heatmap (UK Biobank columns): all CNVs and four carrier
                    # removals, each with 24 single-unit removals + minor_all + no_loci and matched nulls. Every run loads the UKBB shift file
                    # made from the same carrier subset (resources/shared_shifts/<trait>_<panels>_<subset>.npy); a missing file is an error.
                    J = C["joint"]; G = C["grid"]
                    sd = f"{out}/joint/strata"
                    if any(r["subset"] != "all" for r in J["runs"]) and not os.path.exists(f"{sd}/union_strata_info.json"):
                        run(["make_union_strata.py", "--pheno", pheno, "--cnv", cnv, "--exons", G["exons"], "--gene_cov", G["gene_cov"],
                             "--clingen", G["clingen"], "--asym_loci", J["asym_loci"], "--out", sd])
                    calls = []
                    for r in J["runs"]:
                        if r["trait"] != t:
                            continue
                        pf = J["perms_file"].format(trait=t, panels=r["panels"], subset=r["subset"])
                        if not os.path.exists(pf):
                            raise SystemExit(f"shared shift file missing: {pf}")
                        ph_ = pheno if r["subset"] == "all" else f"{sd}/pheno_{r['subset']}.tsv"
                        calls.append(["conga_grid_qt.py", "--pheno", ph_, "--trait", T["column"], "--cnv", cnv, "--twas", T["twas"], "--panels", r["panels"],
                                      "--exons", G["exons"], "--gene_cov", G["gene_cov"], "--clingen", G["clingen"], "--loci", G["loci"],
                                      "--loo_runs", G["loo_runs"], "--design", r["design"], "--only", r["variant"], "--nrep", str(J.get("nrep", 2000)),
                                      "--seed", str(G.get("seed", 2026)), "--loo_mode", "all_minor", "--major_loci", J["major_loci"],
                                      "--perms_file", pf, "--numeric", num, "--factors", fac, "--out", f"{out}/joint/{r['name']}/{t}/{safe(r['variant'])}"])
                    jobs = a.jobs or J.get("jobs", 5); thr = max(1, cpu // jobs)
                    with ThreadPoolExecutor(max_workers=jobs) as ex:
                        list(ex.map(lambda c: run(c, threads=thr), calls))
        except subprocess.CalledProcessError as e:
            # one stratum's failure (e.g. a model that cannot be fitted in a small ancestry stratum) skips the rest of that stratum
            # only; the other strata still run, and the pipeline exits non-zero at the end listing the failed strata
            failed.append(sname); print(f"== stratum {sname} FAILED ({' '.join(map(str, e.cmd[1:3]))}): rest of this stratum skipped", flush=True)

    if failed:
        raise SystemExit(f"failed strata (rerun them with --strata {','.join(failed)}): {failed}; export not run")
    if "export" in steps:
        run(["export_aggregate.py", "--results", root, "--export", C["export_dir"], "--min_cell", str(C["min_cell"])]
            + (["--strata", ",".join(C["strata"])] if C.get("strata") else [])
            + (["--partitions", json.dumps(C["partitions"])] if C.get("partitions") else []))


if __name__ == "__main__":
    main()
