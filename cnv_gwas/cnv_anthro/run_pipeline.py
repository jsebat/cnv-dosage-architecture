"""Single entry point. One code base, one JSON config per cohort (configs/ukbb_expanse.json, configs/aou_workbench.json).

  python run_pipeline.py --config configs/<cohort>.json --steps prep,gwas,plots,alleles,loci,axes,grid,export
        [--only_trait height] [--only_panels gtex_nonbrain] [--strata all,eur,female] [--jobs 8]

Steps
  prep     cohort adapters (UKBB: prep_ukbb.py; AoU: aou_sv_to_cnv.py + aou_pheno_extract.py)
  gwas     probe-level CNV GWAS (cnv_gwas_qt.py); with gwas.ref_probes also at the reference cohort's probe positions (meta-analysis)
  plots    Miami plot per trait
  alleles  Table S9 allele tests (alleles.py)
  loci     per-locus DEL/DUP effects at the growth loci (locus_effects.py; Sacks recurrent intervals + UKBB-discovered loci)
  axes     ConGA consensus2 axes (conga_axes_qt.py)
  grid     genome-wide univariate grid, call-size split and joint models (conga_grid_qt.py), shared shifts for meta-analysis
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
    for sname, filt in strata.items():
        out = os.path.join(root, sname) if len(C.get("strata", {})) > 0 else root
        os.makedirs(out, exist_ok=True)
        pheno = os.path.join(out, "stratum_pheno.tsv") if filt else C["pheno_table"]
        if filt:
            n = make_stratum(C["pheno_table"], filt, raw_of, pheno)
            print(f"== stratum {sname}: N = {n}", flush=True)
        cnv = C["cnv_table"]
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
                G = C["grid"]; jobs = a.jobs or G.get("jobs", 1); thr = max(1, cpu // jobs)
                calls = []
                for D in G["designs"]:
                    if t not in D.get("traits", [t]):
                        continue
                    vs = [v["name"] for v in json.load(open(os.path.join(HERE, D["file"])))]
                    vs = [v for v in vs if not D.get("only") or v in D["only"]]
                    for v in vs:
                        args = ["conga_grid_qt.py", "--pheno", pheno, "--trait", T["column"], "--cnv", cnv, "--twas", T["twas"],
                                "--panels", G.get("panels", "gtex_nonbrain"), "--exons", G["exons"], "--gene_cov", G["gene_cov"],
                                "--clingen", G["clingen"], "--loci", G["loci"], "--loo_runs", G["loo_runs"], "--design", D["file"], "--only", v,
                                "--nrep", str(D.get("nrep", G.get("nrep", 1000))), "--seed", str(G.get("seed", 2026)),
                                "--loo_mode", G.get("loo_mode", "all"), "--numeric", num, "--factors", fac,
                                "--out", f"{out}/grid/{os.path.splitext(os.path.basename(D['file']))[0]}/{t}/{safe(v)}"]
                        if G.get("major_loci"):
                            args += ["--major_loci", G["major_loci"]]
                        if G.get("perms_file"):
                            args += ["--perms_file", G["perms_file"].replace("{trait}", t)]
                        calls.append(args)
                # the first call creates the shared shift file if it does not exist yet; the rest then load it
                if calls and G.get("perms_file") and not os.path.exists(G["perms_file"].replace("{trait}", t)):
                    run(calls.pop(0), threads=cpu)
                with ThreadPoolExecutor(max_workers=jobs) as ex:
                    list(ex.map(lambda c: run(c, threads=thr), calls))
    if "export" in steps:
        run(["export_aggregate.py", "--results", root, "--export", C["export_dir"], "--min_cell", str(C["min_cell"])]
            + (["--strata", ",".join(C["strata"])] if C.get("strata") else [])
            + (["--partitions", json.dumps(C["partitions"])] if C.get("partitions") else []))


if __name__ == "__main__":
    main()
