"""Single entry point. One code base, one JSON config per cohort (configs/ukbb_expanse.json, configs/aou_workbench.json).

  python run_pipeline.py --config configs/<cohort>.json --steps prep,gwas,alleles,axes,plots,export [--only_trait height] [--only_panels gtex_nonbrain]

The cluster (SLURM on Expanse) or the Workbench (Jupyter cell / dsub task) only decides *where* this runs; nothing in
the steps knows about either. Paths in the config may use ${ENV} variables (e.g. ${WORKSPACE_BUCKET}, ${HOME}).
"""
import argparse
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def expand(o):
    if isinstance(o, str):
        return os.path.expandvars(o)
    if isinstance(o, dict):
        return {k: expand(v) for k, v in o.items()}
    if isinstance(o, list):
        return [expand(v) for v in o]
    return o


def run(args):
    cmd = [sys.executable, *args]
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, cwd=HERE)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True); ap.add_argument("--steps", default="prep,gwas,alleles,axes,plots,export")
    ap.add_argument("--only_trait", default=None); ap.add_argument("--only_panels", default=None)
    a = ap.parse_args()
    C = expand(json.load(open(a.config)))
    out = C["results_dir"]; os.makedirs(out, exist_ok=True)
    steps = a.steps.split(",")
    traits = {k: v for k, v in C["traits"].items() if a.only_trait in (None, k)}
    pheno, cnv = C["pheno_table"], C["cnv_table"]
    if "prep" in steps:
        for cmd in C["prep_commands"]:                      # cohort adapters (UKBB: prep_ukbb.py; AoU: aou_sv_to_cnv.py + aou_pheno_extract.py)
            run(cmd)
    for t, T in traits.items():
        num = ",".join(T["numeric"]); fac = ",".join(T["factors"])
        if "gwas" in steps:
            g = C["gwas"]
            args = ["cnv_gwas_qt.py", "--pheno", pheno, "--trait", T["column"], "--cnv", cnv, "--numeric", num, "--factors", fac,
                    "--chroms", g["chroms"], "--out", f"{out}/taskA_gwas/{t}", "--label", f"{C['cohort']} {t}",
                    "--min_ncnv", str(g["min_ncnv"]), "--freq_max", str(g["freq_max"])]
            for opt in ["bp_dir", "known_loci", "segdups", "filtered_cnv_beds"]:
                if g.get(opt):
                    args += [f"--{opt}", g[opt]]
            run(args)
        if "plots" in steps:
            run(["plots.py", "--gwas_dir", f"{out}/taskA_gwas/{t}", "--title", f"{C['cohort']} {t} (RINT within sex)",
                 "--out", f"{out}/taskA_gwas/{t}/miami_{t}.png", "--min_ncnv_plot", str(C["gwas"]["min_ncnv"])])
        if "alleles" in steps:
            run(["alleles.py", "--pheno", pheno, "--cnv", cnv, "--resources", C["alleles"]["resources"], "--traits", T["column"],
                 "--numeric", num, "--factors", fac, "--out", f"{out}/taskA_alleles/{t}"]
                + (["--omar_matrix", C["alleles"]["omar_matrix"]] if C["alleles"].get("omar_matrix") else []))
        if "axes" in steps:
            X = C["axes"]
            for pan in X["panels"]:
                if a.only_panels and pan not in a.only_panels.split(","):
                    continue
                run(["conga_axes_qt.py", "--pheno", pheno, "--trait", T["column"], "--cnv", cnv, "--twas", T["twas"], "--panels", pan,
                     "--exons", X["exons"], "--gene_cov", X["gene_cov"], "--loci", X["loci"], "--loo_runs", X["loo_runs"],
                     "--nrep", str(X["nrep"]), "--seed", str(X["seed"]), "--workers", str(X.get("workers", 1)), "--numeric", num, "--factors", fac,
                     "--out", f"{out}/taskB_axes/{t}_{pan}"])
    if "export" in steps:
        run(["export_aggregate.py", "--results", out, "--export", C["export_dir"], "--min_cell", str(C["min_cell"])])


if __name__ == "__main__":
    main()
