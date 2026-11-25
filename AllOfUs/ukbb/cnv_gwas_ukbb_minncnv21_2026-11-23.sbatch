#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH -t 04:00:00
#SBATCH --nice=100
#SBATCH -J cnv_gwas_ukbb21
#SBATCH -o /expanse/projects/sebat1/jsebat/sacks2026_aou/fig3_2026-11-23/logs/%x.%j.log
# UK Biobank probe-level CNV GWAS with the carrier rule of All of Us: a probe needs at least 21 carriers (Jon 2026-11-23: the 21 rule is the
# GWAS rule, in every cohort; UK Biobank used 12, Omar's PGC value). Same code and inputs as the 2026-11-15 run (pipeline_ukbb_strata), only
# gwas.min_ncnv 12 -> 21; every stratum, height and BMI; probe GWAS + Miami plots. The 2026-11-15 results are not touched.
# Then a table of probes kept, independent tests, genome-wide threshold and significant groups, old (12) vs new (21).
# null model: none (OLS t-test per probe; genome-wide threshold = 0.05 / independent tests, Omar's greedy Jaccard merge of adjacent probes).
# Unique output dir: $B/cnv_gwas_ukbb_minncnv21_$SLURM_JOB_ID (code snapshot inside).
set -e
B=/expanse/projects/sebat1/jsebat/sacks2026_aou/fig3_2026-11-23; OLD=/expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/results/pipeline_ukbb_strata
PY=${PY:-/home/jsebat/micromamba/envs/dnmt3a-py/bin/python}; REPO=${REPO:-$B/code/AllOfUs}; export PATH=$(dirname $PY):$PATH
W=$B/cnv_gwas_ukbb_minncnv21_$SLURM_JOB_ID; mkdir -p $W; S=$W/code_snapshot; mkdir -p $S; cp -r $REPO/cnv_anthro $S/; cd $S/cnv_anthro
$PY - <<EOF
import json
C = json.load(open("configs/ukbb_expanse.json")); assert C["gwas"]["min_ncnv"] == 21
C["results_dir"] = "$W/results"; C["export_dir"] = "$W/export"; C["cohort"] = "UKBB"
json.dump(C, open("$W/config.json", "w"), indent=1)
EOF
$PY run_pipeline.py --config $W/config.json --steps gwas,plots
$PY - <<EOF
import json, glob, os
print("stratum trait | N | DEL kept 12 -> 21 | DEL indep tests | DEL GWS groups | DUP kept | DUP indep tests | DUP GWS groups")
for f in sorted(glob.glob("$W/results/*/taskA_gwas/*/run_info.json")):
    st, t = f.split("/")[-4], f.split("/")[-2]; n = json.load(open(f)); o = "$OLD/%s/taskA_gwas/%s/run_info.json" % (st, t)
    o = json.load(open(o)) if os.path.exists(o) else {}
    g = lambda d, k: d.get(k, "NA")
    print(st, t, "|", n["N"], "|", " | ".join(f"{g(o, k)} -> {g(n, k)}" for k in ("DEL_probes_kept", "DEL_independent_tests", "DEL_gws_groups", "DUP_probes_kept", "DUP_independent_tests", "DUP_gws_groups")))
EOF
chmod -R g+rX $W
