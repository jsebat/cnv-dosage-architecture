#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH -J fig3_selftest
#SBATCH -o /expanse/projects/sebat1/jsebat/sacks2026_aou/fig3_2026-11-23/logs/%x.%j.log
# Self-test of the multi-cohort path of meta_fig3.py (code check only, NOT a result): UK Biobank fits made with the current code are exported
# twice, without suppression and under the All of Us small-cell rule (min_cell 20: terms, genes and bins with 1-20 carriers withheld), and the
# two copies are meta-analysed as if they were two cohorts. Exercises the export suppression and its audit, the cohort checks, the sqrt(N)
# Stouffer with withheld cells, the gene mirror with withheld genes and the inverse-variance bins. Expected: statistics about sqrt(2) x the
# UK Biobank value where both copies report them.
# null model: as meta_fig3.py (shared genome shifts; gene-circle shift); nothing is refitted.
#   sbatch --export=ALL,J=<aou_path_test_<id>/results/all/fig3>,U=<ukbb_export_<id>/results/all/fig3>,W=<ukbb_export_<id>> fig3_two_cohort_selftest_2026-11-23.sbatch
# J: joint fits + gene mirror refitted with the current code (they carry n_carriers_term); U: univariate fits and dose bins. Output: $B/selftest_two_cohort_$SLURM_JOB_ID
set -e
B=/expanse/projects/sebat1/jsebat/sacks2026_aou/fig3_2026-11-23; J=${J:?set J}; U=${U:?set U}; W=${W:?set W=<UKBB export work dir (config, code snapshot)>}
PY=${PY:-/home/jsebat/micromamba/envs/dnmt3a-py/bin/python}; O=$B/selftest_two_cohort_$SLURM_JOB_ID; R=$O/results/all/fig3; mkdir -p $R
for d in size2t_mirror all2t_mirror totdose gene_mirror; do ln -s $J/$d $R/$d; done
for d in uni_N1_inf uni_N1_1 uni_N2_inf dose_bins; do ln -s $U/$d $R/$d; done
cd $W/code_snapshot/cnv_anthro; export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
$PY - <<EOF
import json
C = json.load(open("$W/ukbb_fig3_config.json")); C["results_dir"] = "$O/results"
for k, name in ((0, "A"), (20, "B_mincell20")):
    C["min_cell"] = k; C["fig3"]["export_dir"] = "$O/export_" + name; C["cohort"] = "UKBB_copy_" + name
    json.dump(C, open("$O/config_" + name + ".json", "w"), indent=1)
EOF
$PY run_fig3.py --config $O/config_A.json --steps export
$PY run_fig3.py --config $O/config_B_mincell20.json --steps export
$PY ../meta/meta_fig3.py --cohort UKBB=$O/export_A/all/fig3 --out $O/meta_one --code $W/code_snapshot/cnv_anthro --parts joint,dose_bins
$PY ../meta/meta_fig3.py --cohort UKBB=$O/export_A/all/fig3 --cohort UKBBcopy20=$O/export_B_mincell20/all/fig3 --out $O/meta_two --code $W/code_snapshot/cnv_anthro
cat $O/meta_two/meta_checks.tsv
$PY - <<EOF
import pandas as pd, numpy as np
k = ["trait", "variant", "size", "stratum", "stat"]
a = pd.read_csv("$O/meta_one/size2t_mirror/mirror4_cross_summary.tsv", sep="\t"); b = pd.read_csv("$O/meta_two/size2t_mirror/mirror4_cross_summary.tsv", sep="\t")
m = a.merge(b, on=k, suffixes=("_one", "_two")); m = m[m.stat == "M4"]
print("M4 obs ratio two-copy / one (expect ~1.414 where no cell is withheld):", (m.obs_two / m.obs_one).round(3).describe()[["min", "50%", "max"]].to_dict())
kb = ["trait", "run", "partition", "size", "cls", "side", "bin"]
a = pd.read_csv("$O/meta_one/dose_bins/height/fig3_dose_bins.tsv", sep="\t"); b = pd.read_csv("$O/meta_two/dose_bins/height/fig3_dose_bins.tsv", sep="\t")
m = a.merge(b, on=kb, suffixes=("_one", "_two")); m = m[(m.side != "none") & np.isfinite(pd.to_numeric(m.SE_one, errors="coerce"))]
print("bins: beta ratio (expect 1):", (m.beta_two / m.beta_one).round(4).describe()[["min", "max"]].to_dict(), "| SE ratio (expect 0.707):", (m.SE_two / m.SE_one).round(4).describe()[["min", "max"]].to_dict())
g = pd.read_csv("$O/meta_two/gene_mirror/height/gene_mirror_summary.tsv", sep="\t"); print(g[(g.stratum == "all") & (g.stat == "total")][["size", "obs", "z", "P", "n_genes", "loo_held", "n_loo"]].to_string(index=False))
EOF
chmod -R g+rX $O
