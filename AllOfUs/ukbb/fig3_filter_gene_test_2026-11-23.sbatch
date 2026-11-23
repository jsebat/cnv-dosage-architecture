#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH -t 02:00:00
#SBATCH --nice=100
#SBATCH -J fig3_filter_test
#SBATCH -o /expanse/projects/sebat1/jsebat/sacks2026_aou/fig3_2026-11-23/logs/%x.%j.log
# Test of the rare-CNV filter step and of the gene-based test with it, on UK Biobank, through Milad's entry points (Jon 2026-11-23). NOT the
# UK Biobank Fig. 3 rerun (that needs Jon's go): only the filter, the two cheapest Fig. 3 pieces on the filtered calls (totdose fit, dose
# bins) and the gene-based tables (every call; rare calls), stratum 'all', then their exports, a one-cohort meta and the guard that refuses
# to meta-analyse a filtered with an unfiltered cohort (UNF = the unfiltered UK Biobank Fig. 3 export).
# null model: totdose: genome shift, 2,000 shared shifts (its total-genes term is analytic); bins and gene tests: analytic OLS.
# Unique output dir: $B/filter_gene_test_$SLURM_JOB_ID (code snapshot inside).
#   sbatch --export=ALL,UNF=<ukbb_export_<id>/export/all/fig3> fig3_filter_gene_test_2026-11-23.sbatch
set -e
B=/expanse/projects/sebat1/jsebat/sacks2026_aou/fig3_2026-11-23; UNF=${UNF:?set UNF}
PY=${PY:-/home/jsebat/micromamba/envs/dnmt3a-py/bin/python}; REPO=${REPO:-$B/code/AllOfUs}; export PATH=$(dirname $PY):$PATH
W=$B/filter_gene_test_$SLURM_JOB_ID; mkdir -p $W; S=$W/code_snapshot; mkdir -p $S; cp -r $REPO/cnv_anthro $REPO/meta $REPO/ukbb $S/; cd $S/cnv_anthro
$PY - <<EOF
import json
C = json.load(open("configs/ukbb_expanse.json"))
C["results_dir"] = "$W/results"; C["export_dir"] = "$W/export"; C["min_cell"] = 0; C["cohort"] = "UKBB"
json.dump(C, open("$W/config.json", "w"), indent=1)
C["min_cell"] = 20; C["fig3"]["export_dir"] = "$W/export_fig3_mincell20"; C["genes"]["export_dir"] = "$W/export_genes_mincell20"
json.dump(C, open("$W/config_mincell20.json", "w"), indent=1)
EOF
$PY run_fig3.py --config $W/config.json --steps fit,export --only_fit totdose,dose_bins --jobs 4
$PY run_fig3.py --config $W/config_mincell20.json --steps export --only_fit totdose,dose_bins
$PY run_pipeline.py --config $W/config.json --steps genes,genes_export --strata all
$PY run_pipeline.py --config $W/config_mincell20.json --steps genes_export --strata all
cat $W/results/all/fig3/logs/cnv_filter.log | tail -8
$PY ../meta/meta_fig3.py --cohort UKBB=$W/export_fig3/all/fig3 --out $W/meta_fig3_one --code $S/cnv_anthro --parts joint,dose_bins
set +e
$PY ../meta/meta_fig3.py --cohort UKBB=$UNF --cohort UKBBrare=$W/export_fig3/all/fig3 --out $W/meta_should_fail --code $S/cnv_anthro --parts joint > $W/guard.txt 2>&1
echo "guard (must be non-zero): exit $?"; tail -2 $W/guard.txt | cut -c1-300
set -e
for s in gene_assoc gene_assoc_rare; do
  $PY ../meta/meta_gene_assoc.py --cohort UKBB=$W/export_genes/all/$s --out $W/meta_${s}_one
  $PY ../meta/meta_gene_assoc.py --cohort UKBB=$W/export_genes/all/$s --cohort UKBBcopy20=$W/export_genes_mincell20/all/$s --out $W/meta_${s}_two_selftest
done
set +e
$PY ../ukbb/check_gene_assoc.py --config $W/config.json --table $W/results/all/gene_assoc/height/gene_assoc.tsv.gz --trait height --meta_one $W/meta_gene_assoc_one --meta_two $W/meta_gene_assoc_two_selftest > $W/checks_gene_assoc.txt 2>&1
echo "gene checks (all calls) exit $?"; cat $W/checks_gene_assoc.txt
$PY ../ukbb/check_gene_assoc.py --config $W/config.json --cnv $W/results/all/fig3/cnv_filter/cnv_rare.cnv --table $W/results/all/gene_assoc_rare/height/gene_assoc.tsv.gz --trait height --meta_one $W/meta_gene_assoc_rare_one --meta_two $W/meta_gene_assoc_rare_two_selftest > $W/checks_gene_assoc_rare.txt 2>&1
echo "gene checks (rare calls) exit $?"; cat $W/checks_gene_assoc_rare.txt
$PY - <<EOF
import pandas as pd, json
a = pd.read_csv("$UNF/dose_bins/height/fig3_dose_bins.tsv", sep="\t"); b = pd.read_csv("$W/export_fig3/all/fig3/dose_bins/height/fig3_dose_bins.tsv", sep="\t")
k = ["run", "partition", "size", "cls", "side", "bin"]; m = a.merge(b, on=k, suffixes=("_all", "_rare")); m = m[(m.partition == "uni") & (m.run == "full") & (m["size"] == "N1_inf") & (m.cls == "all")]
print("height, total genes, all CNVs: bins with every call vs rare calls only"); print(m[["side", "bin", "n_all", "beta_all", "SE_all", "n_rare", "beta_rare", "SE_rare"]].to_string(index=False))
g = pd.read_csv("$W/export_genes/all/gene_assoc_rare/height/gene_assoc.tsv.gz", sep="\t"); print("gene table columns:", list(g.columns)); print(g.max_pos_freq.describe().round(5).to_dict())
print(json.load(open("$W/export_fig3/all/fig3/FIG3_EXPORT_INFO.json"))["cnv_filter"])
EOF
chmod -R g+rX $W; chmod o-rwx $W/results/all/fig3/cnv_filter/cnv_rare.cnv
