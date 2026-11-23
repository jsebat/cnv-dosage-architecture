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
#SBATCH -J fig3_callset_test
#SBATCH -o /expanse/projects/sebat1/jsebat/sacks2026_aou/fig3_2026-11-23/logs/%x.%j.log
# Test of the call-set step and the gene-based tables on UK Biobank, through Milad's entry points (Jon 2026-11-23, final decision: Fig. 3
# uses the PGC freeze filters in every cohort; the per-position frequency is recorded, no position filter). Runs the recording step, the two
# cheapest Fig. 3 pieces (totdose fit, dose bins) and the gene tables (stratum 'all'), then: (1) totdose and bins must EQUAL the finished
# UK Biobank export REF (same calls, so no UK Biobank rerun is needed); (2) the new export must meta-analyse with REF (the call-filter guard
# accepts PGC-filtered cohorts); (3) gene z = direct regression; exports under min_cell 0 and 20.
# null model: totdose: genome shift, 2,000 shared shifts (its total-genes term is analytic); bins and gene tests: analytic OLS.
# Unique output dir: $B/callset_gene_test_$SLURM_JOB_ID (code snapshot inside).
#   sbatch --export=ALL,REF=<ukbb_export_<id>/export/all/fig3> fig3_filter_gene_test_2026-11-23.sbatch
set -e
B=/expanse/projects/sebat1/jsebat/sacks2026_aou/fig3_2026-11-23; REF=${REF:?set REF}
PY=${PY:-/home/jsebat/micromamba/envs/dnmt3a-py/bin/python}; REPO=${REPO:-$B/code/AllOfUs}; export PATH=$(dirname $PY):$PATH
W=$B/callset_gene_test_$SLURM_JOB_ID; mkdir -p $W; S=$W/code_snapshot; mkdir -p $S; cp -r $REPO/cnv_anthro $REPO/meta $REPO/ukbb $S/; cd $S/cnv_anthro
$PY - <<EOF
import json
C = json.load(open("configs/ukbb_expanse.json"))
C["results_dir"] = "$W/results"; C["export_dir"] = "$W/export"; C["min_cell"] = 0; C["cohort"] = "UKBB"
json.dump(C, open("$W/config.json", "w"), indent=1)
C["min_cell"] = 20; C["fig3"]["export_dir"] = "$W/export_fig3_mincell20"; C["genes"]["export_dir"] = "$W/export_genes_mincell20"
json.dump(C, open("$W/config_mincell20.json", "w"), indent=1)
EOF
$PY run_fig3.py --config $W/config.json --steps filter,fit,export --only_fit totdose,dose_bins --jobs 4
$PY run_fig3.py --config $W/config_mincell20.json --steps export --only_fit totdose,dose_bins
$PY run_pipeline.py --config $W/config.json --steps genes,genes_export --strata all
$PY run_pipeline.py --config $W/config_mincell20.json --steps genes_export --strata all
tail -7 $W/results/all/fig3/logs/cnv_filter.log
$PY ../meta/meta_fig3.py --cohort UKBB=$W/export_fig3/all/fig3 --out $W/meta_fig3_one --code $S/cnv_anthro --parts joint,dose_bins
$PY ../meta/meta_fig3.py --cohort UKBB=$REF --cohort UKBBcopy=$W/export_fig3/all/fig3 --out $W/meta_with_ref_selftest --code $S/cnv_anthro --parts joint,dose_bins
echo "call-filter guard: PGC-filtered cohorts combine (exit $?)"
$PY ../meta/meta_fig3.py --cohort UKBB=$REF --out $W/meta_ref_one --code $S/cnv_anthro --parts joint,dose_bins
set +e
$PY ../meta/compare_fig3_tables.py --a totdose=$W/meta_ref_one/totdose,dose_bins=$W/meta_ref_one/dose_bins --b totdose=$W/meta_fig3_one/totdose,dose_bins=$W/meta_fig3_one/dose_bins > $W/compare_with_ref.txt 2>&1
echo "totdose + bins equal the finished UK Biobank export: exit $?"; cut -c1-200 $W/compare_with_ref.txt
$PY ../meta/meta_gene_assoc.py --cohort UKBB=$W/export_genes/all/gene_assoc --out $W/meta_gene_assoc_one
$PY ../meta/meta_gene_assoc.py --cohort UKBB=$W/export_genes/all/gene_assoc --cohort UKBBcopy20=$W/export_genes_mincell20/all/gene_assoc --out $W/meta_gene_assoc_two_selftest
$PY ../ukbb/check_gene_assoc.py --config $W/config.json --table $W/results/all/gene_assoc/height/gene_assoc.tsv.gz --trait height --meta_one $W/meta_gene_assoc_one --meta_two $W/meta_gene_assoc_two_selftest > $W/checks_gene_assoc.txt 2>&1
echo "gene checks exit $?"; cat $W/checks_gene_assoc.txt
$PY - <<EOF
import pandas as pd, json, os
g = pd.read_csv("$W/export_genes/all/gene_assoc/height/gene_assoc.tsv.gz", sep="\t"); print("gene table columns:", list(g.columns)); print("max_pos_freq:", g.max_pos_freq.describe().round(5).to_dict())
print("gene sets exported:", sorted(os.listdir("$W/export_genes/all"))); print(json.load(open("$W/export_fig3/all/fig3/FIG3_EXPORT_INFO.json"))["cnv_filter"])
print(sorted(os.listdir("$W/export_fig3/all/fig3/cnv_filter")), sorted(os.listdir("$W/results/all/fig3/cnv_filter")))
EOF
chmod -R g+rX $W
