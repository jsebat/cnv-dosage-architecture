#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH -t 00:45:00
#SBATCH -J gene_export21
#SBATCH -o /expanse/projects/sebat1/jsebat/sacks2026_aou/fig3_2026-11-23/logs/%x.%j.log
# Gene-based CNV association, UK Biobank: export with the floor of 21 carriers that every cohort uses (Jon 2026-11-23: "ok 21"), from the
# finished in-cohort tables (SRC/results: every gene with >= 1 carrier; nothing is refitted), then the one-cohort meta, a two-copy self-test
# and the checks (gene z vs a direct regression; genes kept, lambda and genes at P < 2.5e-6 per trait and dosage).
# null model: none (analytic OLS). Unique output dir: $B/gene_export_ukbb_floor21_$SLURM_JOB_ID (code snapshot inside).
#   sbatch --export=ALL,SRC=<gene_assoc_ukbb_<id>> gene_export_ukbb_floor21_2026-11-23.sbatch
set -e
B=/expanse/projects/sebat1/jsebat/sacks2026_aou/fig3_2026-11-23; SRC=${SRC:?set SRC=<gene_assoc_ukbb_<jobid>>}
PY=${PY:-/home/jsebat/micromamba/envs/dnmt3a-py/bin/python}; REPO=${REPO:-$B/code/AllOfUs}; export PATH=$(dirname $PY):$PATH
W=$B/gene_export_ukbb_floor21_$SLURM_JOB_ID; mkdir -p $W; S=$W/code_snapshot; mkdir -p $S; cp -r $REPO/cnv_anthro $REPO/meta $REPO/ukbb $S/; cd $S/cnv_anthro
$PY - <<EOF
import json
C = json.load(open("configs/ukbb_expanse.json")); assert C["genes"]["min_carriers"] == 21
C["results_dir"] = "$SRC/results"; C["export_dir"] = "$W/export"; C["min_cell"] = 0; C["cohort"] = "UKBB"
C["partitions"] = {"sex": ["female", "male"], "ancestry": ["EUR", "AFR", "ASN", "MIX"]}
json.dump(C, open("$W/config.json", "w"), indent=1)
C["min_cell"] = 20; C["genes"]["export_dir"] = "$W/export_genes_mincell20"; json.dump(C, open("$W/config_mincell20.json", "w"), indent=1)
EOF
$PY run_pipeline.py --config $W/config.json --steps genes_export --strata all
$PY run_pipeline.py --config $W/config_mincell20.json --steps genes_export --strata all
(cd $W && tar czf ukbb_gene_assoc_export.tar.gz export_genes && sha256sum ukbb_gene_assoc_export.tar.gz > ukbb_gene_assoc_export.tar.gz.sha256 && ls -la ukbb_gene_assoc_export.tar.gz)
$PY ../meta/meta_gene_assoc.py --cohort UKBB=$W/export_genes/all/gene_assoc --out $W/meta_one
$PY ../meta/meta_gene_assoc.py --cohort UKBB=$W/export_genes/all/gene_assoc --cohort UKBBcopy20=$W/export_genes_mincell20/all/gene_assoc --out $W/meta_two_selftest
set +e
for t in height bmi; do
  $PY ../ukbb/check_gene_assoc.py --config $W/config.json --table $SRC/results/all/gene_assoc/$t/gene_assoc.tsv.gz --trait $t --meta_one $W/meta_one --meta_two $W/meta_two_selftest > $W/checks_$t.txt 2>&1
  echo "checks $t exit $?"; cat $W/checks_$t.txt
done
chmod -R g+rX $W
