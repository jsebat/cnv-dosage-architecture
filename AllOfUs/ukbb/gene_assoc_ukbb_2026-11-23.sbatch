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
#SBATCH -J gene_assoc_ukbb
#SBATCH -o /expanse/projects/sebat1/jsebat/sacks2026_aou/fig3_2026-11-23/logs/%x.%j.log
# Gene-based CNV association (Jon 2026-11-23: "add a gene based test so Milad can output gene based sum stats"), UK Biobank, through the same
# entry point Milad uses: run_pipeline.py --steps genes,genes_export. Per gene and dosage: trait ~ covariates + carrier indicator (DEL = any
# exon overlap >= 1 bp, DUP = full-length transcript inside the call), height and BMI, strata all / EUR / AFR / ASN / MIX / female / male.
# Then: UKBB export (min_cell 0: rows with >= 5 carriers), a second export under the All of Us rule (min_cell 20, with the cross-stratum
# rule) to exercise the suppression and its audit, a one-cohort meta, a two-copy meta self-test (code check, not a result), and a check of the
# gene z against a direct per-gene regression for 20 genes.
# null model: none (analytic OLS). Unique output dir: $B/gene_assoc_ukbb_$SLURM_JOB_ID (code snapshot inside).
set -e
B=/expanse/projects/sebat1/jsebat/sacks2026_aou/fig3_2026-11-23
PY=${PY:-/home/jsebat/micromamba/envs/dnmt3a-py/bin/python}; REPO=${REPO:-$B/code/AllOfUs}; export PATH=$(dirname $PY):$PATH
W=$B/gene_assoc_ukbb_$SLURM_JOB_ID; mkdir -p $W; S=$W/code_snapshot; mkdir -p $S; cp -r $REPO/cnv_anthro $REPO/meta $REPO/ukbb $S/; cd $S/cnv_anthro
$PY - <<EOF
import json
C = json.load(open("configs/ukbb_expanse.json"))
C["results_dir"] = "$W/results"; C["export_dir"] = "$W/export"; C["min_cell"] = 0; C["cohort"] = "UKBB"
C["partitions"] = {"sex": ["female", "male"], "ancestry": ["EUR", "AFR", "ASN", "MIX"]}
json.dump(C, open("$W/config.json", "w"), indent=1)
C["min_cell"] = 20; C["genes"]["export_dir"] = "$W/export_genes_mincell20"; json.dump(C, open("$W/config_mincell20.json", "w"), indent=1)
EOF
$PY run_pipeline.py --config $W/config.json --steps genes,genes_export
$PY run_pipeline.py --config $W/config_mincell20.json --steps genes_export --strata all
(cd $W && tar czf ukbb_gene_assoc_export.tar.gz export_genes && sha256sum ukbb_gene_assoc_export.tar.gz > ukbb_gene_assoc_export.tar.gz.sha256 && ls -la ukbb_gene_assoc_export.tar.gz)
$PY ../meta/meta_gene_assoc.py --cohort UKBB=$W/export_genes/all/gene_assoc --out $W/meta_one
$PY ../meta/meta_gene_assoc.py --cohort UKBB=$W/export_genes/all/gene_assoc --cohort UKBBcopy20=$W/export_genes_mincell20/all/gene_assoc --out $W/meta_two_selftest
$PY ../ukbb/check_gene_assoc.py --config $W/config.json --table $W/results/all/gene_assoc/height/gene_assoc.tsv.gz --trait height --meta_one $W/meta_one --meta_two $W/meta_two_selftest > $W/checks.txt 2>&1
echo "checks exit $?"; cat $W/checks.txt
chmod -R g+rX $W
