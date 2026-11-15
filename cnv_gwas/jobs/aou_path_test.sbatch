#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH -t 02:00:00
#SBATCH --nice=0
#SBATCH --job-name=anthro_aoutest
#SBATCH -o /expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/logs/aoutest_%j.out
#SBATCH -e /expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/logs/aoutest_%j.err
# AoU code path, end to end, on UKBB data (config swap only).
set -euo pipefail
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
PY=/home/jsebat/micromamba/envs/dnmt3a-py/bin/python
P=/expanse/projects/sebat1/jsebat/ukbb_anthro_cnv; C=$P/jobs/snap_aoutest/cnv_anthro
T=$P/data/aou_path_test; mkdir -p $T $P/results/aou_path_test
cd $C
$PY tests/build_ukbb_as_aou.py --pheno $P/data/ukbb_pheno_std.tsv --cnv $P/data/ukbb_pgc_xd.cnv --outdir $T
$PY run_pipeline.py --config configs/ukbb_as_aou_test.json --steps prep,gwas,plots,alleles,axes,export
ls -R $P/results/aou_path_test_export | head -50
