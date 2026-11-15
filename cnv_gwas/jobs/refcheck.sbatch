#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH --nice=0
#SBATCH --job-name=anthro_refcheck
#SBATCH -o /expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/logs/refcheck_%j.out
#SBATCH -e /expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/logs/refcheck_%j.err
# Fast vs reference (naive dense OLS per test) on UKBB chr22 height + axes model; timings recorded.
set -euo pipefail
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
PY=/home/jsebat/micromamba/envs/dnmt3a-py/bin/python
P=/expanse/projects/sebat1/jsebat/ukbb_anthro_cnv; C=$P/jobs/snap_refcheck/cnv_anthro; R=$P/code/resources
cd $C
OUT=$P/results/reference_check; mkdir -p $OUT
/usr/bin/time -v $PY cnv_gwas_qt.py --pheno $P/data/ukbb_pheno_std.tsv --trait height_rint --cnv $P/data/ukbb_pgc_xd.cnv \
  --bp_dir /expanse/projects/sebat1/oshanta/home/IndividualDisorder/Freeze_2024/inputs/bp_files --chroms 22-22 \
  --known_loci $R/omar_alleles/gws_all_cnvs_known_loci.txt --segdups $R/segDups.bed --out $OUT/fast_chr22 2> $OUT/fast_chr22_time.txt
grep -E "Elapsed|Maximum resident" $OUT/fast_chr22_time.txt
$PY reference_check.py --pheno $P/data/ukbb_pheno_std.tsv --cnv $P/data/ukbb_pgc_xd.cnv --gwas_out $OUT/fast_chr22 --chrom 22 --n_probes 200 \
  --axes --twas $R/twas/height.txt --exons $R/conga/hg38_refGene_20200708.exon.txt --gene_cov $R/conga/gene_covariates.tsv --out $OUT
