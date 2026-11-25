#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH -t 04:00:00
#SBATCH -J gene_mirror
#SBATCH -o /expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/logs/gene_mirror.%j.log
# Dosage paper Fig. 3 (Jon 2026-11-22: "The yellow arrow is the total mirror coming from correlation of Del and Dup effects. Orange arrow is the
# portion explained by GWAS. You need to add a total Del vs Dup to this like Shanta"). gene_mirror_qt.py per trait: per-gene DEL and DUP score z,
# total = r(zDEL, zDUP) adjusted for constrained-gene counts, residual = + signed TWAS count, explained = total - residual; all sizes and
# single-gene CNVs; all / less-constrained / constrained genes; major-locus LOO + minor_all + no_loci.
# null model: DUP-vs-DEL gene-circle shift (>= 10 Mb), 10,000 shifts, same shift set for every run/model (as Shanta gene_mirror_summary.R).
set -e
P=/expanse/projects/sebat1/jsebat/ukbb_anthro_cnv; R=$P/code/resources; PY=/home/jsebat/micromamba/envs/dnmt3a-py/bin/python
O=$P/results/gene_mirror_2026-11-22; S=$O/code_snapshot; cd $S; export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
MAJ=16p11_2_BP4_BP5,16p11_2_BP1_BP3,22q11_21_A_D,15q11_2_BP1_BP2,1q21_1,16p13_11_BP1_BP2,17q12,17p12
for t in height bmi; do
  $PY gene_mirror_qt.py --pheno $P/data/ukbb_pheno_std.tsv --trait ${t}_rint --cnv $P/data/ukbb_pgc_xd.cnv --twas $R/twas/$t.txt --panels gtex_nonbrain \
     --exons $R/conga/hg38_refGene_20200708.exon.txt --gene_cov $R/conga/gene_covariates.tsv --loci $R/growth_loci/loci_sacks34.tsv \
     --loo_runs $R/growth_loci/loo_runs_sacks34.tsv --major_loci $MAJ --nrep 10000 --out $O/$t > $O/$t.log 2>&1 &
done
wait
chmod -R g+rX $O
