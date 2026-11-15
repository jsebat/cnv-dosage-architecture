#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH -t 06:00:00
#SBATCH --nice=0
#SBATCH --array=0-5
#SBATCH --job-name=anthro_axes
#SBATCH -o /expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/logs/axes_%A_%a.out
#SBATCH -e /expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/logs/axes_%A_%a.err
# Task B: ConGA axes (Engchuan et al.) on height / BMI; 2,000 genome-shift reps, 18 LOO + no-loci, matched nulls.
# array: trait x TWAS panel set (gtex_nonbrain = primary; brain, all = sensitivity)
set -euo pipefail
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
PY=/home/jsebat/micromamba/envs/dnmt3a-py/bin/python
P=/expanse/projects/sebat1/jsebat/ukbb_anthro_cnv
C=$P/jobs/snap_axes_full/cnv_anthro
RES=$P/code/resources
TRS=(height height height bmi bmi bmi); PAN=(gtex_nonbrain brain all gtex_nonbrain brain all)
T=${TRS[$SLURM_ARRAY_TASK_ID]}; S=${PAN[$SLURM_ARRAY_TASK_ID]}
OUT=$P/results/taskB_axes/${T}_${S}; mkdir -p $OUT
cd $C
$PY conga_axes_qt.py --pheno $P/data/ukbb_pheno_std.tsv --trait ${T}_rint --cnv $P/data/ukbb_pgc_xd.cnv \
  --twas $RES/twas/$T.txt --panels $S --exons $RES/conga/hg38_refGene_20200708.exon.txt \
  --gene_cov $RES/conga/gene_covariates.tsv --loci $RES/conga/omar_loci_all.tsv --loo_runs $RES/conga/loo2_runs.tsv \
  --nrep 2000 --seed 2026 --out $OUT
