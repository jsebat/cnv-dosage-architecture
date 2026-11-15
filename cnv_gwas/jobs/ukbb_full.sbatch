#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH -t 08:00:00
#SBATCH --nice=0
#SBATCH --array=0-7
#SBATCH --job-name=anthro_ukbb
#SBATCH -o /expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/logs/ukbb_full_%A_%a.out
#SBATCH -e /expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/logs/ukbb_full_%A_%a.err
# Full UKBB run through the single runner (configs/ukbb_expanse.json). Tasks 0-1: Task A (GWAS + plots + alleles) for
# height / BMI. Tasks 2-7: Task B axes, trait x TWAS panel set (2,000 genome-shift reps, 18 LOO + no-loci, matched nulls).
set -euo pipefail
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
PY=/home/jsebat/micromamba/envs/dnmt3a-py/bin/python
C=/expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/jobs/snap_ukbb_full/cnv_anthro
cd $C
TR=(height bmi height height height bmi bmi bmi); PN=(- - gtex_nonbrain brain all gtex_nonbrain brain all)
i=$SLURM_ARRAY_TASK_ID
if [ $i -lt 2 ]; then
  $PY run_pipeline.py --config configs/ukbb_expanse.json --steps gwas,plots,alleles --only_trait ${TR[$i]}
else
  $PY run_pipeline.py --config configs/ukbb_expanse.json --steps axes --only_trait ${TR[$i]} --only_panels ${PN[$i]}
fi
