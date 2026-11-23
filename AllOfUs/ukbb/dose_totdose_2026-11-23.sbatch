#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH -t 0:30:00
#SBATCH -J dose_totdose
#SBATCH -o /expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/logs/dose_totdose.%j.log
# Total dose-aligned shift (total genes, all CNV sizes) from the totdose fit (Jon 2026-11-23), analytic.
# null model: genome shift of the gene annotation (existing null_runs.npz); no new fit.
O=/expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/results/totdose_2026-11-23
export OPENBLAS_NUM_THREADS=1
~/micromamba/envs/dnmt3a-py/bin/python /expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/code/dose_summary.py --grid $O/grid --variant totdose \
  --design $O/code_snapshot/designs/class_models_2026-11-22.json --out $O/dose_summary.tsv
