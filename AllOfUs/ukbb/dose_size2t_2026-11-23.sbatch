#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH -t 0:30:00
#SBATCH -J dose_size2t
#SBATCH -o /expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/logs/dose_size2t.%j.log
# Dose-aligned contrasts of the size2t joint model (Jon 2026-11-23): D4, weak-or-no-TWAS HC dose, total-genes dose.
# null model: genome shift of the gene annotation (existing null_runs.npz); no new fit.
O=/expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/results/size2t_mirror_2026-11-23
export OPENBLAS_NUM_THREADS=1
~/micromamba/envs/dnmt3a-py/bin/python /expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/code/dose_summary.py --grid $O/grid --variant size2t \
  --design $O/code_snapshot/designs/class_models_2026-11-22.json --out $O/dose_summary.tsv
