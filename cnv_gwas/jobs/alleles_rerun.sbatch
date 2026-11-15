#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH --nice=0
#SBATCH --job-name=anthro_alleles
#SBATCH -o /expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/logs/alleles_rerun_%j.out
#SBATCH -e /expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/logs/alleles_rerun_%j.err
set -euo pipefail
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
PY=/home/jsebat/micromamba/envs/dnmt3a-py/bin/python
C=/expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/jobs/snap_alleles_rerun/cnv_anthro
cd $C
$PY run_pipeline.py --config configs/ukbb_expanse.json --steps alleles
