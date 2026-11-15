#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -t 00:40:00
#SBATCH --nice=100
#SBATCH -J cs_smoke
#SBATCH -o /expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/logs/cs_smoke.%j.log
# Smoke test: call-size classes in conga_grid_qt.py (height, one variant, 12 reps). Frozen snapshot made by the submit script.
export OPENBLAS_NUM_THREADS=8
bash -c "$(cat /expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/jobs/smoke_callsize_cmd.txt)"
