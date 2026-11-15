#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH -t 00:45:00
#SBATCH --nice=0
#SBATCH --job-name=s9_repl
#SBATCH -o /home/jsebat/ukbb_anthro/replication/s9_%j.out
#SBATCH -e /home/jsebat/ukbb_anthro/replication/s9_%j.err
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
/home/jsebat/micromamba/envs/rclogit/bin/Rscript /home/jsebat/ukbb_anthro/replication/replicate_s9_alleles.R
