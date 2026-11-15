#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -t 00:30:00
#SBATCH --nice=100
#SBATCH -J locus_eff
#SBATCH -o /expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/logs/locus_effects.%j.log
# Per-locus UKBB DEL/DUP effects for the growth locus-level analysis (new UKBB loci + 34 Sacks intervals as calibration). Frozen code snapshot.
set -e
P=/expanse/projects/sebat1/jsebat/ukbb_anthro_cnv
O=$P/results/locus_effects_newloci_2026-11-15; mkdir -p $O
S=$O/code_snapshot; mkdir -p $S; cp $P/code/cnv_anthro/locus_effects.py $P/code/cnv_anthro/common.py $S/
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
PY=/home/jsebat/micromamba/envs/dnmt3a-py/bin/python
cd $S
$PY locus_effects.py --pheno $P/data/ukbb_pheno_std.tsv --cnv $P/data/ukbb_pgc_xd.cnv \
  --loci $O/loci_for_ukbb_effects.tsv --out $O/locus_effects_ukbb.tsv
chmod -R g+rX $O
