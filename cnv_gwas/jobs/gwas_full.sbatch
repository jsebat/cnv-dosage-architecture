#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH --nice=0
#SBATCH --array=0-1
#SBATCH --job-name=anthro_gwas
#SBATCH -o /expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/logs/gwas_%A_%a.out
#SBATCH -e /expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/logs/gwas_%A_%a.err
# Task A: probe-level CNV GWAS (Omar's probes/filters/thresholds, linear model) + allele-level tests, height and BMI.
set -euo pipefail
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
PY=/home/jsebat/micromamba/envs/dnmt3a-py/bin/python
P=/expanse/projects/sebat1/jsebat/ukbb_anthro_cnv
C=$P/jobs/snap_gwas_full/cnv_anthro
RES=$P/code/resources
FD=/expanse/projects/sebat1/oshanta/home/CrossDisorder/Freeze_2024/inputs/data_10perc/filteredCNVs
FILT=$FD/filter.AB.TELE.CENTRO.bed,$FD/filter.CONSENTYPE.bed,$FD/filter.SEGDUP.bed,$FD/filter.FREQUENCY_10percent.bed,$FD/filter.FREQUENCY_CROSS_PF.bed
TR=(height bmi); T=${TR[$SLURM_ARRAY_TASK_ID]}
OUT=$P/results/taskA_gwas/$T; mkdir -p $OUT
cd $C
$PY cnv_gwas_qt.py --pheno $P/data/ukbb_pheno_std.tsv --trait ${T}_rint --cnv $P/data/ukbb_pgc_xd.cnv \
  --bp_dir /expanse/projects/sebat1/oshanta/home/IndividualDisorder/Freeze_2024/inputs/bp_files --chroms 1-23 \
  --known_loci $RES/omar_alleles/gws_all_cnvs_known_loci.txt --segdups $RES/segDups.bed --filtered_cnv_beds $FILT \
  --out $OUT --label "UKBB $T"
$PY plots.py --gwas_dir $OUT --title "UKBB $T (RINT within sex)" --out $OUT/miami_$T.png
