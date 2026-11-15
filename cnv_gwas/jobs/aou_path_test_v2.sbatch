#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH -t 02:00:00
#SBATCH --nice=100
#SBATCH -J aoutest2
#SBATCH -o /expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/logs/aoutest2.%j.log
# AoU code path v2 on 5,000 re-encoded UKBB people: strata + loci + grid (shared UKBB shifts) + export + meta tests. Frozen snapshot.
set -e
P=/expanse/projects/sebat1/jsebat/ukbb_anthro_cnv; PY=/home/jsebat/micromamba/envs/dnmt3a-py/bin/python
S=$P/jobs/snap_aoutest2_$SLURM_JOB_ID; mkdir -p $S; cp -r $P/code/cnv_anthro $S/
cd $S/cnv_anthro
/usr/bin/time -v $PY run_pipeline.py --config configs/ukbb_as_aou_test.json --steps gwas,alleles,loci,grid,export 2> $P/logs/aoutest2_time.$SLURM_JOB_ID.txt
grep -E "Maximum resident|Elapsed" $P/logs/aoutest2_time.$SLURM_JOB_ID.txt
E=$P/results/aou_path_test_v2_export; M=$P/results/aou_path_test_v2_meta; mkdir -p $M
$PY meta_sumstats.py --cohort F=$E/female --cohort Mx=$E/male --trait height --out $M/sumstats_height_female_male
for t in height bmi; do for d in grid_univariate grid_callsize joint_2026-11-15; do
  [ -d $E/female/grid/$d/$t ] && $PY meta_grid.py --cohort F=$E/female/grid/$d/$t --cohort Mx=$E/male/grid/$d/$t --out $M/meta_${d}_${t}.tsv
done; done
ls -la $M; head -5 $M/meta_grid_univariate_height.tsv
