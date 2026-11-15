#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -t 00:40:00
#SBATCH --nice=100
#SBATCH -J grid_smoke
#SBATCH -o /expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/logs/grid_smoke.%j.log
# Smoke test of conga_grid_qt.py (genome-wide univariate grid): height, 2 variants, 5 null reps. Frozen code snapshot.
set -e
P=/expanse/projects/sebat1/jsebat/ukbb_anthro_cnv; R=$P/code/resources
O=$P/results/grid_smoke_$SLURM_JOB_ID; S=$O/code_snapshot; mkdir -p $S
cp $P/code/cnv_anthro/{conga_grid_qt.py,conga_axes_qt.py,common.py} $S/; cp -r $P/code/cnv_anthro/designs $S/
export OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8
cd $S
/home/jsebat/micromamba/envs/dnmt3a-py/bin/python conga_grid_qt.py --pheno $P/data/ukbb_pheno_std.tsv --trait height_rint --cnv $P/data/ukbb_pgc_xd.cnv \
  --twas $R/twas/height.txt --panels gtex_nonbrain --exons $R/conga/hg38_refGene_20200708.exon.txt --gene_cov $R/conga/gene_covariates.tsv \
  --clingen $R/conga/ClinGen_gene_curation_list_GRCh38.tsv --loci $R/conga/omar_loci_all.tsv --loo_runs $R/conga/loo2_runs.tsv \
  --design designs/grid_univariate.json --only "g:all:c3:uns:scaled|adj,pres:n95HI|adj" --nrep 5 --seed 2026 --out $O
