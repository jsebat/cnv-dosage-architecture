#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G
#SBATCH -t 5:00:00
#SBATCH --nice=100
#SBATCH -J gw_size2t
#SBATCH -o /expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/logs/gw_size2t.%j.log
# size2t run (Jon 2026-11-23): per CNV size n_all + up/down x LC/HC + other HC (reference other LC), all four CNV x TWAS pairs.
# Dosage paper Fig. 3 (Jon 2026-11-22: "show the mirror ... a separate grid ... like Shanta Fig. 3E-F"): mirror and dosage-effect grids with
# rows = CNV trait (height, BMI) and columns = TWAS trait (height, BMI). Runs joint6 (six TWAS x constraint classes) and size4 (the same by CNV
# size) for all four CNV-trait x TWAS-trait pairings, 2,000 shared shifts (the CNV trait's shift file), LOO = drivers with fixed major loci.
# Summary: mirror4_summary.py (four-tile mirror = (DUP arm - DEL arm)/sqrt2 of the signed TWAS terms, with the same statistic in every null
# replicate). null model: genome shift of the gene annotation, shared shifts; this file only runs the fits. Unique output dir; code snapshot.
set -e
P=/expanse/projects/sebat1/jsebat/ukbb_anthro_cnv; R=$P/code/resources; PY=/home/jsebat/micromamba/envs/dnmt3a-py/bin/python
O=$P/results/size2t_mirror_2026-11-23; S=$O/code_snapshot; cd $S; export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
MAJ=16p11_2_BP4_BP5,16p11_2_BP1_BP3,22q11_21_A_D,15q11_2_BP1_BP2,1q21_1,16p13_11_BP1_BP2,17q12,17p12
test ! -e $O/tasks.txt || { echo "output dir already used"; exit 1; }
fit() {  # cnv_trait twas_trait variant
  $PY conga_grid_qt.py --pheno $P/data/ukbb_pheno_std.tsv --trait $1_rint --cnv $P/data/ukbb_pgc_xd.cnv --twas $R/twas/$2.txt --panels gtex_nonbrain \
     --exons $R/conga/hg38_refGene_20200708.exon.txt --gene_cov $R/conga/gene_covariates.tsv --clingen $R/conga/ClinGen_gene_curation_list_GRCh38.tsv \
     --loci $R/growth_loci/loci_sacks34.tsv --loo_runs $R/growth_loci/loo_runs_sacks34.tsv --design designs/class_models_2026-11-22.json --only "$3" \
     --nrep 2000 --seed 2026 --loo_mode drivers --major_loci $MAJ --perms_file $R/shared_shifts/$1_shifts.npy --out $O/grid/$1_x_$2/$3
}
export -f fit; export P R PY O MAJ S
T=$O/tasks.txt; : > $T
for c in height bmi; do for t in height bmi; do for v in size2t; do echo "$c $t $v" >> $T; done; done; done
cat $T | xargs -P 8 -L 1 bash -c 'fit $0 $1 $2 > $O/fit_$0_x_$1_$2.log 2>&1'
echo "fits: $(find $O/grid -name grid_summary.tsv | wc -l) of 4"
cd $O && $PY $S/mirror4_summary.py --grid grid --design $S/designs/class_models_2026-11-22.json --out mirror4_cross_summary.tsv > mirror4.log 2>&1
chmod -R g+rX $O
