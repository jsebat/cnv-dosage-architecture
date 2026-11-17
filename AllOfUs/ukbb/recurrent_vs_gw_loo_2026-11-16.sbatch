#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=160G
#SBATCH -t 04:00:00
#SBATCH --nice=100
#SBATCH -J rvg_loo
#SBATCH -o /expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/logs/rvg_loo.%j.log
# Jon 2026-11-16: LOO for Supplementary Note 5 (recurrent vs genome-wide).
#  (1) TWAS-gene enrichment of the 34 intervals: each of the 24 LOO units dropped in turn (matched null = same window draws without
#      the dropped intervals), plus 16p11.2 BP1-3/BP2-3/BP4-5 together, and 16p11.2 + 17q12 together (extra runs, not in the tally).
#  (2) Size-bin mirror grid (segdup filter on): loo_mode all (every unit with its own matched null), no_loci, no_overlap34, and the
#      same two extra removals.
# Frozen snapshot from local clone; smoke first; unique output dir.
set -e
P=/expanse/projects/sebat1/jsebat/ukbb_anthro_cnv; R=$P/code/resources; PY=/home/jsebat/micromamba/envs/dnmt3a-py/bin/python
O=$P/results/recurrent_vs_gw_loo_2026-11-16; S=$O/code_snapshot; cd $S; export OPENBLAS_NUM_THREADS=2
G="--exons $R/conga/hg38_refGene_20200708.exon.txt --gene_cov $R/conga/gene_covariates.tsv"
LR="--loo_runs $R/growth_loci/loo_runs_sacks34.tsv"; XR="--extra_runs $S/designs/extra_runs_16p11_2026-11-16.tsv"
# ---- (1) TWAS enrichment LOO (public annotations only)
$PY twas_enrichment_intervals.py $G --twas_height $R/twas/height.txt --twas_bmi $R/twas/bmi.txt --loci $R/growth_loci/loci_sacks34.tsv \
    --segdup_bed $R/segDups.bed $LR $XR --out $O/twas_enrichment > $O/twas_enrichment.log 2>&1 &
# ---- (2) size bins with LOO
fit() {  # trait panels variant nrep out
  $PY conga_grid_qt.py --pheno $P/data/ukbb_pheno_std.tsv --trait $1_rint --cnv $P/data/ukbb_pgc_xd.cnv --twas $R/twas/$1.txt --panels $2 $G \
     --clingen $R/conga/ClinGen_gene_curation_list_GRCh38.tsv --loci $R/growth_loci/loci_sacks34.tsv $LR $XR \
     --design designs/size_bins_2026-11-16.json --only "$3" --nrep $4 --seed 2026 --loo_mode all --no_overlap34 --segdup_bed $R/segDups.bed --out $5
}
fit bmi gtex_nonbrain bp:g:HC:c3:uns:mean 20 $O/smoke > $O/smoke.log 2>&1
echo "smoke OK"; grep -h -E '"N"|genes_excluded_segdup|minutes' $O/smoke/run_info.json; head -c 600 $O/smoke/grid_summary.tsv; echo; rm -rf $O/smoke
T=$O/tasks.txt; : > $T
for v in g:HC:sum:uns:scaled g:HC:c3:sig:mean; do for p in bp ng; do echo "height gtex_nonbrain $p:$v" >> $T; done; done
for v in g:HC:c3:uns:mean g:all:sum:uns:scaled g:HC:c4:uns:mean; do for p in bp ng; do echo "bmi gtex_nonbrain $p:$v" >> $T; done; done
for v in g:HC:c4:uns:mean g:HC:c3:uns:mean; do for p in bp ng; do echo "bmi all $p:$v" >> $T; done; done
export -f fit; export P R PY O G LR XR S
cat $T | xargs -P 14 -L 1 bash -c 'fit $0 $1 "$2" 1000 $O/size_bins/$0_$1/${2//:/_} > $O/size_bins_$0_$1_${2//:/_}.log 2>&1'
wait
echo "tasks: $(wc -l < $T); summaries: $(find $O/size_bins -name grid_summary.tsv | wc -l)"
chmod -R g+rX $O
