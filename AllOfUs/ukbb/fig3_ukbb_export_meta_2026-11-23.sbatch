#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH -t 04:00:00
#SBATCH --nice=100
#SBATCH -J fig3_ukbb_export
#SBATCH -o /expanse/projects/sebat1/jsebat/sacks2026_aou/fig3_2026-11-23/logs/%x.%j.log
# Sacks 2026 dosage paper, Fig. 3 joint model (Jon 2026-11-23): UK Biobank meta-ready export from the finished runs (size2t, all2t, totdose,
# gene mirror; no refit) plus, fitted here, the univariate tests of each class per CNV-size stratum (uni_*: genome-shift null, 2,000 shared
# shifts, LOO as the joint fits) and the binned dose response of every Fig. 3 stratum (fig3_dose_bins.py, computed here: analytic OLS), then the
# UKBB-only meta-analysis and its comparison with the UK Biobank tables the current figure reads; the all-CNV bins are also compared with the
# Fig. 2 tables (class_bins_v5: partitions all and twas_x_con).
# null model: genome shift of the gene annotation (2,000 shared shifts, existing null_runs.npz) and gene-circle shift (gene mirror, 10,000);
# nothing is refitted here: statistics are recomputed from the exported t / per-gene z.
# Unique output dirs: $B/ukbb_export_$SLURM_JOB_ID and $B/meta_ukbb_only_$SLURM_JOB_ID (code snapshot inside the export folder).
set -e
B=/expanse/projects/sebat1/jsebat/sacks2026_aou/fig3_2026-11-23; R=/expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/results
PY=${PY:-/home/jsebat/micromamba/envs/dnmt3a-py/bin/python}; REPO=${REPO:-$B/code/AllOfUs}
W=$B/ukbb_export_$SLURM_JOB_ID; M=$B/meta_ukbb_only_$SLURM_JOB_ID; export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
mkdir -p $W; S=$W/code_snapshot; mkdir -p $S; cp -r $REPO/cnv_anthro $REPO/meta $REPO/ukbb $S/; cd $S/cnv_anthro
$PY ../ukbb/assemble_ukbb_fig3.py --results $R --config configs/ukbb_expanse.json --work $W
$PY run_fig3.py --config $W/ukbb_fig3_config.json --steps fit --only_fit dose_bins,uni_N1_inf,uni_N1_1,uni_N2_inf --jobs 8
$PY run_fig3.py --config $W/ukbb_fig3_config.json --steps export --min_cell 0
(cd $W && tar czf ukbb_fig3_export.tar.gz export && sha256sum ukbb_fig3_export.tar.gz > ukbb_fig3_export.tar.gz.sha256 && ls -la ukbb_fig3_export.tar.gz)
$PY ../meta/meta_fig3.py --cohort UKBB=$W/export/all/fig3 --out $M --code $S/cnv_anthro
set +e
$PY ../meta/compare_fig3_tables.py --a size2t_mirror=$R/size2t_mirror_2026-11-23,all2t_mirror=$R/all2t_mirror_2026-11-23,totdose=$R/totdose_2026-11-23,gene_mirror=$R/gene_mirror_2026-11-22 --b $M > $M/compare_with_ukbb_runs.txt 2>&1
echo "compare exit $?"; cat $M/compare_with_ukbb_runs.txt
$PY ../ukbb/check_bins_vs_fig2.py --new $W/export/all/fig3/dose_bins --old $R/genomewide_dosage_2026-11-22/class_bins_v5 > $M/bins_vs_fig2_class_bins_v5.txt 2>&1
echo "bins vs Fig. 2 exit $?"; cat $M/bins_vs_fig2_class_bins_v5.txt
chmod -R g+rX $W $M
