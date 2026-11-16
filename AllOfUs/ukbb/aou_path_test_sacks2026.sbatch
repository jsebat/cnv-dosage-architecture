#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH -t 06:00:00
#SBATCH --nice=100
#SBATCH -J aou_path_test_s26
#SBATCH -o %x.%j.log
# Plumbing test of the AoU config before Milad runs it in the Workbench: 5,000 UKBB people re-encoded as AoU inputs, every step of
# configs/aou_sacks2026.json (small nrep, verified UKBB shift files), export, then the meta-analysis against the UKBB reference export.
# Needs make_ukbb_shifts.sbatch and assemble_ukbb_reference.sbatch done under the same OUT_ROOT.
#   sbatch --export=ALL,REPO=<clone>/AllOfUs,OUT_ROOT=<dir>,PY=<python> aou_path_test_sacks2026.sbatch
set -e
PY=${PY:-python3}; REPO=${REPO:?set REPO=<clone>/AllOfUs}; export OUT_ROOT=${OUT_ROOT:-/expanse/projects/sebat1/$USER/sacks2026_aou}
T=$OUT_ROOT/aou_path_test_sacks2026; if [ -e $T ]; then echo "$T exists; use a new OUT_ROOT or remove it"; exit 1; fi
S=$OUT_ROOT/aou_path_test_code_snapshot_$SLURM_JOB_ID; mkdir -p $S; cp -r $REPO/cnv_anthro $REPO/meta $S/
cd $S/cnv_anthro; export OPENBLAS_NUM_THREADS=4
$PY run_pipeline.py --config configs/ukbb_as_aou_sacks2026_test.json --steps gwas,plots,alleles,loci,burden,additivity,scores,grid,joint,export --jobs 8 --threads 32
$PY ../meta/meta_sacks2026.py --ukbb $OUT_ROOT/ukbb_reference_sacks2026/export/all --aou $OUT_ROOT/aou_path_test_sacks2026_export/all --out $OUT_ROOT/aou_path_test_sacks2026_meta
find $OUT_ROOT/aou_path_test_sacks2026_meta -name "*.tsv" | wc -l; head -3 $OUT_ROOT/aou_path_test_sacks2026_meta/joint/all/height/grid_summary.tsv
rm -f $T/*/stratum_pheno.tsv $T/*/joint/strata/pheno_*.tsv      # person-level subsets of the test data
chmod -R g+rX $OUT_ROOT
