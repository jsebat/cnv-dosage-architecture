#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -t 03:00:00
#SBATCH --nice=100
#SBATCH -J ukbb_shifts
#SBATCH -o %x.%j.log
# UKBB shift lists for the UKBB + All of Us meta-analysis (Sacks 2026), verified against the perms_sha1 of every UKBB run they belong to.
# Frozen code snapshot; unique output dir; person-level strata files deleted at the end. Reads the UKBB data in jsebat's project folder
# (group-readable); writes under OUT_ROOT (default: your own folder). Submit from the Expanse login node, in the folder for the log:
#   sbatch --test-only --export=ALL,REPO=<clone>/AllOfUs,OUT_ROOT=<dir>,PY=<python> make_ukbb_shifts.sbatch
#   sbatch --export=ALL,REPO=<clone>/AllOfUs,OUT_ROOT=<dir>,PY=<python> make_ukbb_shifts.sbatch
set -e
P=/expanse/projects/sebat1/jsebat/ukbb_anthro_cnv; PY=${PY:-python3}      # python 3 with numpy, scipy, pandas
REPO=${REPO:?set REPO=<clone of jsebat/cnv-growth-bmi, branch cnv-gwas>/AllOfUs}; OUT_ROOT=${OUT_ROOT:-/expanse/projects/sebat1/$USER/sacks2026_aou}
O=$OUT_ROOT/shared_shifts_sacks2026_2026-11-16; S=$O/code_snapshot; mkdir -p $S
cp $REPO/cnv_anthro/*.py $S/; cp -r $REPO/cnv_anthro/designs $S/; cp $REPO/ukbb/make_ukbb_shifts.py $S/
export OPENBLAS_NUM_THREADS=8
if [ -e $O/shared_shifts/SHIFTS_MANIFEST.tsv ]; then echo "$O already has results; use a new OUT_ROOT"; exit 1; fi
$PY $S/make_ukbb_shifts.py --code $S --out $O/shared_shifts --work $O/work --python $PY
chmod -R g+rX $O
