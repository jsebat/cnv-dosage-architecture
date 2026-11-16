#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH --nice=100
#SBATCH -J meta_sacks2026
#SBATCH -o %x.%j.log
# Federated meta-analysis at UCSD: the All of Us summary statistics (export_aggregate tarball egressed from the Workbench, with its
# null_runs.npz test-statistic arrays) are meta-analysed here with the UKBB reference export. Nothing individual-level is read.
#   sbatch --export=ALL,REPO=<clone>/AllOfUs,AOU_TAR=<aou_sacks2026_export_<date>.tar.gz>,PY=<python> meta_sacks2026_expanse.sbatch
# Optional: UKBB=<UKBB export stratum folder> (default: the A2 export), OUT=<new meta folder>.
set -e
PY=${PY:-python3}; REPO=${REPO:?set REPO=<clone>/AllOfUs}; AOU_TAR=${AOU_TAR:?set AOU_TAR=<AoU export tarball>}
B=/expanse/projects/sebat1/jsebat/sacks2026_aou
UKBB=${UKBB:-$B/ukbb_reference_sacks2026/export/all}
TAG=$(basename $AOU_TAR .tar.gz); A=$(dirname $AOU_TAR)/$TAG; OUT=${OUT:-$(dirname $AOU_TAR)/meta_$TAG}
if [ -e $OUT ]; then echo "$OUT exists; set a new OUT"; exit 1; fi
[ -f $AOU_TAR.sha256 ] && (cd $(dirname $AOU_TAR) && sha256sum -c $(basename $AOU_TAR).sha256)
mkdir -p $A; tar xzf $AOU_TAR -C $A
AOU=$(find $A -maxdepth 3 -type d -path '*export_aggregate/all' | head -1); [ -n "$AOU" ] || { echo "no export_aggregate/all in $AOU_TAR"; exit 1; }
S=$(dirname $OUT)/meta_code_snapshot_$SLURM_JOB_ID; mkdir -p $S; cp -r $REPO/cnv_anthro $REPO/meta $S/
cd $S/cnv_anthro; export OPENBLAS_NUM_THREADS=8
$PY ../meta/meta_sacks2026.py --ukbb $UKBB --aou $AOU --out $OUT
echo "tables: $(find $OUT -name '*.tsv' | wc -l)"
echo "grid rows with shared_shifts False (must be 0):"; for f in $(find $OUT -name 'grid_summary.tsv'); do awk -F'\t' -v f=$f 'NR==1{for(i=1;i<=NF;i++) if($i=="shared_shifts") c=i; next} c && $c!="True"{n++} END{if(n) print f, n}' $f; done
rm -rf $S; chmod -R g+rX $OUT $A
