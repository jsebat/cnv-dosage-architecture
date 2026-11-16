#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH --nice=100
#SBATCH -J ukbb_reference
#SBATCH -o %x.%j.log
# UKBB reference export for the meta-analysis: symlinked results tree (no copies of the UKBB results) + aggregate export + tarball.
#   sbatch --export=ALL,REPO=<clone>/AllOfUs,OUT_ROOT=<dir>,PY=<python> assemble_ukbb_reference.sbatch
set -e
PY=${PY:-python3}; REPO=${REPO:?set REPO=<clone>/AllOfUs}; OUT_ROOT=${OUT_ROOT:-/expanse/projects/sebat1/$USER/sacks2026_aou}
O=$OUT_ROOT/ukbb_reference_sacks2026; if [ -e $O/export ]; then echo "$O/export exists; use a new OUT_ROOT or remove it"; exit 1; fi
S=$O/code_snapshot; mkdir -p $S; cp $REPO/cnv_anthro/*.py $S/; cp -r $REPO/cnv_anthro/designs $S/; cp $REPO/ukbb/assemble_ukbb_reference.py $S/
$PY $S/assemble_ukbb_reference.py --code $S --out $O --python $PY
tar czf $O/ukbb_reference_sacks2026_export.tar.gz -C $O/export all && sha256sum $O/ukbb_reference_sacks2026_export.tar.gz > $O/ukbb_reference_sacks2026_export.sha256
ls -la $O; chmod -R g+rX $O
