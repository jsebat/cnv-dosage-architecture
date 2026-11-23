#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --qos=ind-shared-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH -J meta_fig3
#SBATCH -o /expanse/projects/sebat1/jsebat/sacks2026_aou/fig3_2026-11-23/logs/%x.%j.log
# Sacks 2026 dosage paper, Fig. 3 joint model: federated meta-analysis at UCSD of the UK Biobank export and one or more All of Us Fig. 3
# exports (aggregate tarballs egressed from the Workbench: t of every term in every run and shift replicate, per-gene z; no participant data).
#   sbatch --export=ALL,AOU_TARS="<aou_fig3 tar.gz>[ <array tar.gz>]",LABELS="AoU_WGS[ AoU_array]",UKBB=<ukbb_export_<id>/export/all/fig3> meta_fig3_expanse.sbatch
# Output: $B/meta_<labels>_$SLURM_JOB_ID (unique). Figure: GW2_FIG3_META=<that folder> python figs_dosage_story.py D
# null model: genome shift of the gene annotation, the same 2,000 UK Biobank shifts in every cohort (meta formed replicate by replicate);
# gene mirror: gene-circle shift on the meta per-gene z; total-genes terms analytic. meta_fig3.py stops if cohorts did not share shifts.
set -e
B=/expanse/projects/sebat1/jsebat/sacks2026_aou/fig3_2026-11-23
PY=${PY:-/home/jsebat/micromamba/envs/dnmt3a-py/bin/python}; REPO=${REPO:-$B/code/AllOfUs}
UKBB=${UKBB:?set UKBB=<UKBB fig3 export>/all/fig3}; AOU_TARS=${AOU_TARS:?set AOU_TARS}; LABELS=${LABELS:?set LABELS (one per tarball)}
read -ra TS <<< "$AOU_TARS"; read -ra LS <<< "$LABELS"; [ ${#TS[@]} -eq ${#LS[@]} ] || { echo "AOU_TARS and LABELS differ in length"; exit 1; }
OUT=$B/meta_UKBB_$(echo $LABELS | tr ' ' '_')_$SLURM_JOB_ID; X=$B/aou_extracted_$SLURM_JOB_ID; mkdir -p $X
S=$X/code_snapshot; mkdir -p $S; cp -r $REPO/cnv_anthro $REPO/meta $S/
ARGS="--cohort UKBB=$UKBB"
for i in "${!TS[@]}"; do
  T=${TS[$i]}; L=${LS[$i]}
  [ -f $T.sha256 ] && (cd $(dirname $T) && sha256sum -c $(basename $T).sha256)
  mkdir -p $X/$L; tar xzf $T -C $X/$L
  if find $X/$L -type f | grep -i -E "pheno|genotype|\.cnv$"; then echo "person-level file name in $T: stop"; exit 1; fi
  D=$(find $X/$L -type d -path '*all/fig3' | head -1); [ -n "$D" ] || { echo "no all/fig3 in $T"; exit 1; }
  ARGS="$ARGS --cohort $L=$D"
done
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
$PY $S/meta/meta_fig3.py $ARGS --out $OUT --code $S/cnv_anthro
cat $OUT/meta_checks.tsv; chmod -R g+rX $OUT $X
