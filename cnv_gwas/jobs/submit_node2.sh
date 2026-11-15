#!/bin/bash
# usage (WSL): bash submit_node2.sh <tasks file> <snapshot name> <FROM> <TO> <SLOTS>
set -e
P=/expanse/projects/sebat1/jsebat/ukbb_anthro_cnv
L="/mnt/c/Users/Jonathan Sebat/sebatlab Dropbox/Jon/Claude/UKBB_anthro_CNV"
scp -q "$L"/code/cnv_anthro/*.py expanse:$P/code/cnv_anthro/
scp -q "$L"/code/cnv_anthro/configs/*.json expanse:$P/code/cnv_anthro/configs/
scp -q "$L"/jobs/node2_tasks.sbatch "$L"/jobs/$1 expanse:$P/jobs/
ssh expanse "cd $P && sed -i 's/\r\$//' code/cnv_anthro/*.py code/cnv_anthro/configs/*.json jobs/node2_tasks.sbatch jobs/$1 && rm -rf jobs/$2 && mkdir -p jobs/$2 && cp -r code/cnv_anthro jobs/$2/ && echo PENDING_ddp195: \$(/cm/shared/apps/slurm/current/bin/squeue -A ddp195 -t PD -h | wc -l) && /cm/shared/apps/slurm/current/bin/sbatch --parsable --export=ALL,TASKS=$P/jobs/$1,FROM=$3,TO=$4,SLOTS=$5 jobs/node2_tasks.sbatch"
