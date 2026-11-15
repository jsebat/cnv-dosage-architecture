#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-compute
#SBATCH --qos=ind-normal
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=240G
#SBATCH -t 08:00:00
#SBATCH --nice=100
#SBATCH -J anthro_node2
#SBATCH -o /expanse/projects/sebat1/jsebat/ukbb_anthro_cnv/logs/node2.%j.log
# Lab job shape (Jon): 2 whole ind-compute nodes, one srun worker per node; each worker runs SLOTS task lines in
# parallel from TASKS (one shell command per line), lines FROM..TO, claimed atomically with mkdir claims/<line>.
# Submit: sbatch --export=ALL,TASKS=<file>,FROM=1,TO=<n>,SLOTS=<k> node2_tasks.sbatch   (check squeue -A ddp195 -t PD first)
export TASKS FROM=${FROM:-1} TO=${TO:-1000000} SLOTS=${SLOTS:-8}
CL=$(dirname $TASKS)/claims.$(basename $TASKS); mkdir -p $CL
export CL
srun --ntasks=2 --ntasks-per-node=1 --cpus-per-task=128 bash -c '
  export OPENBLAS_NUM_THREADS=$((128 / SLOTS)) OMP_NUM_THREADS=$((128 / SLOTS))
  seq $FROM $TO | xargs -P $SLOTS -I{} bash -c "
    L=\$(sed -n {}p $TASKS); [ -z \"\$L\" ] && exit 0
    mkdir $CL/{} 2>/dev/null || exit 0
    bash -c \"\$L\" > $CL/{}/log 2>&1 && touch $CL/{}/DONE || echo FAILED line {} \$(hostname)"
  echo "$(date +%H:%M) $(hostname) worker finished $FROM-$TO"'
