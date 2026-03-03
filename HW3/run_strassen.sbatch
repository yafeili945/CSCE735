#!/bin/bash
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --job-name=strassen
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --time=0:30:00
#SBATCH --partition=compute
#SBATCH --mem=16G
#SBATCH --output=output.%j
#SBATCH --exclusive

module load intel

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=true
export OMP_PLACES=cores

# 例子：扫 k=10..14，q=4..9（你可改范围）
for k in 10 11 12 13 14; do
  for q in 4 5 6 7 8 9; do
    echo "k=$k q=$q threads=$OMP_NUM_THREADS"
    ./strassen_omp.exe $k $q
  done
done
