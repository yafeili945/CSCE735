#!/bin/bash
#SBATCH -J strassen
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH -t 00:30:00
#SBATCH -p compute
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
