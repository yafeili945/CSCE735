#!/bin/bash
#SBATCH --export=NONE
#SBATCH --job-name=strassen_t48
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=64G
#SBATCH --output=output.%j

module load intel

# 绑核：减少抖动
export OMP_PROC_BIND=true
export OMP_PLACES=cores

# 只跑 48 threads（用满本节点）
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-48}

echo "=============================="
echo "RUN: OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo "=============================="

# 输出到单独文件
out="result_threads_test_${OMP_NUM_THREADS}.txt"
echo "Writing to $out"
echo "SLURM_JOB_ID=$SLURM_JOB_ID" >  "$out"
echo "SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}" >> "$out"
echo "OMP_PROC_BIND=$OMP_PROC_BIND" >> "$out"
echo "OMP_PLACES=$OMP_PLACES" >> "$out"
echo "" >> "$out"

for k in 10 11 12 13 14; do
  for q in 4 5 6 7 8 9; do
    echo "k=$k q=$q threads=$OMP_NUM_THREADS" | tee -a "$out"
    ./strassen_omp.exe $k $q | tee -a "$out"
  done
done