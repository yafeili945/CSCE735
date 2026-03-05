#!/bin/bash
#SBATCH --export=NONE
#SBATCH --job-name=strassen_t1_resume
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G
#SBATCH --output=output.%j

module load intel

# 绑核：减少抖动
export OMP_PROC_BIND=true
export OMP_PLACES=cores

exe=./strassen_omp.exe

# ==============================
# Part A: 续跑 1 thread（补齐 k=13 q=8..9 + k=14 q=4..9）
# ==============================
export OMP_NUM_THREADS=1
out1="result_threads_1_resume.txt"

echo "==============================" | tee  "$out1"
echo "RESUME RUN: OMP_NUM_THREADS=1" | tee -a "$out1"
echo "==============================" | tee -a "$out1"
echo "SLURM_JOB_ID=$SLURM_JOB_ID" >> "$out1"
echo "SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}" >> "$out1"
echo "OMP_PROC_BIND=$OMP_PROC_BIND" >> "$out1"
echo "OMP_PLACES=$OMP_PLACES" >> "$out1"
echo "" >> "$out1"

# 你已经跑到 k=13 q=7，所以从 q=8 开始
k=13
for q in 7 8 9; do
  echo "k=$k q=$q threads=$OMP_NUM_THREADS" | tee -a "$out1"
  $exe $k $q | tee -a "$out1"
done

# k=14 全部 q=4..9
k=14
for q in 4 5 6 7 8 9; do
  echo "k=$k q=$q threads=$OMP_NUM_THREADS" | tee -a "$out1"
  $exe $k $q | tee -a "$out1"
done

# ==============================
# Part B: 跑 48 threads 全部（k=10..14, q=4..9）
# ==============================
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-48}
out48="result_threads_${OMP_NUM_THREADS}.txt"

echo "==============================" | tee  "$out48"
echo "FULL RUN: OMP_NUM_THREADS=$OMP_NUM_THREADS" | tee -a "$out48"
echo "==============================" | tee -a "$out48"
echo "SLURM_JOB_ID=$SLURM_JOB_ID" >> "$out48"
echo "SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}" >> "$out48"
echo "OMP_PROC_BIND=$OMP_PROC_BIND" >> "$out48"
echo "OMP_PLACES=$OMP_PLACES" >> "$out48"
echo "" >> "$out48"

for k in 10 11 12 13 14; do
  for q in 4 5 6 7 8 9; do
    echo "k=$k q=$q threads=$OMP_NUM_THREADS" | tee -a "$out48"
    $exe $k $q | tee -a "$out48"
  done
done