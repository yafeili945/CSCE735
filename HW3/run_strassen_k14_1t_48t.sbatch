#!/bin/bash
#SBATCH --export=NONE
#SBATCH --job-name=strassen_k14_1t_48t
#SBATCH --time=12:00:00
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

# =========================
# Part A: 1 thread, k=13 q=9
# =========================
export OMP_NUM_THREADS=1
out1_139="result_threads_1_k13_q9.txt"

echo "==============================" | tee  "$out1_139"
echo "RUN: OMP_NUM_THREADS=$OMP_NUM_THREADS (k=13,q=9)" | tee -a "$out1_139"
echo "==============================" | tee -a "$out1_139"
echo "SLURM_JOB_ID=$SLURM_JOB_ID" >> "$out1_139"
echo "SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}" >> "$out1_139"
echo "OMP_PROC_BIND=$OMP_PROC_BIND" >> "$out1_139"
echo "OMP_PLACES=$OMP_PLACES" >> "$out1_139"
echo "" >> "$out1_139"

echo "k=13 q=9 threads=$OMP_NUM_THREADS" | tee -a "$out1_139"
$exe 13 9 | tee -a "$out1_139"


# =========================
# Part B: k=14 full sweep for BOTH 1 thread and 48 threads
# =========================

# ---- B1) 1 thread for k=14, q=4..9 ----
export OMP_NUM_THREADS=1
out1_14="result_threads_1_k14_allq.txt"

echo "==============================" | tee  "$out1_14"
echo "RUN: OMP_NUM_THREADS=$OMP_NUM_THREADS (k=14,q=4..9)" | tee -a "$out1_14"
echo "==============================" | tee -a "$out1_14"
echo "SLURM_JOB_ID=$SLURM_JOB_ID" >> "$out1_14"
echo "SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}" >> "$out1_14"
echo "OMP_PROC_BIND=$OMP_PROC_BIND" >> "$out1_14"
echo "OMP_PLACES=$OMP_PLACES" >> "$out1_14"
echo "" >> "$out1_14"

for q in 4 5 6 7 8 9; do
  echo "k=14 q=$q threads=$OMP_NUM_THREADS" | tee -a "$out1_14"
  $exe 14 $q | tee -a "$out1_14"
done

# ---- B2) 48 threads for k=14, q=4..9 ----
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-48}
out48_14="result_threads_${OMP_NUM_THREADS}_k14_allq.txt"

echo "==============================" | tee  "$out48_14"
echo "RUN: OMP_NUM_THREADS=$OMP_NUM_THREADS (k=14,q=4..9)" | tee -a "$out48_14"
echo "==============================" | tee -a "$out48_14"
echo "SLURM_JOB_ID=$SLURM_JOB_ID" >> "$out48_14"
echo "SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}" >> "$out48_14"
echo "OMP_PROC_BIND=$OMP_PROC_BIND" >> "$out48_14"
echo "OMP_PLACES=$OMP_PLACES" >> "$out48_14"
echo "" >> "$out48_14"

for q in 4 5 6 7 8 9; do
  echo "k=14 q=$q threads=$OMP_NUM_THREADS" | tee -a "$out48_14"
  $exe 14 $q | tee -a "$out48_14"
done