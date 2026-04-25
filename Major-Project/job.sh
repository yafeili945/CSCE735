#!/bin/bash

#SBATCH --job-name=cholesky_full
#SBATCH --output=exp_results.out
#SBATCH --error=exp_results.err
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=03:00:00

# 1. Load environment
module load intel CUDA

# 2. Set executable permissions
chmod +x original.exe parallel_8.exe parallel_16.exe parallel_32.exe

# 3. Experiment: full comparison of parallel versions (512 - 4096)
# This allows plotting three curves (8x8, 16x16, 32x32) to compare performance by matrix size
for size in 512 1024 2048 4096
do
    echo "=========================================================="
    echo "MATRIX SIZE: $size"
    echo "=========================================================="
    
    echo "Running Parallel 8x8..."
    ./parallel_8.exe $size
    
    echo "Running Parallel 16x16..."
    ./parallel_16.exe $size
    
    echo "Running Parallel 32x32..."
    ./parallel_32.exe $size
    echo "----------------------------------------------------------"
done

# 4. Experiment: baseline original-version comparison (512 - 2048)
echo ""
echo "=========================================================="
echo "EXPERIMENT: Original Version (Single Thread reference)"
echo "=========================================================="
for size in 512 1024 2048
do
    echo "Original Matrix Size: $size"
    ./original.exe $size
    echo "----------------------------------------------------------"
done

echo "Full experiments finished."