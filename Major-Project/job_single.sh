#!/bin/bash

#SBATCH --job-name=single_4096
#SBATCH --output=single_4096.out
#SBATCH --error=single_4096.err
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=12:00:00

# 1. 环境加载
module load intel CUDA

# 2. 权限设置
chmod +x original.exe parallel_8.exe parallel_16.exe parallel_32.exe

# 4. 实验：原版对比测试 (512 - 2048)
echo ""
echo "=========================================================="
echo "EXPERIMENT: Original Version (Single Thread reference)"
echo "=========================================================="
for size in 4096
do
    echo "Original Matrix Size: $size"
    ./original.exe $size
    echo "----------------------------------------------------------"
done

echo "Full experiments finished."