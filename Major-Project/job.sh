#!/bin/bash

#SBATCH --job-name=cholesky_run
#SBATCH --output=exp_results.out
#SBATCH --error=exp_results.err
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1               # 默认会分配可用的 GPU，通常是 A100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=03:00:00

# 1. 环境加载
module load intel CUDA

# 2. 确保文件有执行权限（保险起见）
chmod +x original.exe parallel.exe

# 3. 运行实验：并行版扩展性测试 (512 - 4096)
echo "=== Starting Parallel Version Experiments ==="
for size in 512 1024 2048 4096
do
    echo "Matrix Size: $size"
    ./parallel.exe $size
    echo "---------------------------------------"
done

# 4. 运行实验：原版对比测试 (512 - 2048)
echo "=== Starting Original Version Experiments ==="
for size in 512 1024 2048
do
    echo "Matrix Size: $size"
    ./original.exe $size
    echo "---------------------------------------"
done

echo "All experiments finished."