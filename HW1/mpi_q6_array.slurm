#!/bin/bash
#SBATCH --job-name=mpi_q6_ntpn16
#SBATCH --time=0:30:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --mem=16G
#SBATCH --output=output.%j
#SBATCH --export=NONE
#SBATCH --get-user-env=L

module load intel

echo "===== Q6: n=1e10, p=64 ====="
echo "nodes=4, ntasks-per-node=16"
mpirun -np 64 ./compute_pi_mpi.exe 10000000000