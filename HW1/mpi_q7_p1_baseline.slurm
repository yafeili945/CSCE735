#!/bin/bash
#SBATCH --job-name=mpi_q7_p1
#SBATCH --time=0:10:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH --output=output.%j
#SBATCH --export=NONE
#SBATCH --get-user-env=L

module load intel

echo "===== Q7 baseline (p=1): n=1e2,1e4,1e6,1e8 ====="
mpirun -np 1 ./compute_pi_mpi.exe 100
mpirun -np 1 ./compute_pi_mpi.exe 10000
mpirun -np 1 ./compute_pi_mpi.exe 1000000
mpirun -np 1 ./compute_pi_mpi.exe 100000000