#!/bin/bash
#SBATCH --job-name=mpi_q5_q7
#SBATCH --time=0:30:00
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=4
#SBATCH --mem=8G
#SBATCH --output=output.%j
#SBATCH --export=NONE
#SBATCH --get-user-env=L

module load intel

echo "===== Q5: n=1e8, p=1,2,4,8,16,32,64 (ntasks-per-node=4) ====="
mpirun -np 1  ./compute_pi_mpi.exe 100000000
mpirun -np 2  ./compute_pi_mpi.exe 100000000
mpirun -np 4  ./compute_pi_mpi.exe 100000000
mpirun -np 8  ./compute_pi_mpi.exe 100000000
mpirun -np 16 ./compute_pi_mpi.exe 100000000
mpirun -np 32 ./compute_pi_mpi.exe 100000000
mpirun -np 64 ./compute_pi_mpi.exe 100000000

echo "===== Q7 (p=64): n=1e2,1e4,1e6,1e8 (ntasks-per-node=4) ====="
mpirun -np 64 ./compute_pi_mpi.exe 100
mpirun -np 64 ./compute_pi_mpi.exe 10000
mpirun -np 64 ./compute_pi_mpi.exe 1000000
mpirun -np 64 ./compute_pi_mpi.exe 100000000