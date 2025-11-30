#!/bin/bash
#SBATCH --job-name=heat_benchmark
#SBATCH --output=slurm-out/heat_benchmark-%j.out
#SBATCH --error=slurm-out/heat_benchmark-%j.err
#SBATCH --partition=all
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=96       # Total MPI ranks: 192
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=2G

module load mpi
module load gcc   # if needed for OpenMP

echo "Starting benchmark (MPI + OpenMP)..."

# Run the benchmark executable (handles both MPI and OpenMP internally)
# MPI stage: uses all ranks
# OpenMP stage: only rank 0 runs with OMP_NUM_THREADS
export OMP_NUM_THREADS=96   # max threads rank 0 will use
srun -n 192 ./out/version3
