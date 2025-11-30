#!/bin/bash
#SBATCH --job-name=heat_benchmark
#SBATCH --output=slurm-out/heat_benchmark-%j.out
#SBATCH --error=slurm-out/heat_benchmark-%j.err
#SBATCH --partition=all          # choose a partition that has enough CPUs
#SBATCH --nodes=2                # 2 nodes for MPI
#SBATCH --ntasks-per-node=64     # MPI ranks per node
#SBATCH --cpus-per-task=1        # 1 CPU per MPI rank
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=2G

module load mpi
module load gcc

echo "Starting job on $SLURM_JOB_NUM_NODES nodes, $SLURM_NTASKS tasks total"

export OMP_NUM_THREADS=128
srun -n 128 ./out/version3
