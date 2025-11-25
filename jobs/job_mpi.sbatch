#!/bin/bash
#SBATCH --job-name=heat_mpi
#SBATCH --output=slurm-out/heat_mpi-%j.out
#SBATCH --error=slurm-out/heat_mpi-%j.err
#SBATCH --partition=normal
#SBATCH --nodes=2                     # number of nodes
#SBATCH --ntasks-per-node=8           # MPI ranks per node
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=2G

# module load mpi  # if needed

TOTAL_RANKS=$(( SLURM_NNODES * SLURM_NTASKS_PER_NODE ))
echo "Nodes: $SLURM_NNODES, ranks per node: $SLURM_NTASKS_PER_NODE, total ranks: $TOTAL_RANKS"

# srun will launch one MPI rank per task on allocated nodes
srun -n $TOTAL_RANKS ./out/version2  # or add program args
