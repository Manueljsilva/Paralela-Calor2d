#!/bin/bash
#SBATCH --job-name=heat_hybrid
#SBATCH --output=slurm-out/heat_hybrid-%j.out
#SBATCH --error=slurm-out/heat_hybrid-%j.err
#SBATCH --partition=normal
#SBATCH --nodes=2                  # number of nodes
#SBATCH --ntasks-per-node=4        # MPI ranks per node
#SBATCH --cpus-per-task=8          # cores reserved per MPI rank (for OpenMP threads)
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=2G

# optional module loads
# module load gcc openmpi

# Set OpenMP threads per MPI rank
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

TOTAL_RANKS=$(( SLURM_NNODES * SLURM_NTASKS_PER_NODE ))
echo "Nodes: $SLURM_NNODES, ranks per node: $SLURM_NTASKS_PER_NODE, cpus per rank: $SLURM_CPUS_PER_TASK"
echo "Total MPI ranks: $TOTAL_RANKS, OMP_NUM_THREADS=$OMP_NUM_THREADS"

# Launch: srun ensures correct mapping and cpu-binding
srun -n $TOTAL_RANKS ./out/version3  # pass args if needed
