#!/bin/bash
#SBATCH --job-name=benchmark
#SBATCH --output=slurm-out/benchmark-%j.out
#SBATCH --error=slurm-out/benchmark-%j.err
#SBATCH --partition=all           # choose a big partition with more nodes
#SBATCH --nodes=5                 # number of nodes to use
#SBATCH --ntasks-per-node=8       # MPI ranks per node
#SBATCH --time=02:00:00           # estimated runtime
#SBATCH --mem-per-cpu=4G          # memory per CPU
#SBATCH --exclusive               # optional: no other jobs on these nodes

# Load modules if needed
# module load mpi

# Total MPI ranks
TOTAL_RANKS=$(( SLURM_NNODES * SLURM_NTASKS_PER_NODE ))

# Set OpenMP threads per rank
export OMP_NUM_THREADS=4

echo "Running hybrid MPI+OpenMP benchmark:"
echo "Nodes: $SLURM_NNODES, MPI ranks/node: $SLURM_NTASKS_PER_NODE, Total ranks: $TOTAL_RANKS"
echo "OpenMP threads per rank: $OMP_NUM_THREADS"

# Run the benchmark
srun -n $TOTAL_RANKS ./out/version3
