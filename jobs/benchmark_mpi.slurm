#!/bin/bash
#SBATCH --job-name=mpi_bench
#SBATCH --output=logs/mpi_bench_%j.out
#SBATCH --error=logs/mpi_bench_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1-8
#SBATCH --ntasks-per-node=16
#SBATCH --partition=all

mkdir -p plots logs
rm -f plots/mpi_results.csv

CORES_PER_NODE=16

for np in 1 2 4 8 16 32 64 128; do
    echo "--- Processes = $np ---"

    if [ $np -le $CORES_PER_NODE ]; then
        TASKS_PER_NODE=$np
    else
        TASKS_PER_NODE=$CORES_PER_NODE
    fi

    echo "Running with $np tasks, $TASKS_PER_NODE tasks per node"
    mpirun -np $np --map-by ppr:$TASKS_PER_NODE:node ./out/mpi_bench
done

echo "✓ MPI CSV generated: plots/mpi_bench.csv"
