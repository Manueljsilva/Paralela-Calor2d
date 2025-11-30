#!/bin/bash
#SBATCH --job-name=mpi_bench
#SBATCH --output=logs/mpi_bench_%j.out
#SBATCH --error=logs/mpi_bench_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=2
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=16
#SBATCH --partition=standard

mkdir -p plots logs
rm -f plots/mpi_results.csv

for np in 1 2 4 8 16 32; do
    echo "--- Processes = $np ---"
    # Adjust tasks per node proportionally for smaller np
    if [ $np -le 16 ]; then
        srun --ntasks=$np --ntasks-per-node=$np ./out/mpi_bench
    else
        srun --ntasks=$np --ntasks-per-node=16 ./out/mpi_bench
    fi
done

echo "✓ MPI CSV generated: plots/mpi_bench.csv"
