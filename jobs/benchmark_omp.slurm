#!/bin/bash
#SBATCH --job-name=omp_bench
#SBATCH --output=logs/omp_bench_%j.out
#SBATCH --error=logs/omp_bench_%j.err
#SBATCH --time=02:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64        # restricted to 64 cores
#SBATCH --partition=all

#SBATCH --mail-user=alejandro.ore@utec.edu.pe
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT_50

set -e

mkdir -p plots logs
rm -f plots/omp_results.csv

# Thread binding for better performance
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

for t in 1 2 4 8 16 32 64; do
    echo "--- Threads = $t ---"
    export OMP_NUM_THREADS=$t
    ./out/version1
done

echo "✓ OpenMP CSV generated: plots/omp_bench.csv"
