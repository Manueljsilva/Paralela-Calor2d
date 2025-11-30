#!/bin/bash
#SBATCH --job-name=omp_bench
#SBATCH --output=logs/omp_bench_%j.out
#SBATCH --error=logs/omp_bench_%j.err
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --partition=standard

mkdir -p plots logs
rm -f plots/omp_results.csv

for t in 1 2 4 8 16 32; do
    echo "--- Threads = $t ---"
    export OMP_NUM_THREADS=$t
    ./out/omp_bench
done

echo "✓ OpenMP CSV generated: plots/omp_bench.csv"
