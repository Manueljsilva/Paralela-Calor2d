#!/bin/bash
#SBATCH --job-name=seq_bench
#SBATCH --output=logs/seq_bench_%j.out
#SBATCH --error=logs/seq_bench_%j.err
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=standard

mkdir -p plots logs
echo "Running SEQUENTIAL benchmark..."
rm -f plots/seq_results.csv
./out/seq_bench
echo "✓ SEQUENTIAL CSV generated: plots/seq_bench.csv"
