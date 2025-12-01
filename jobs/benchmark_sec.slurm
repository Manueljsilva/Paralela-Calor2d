#!/bin/bash
#SBATCH --job-name=sec_bench
#SBATCH --output=logs/sec_bench_%j.out
#SBATCH --error=logs/sec_bench_%j.err
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=all

#SBATCH --mail-user=alejandro.ore@utec.edu.pe
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT_50

mkdir -p plots logs
echo "Running SEQUENTIAL benchmark..."
rm -f plots/sec_results.csv
./out/secuencial
echo "✓ SEQUENTIAL CSV generated: plots/sec_bench.csv"
