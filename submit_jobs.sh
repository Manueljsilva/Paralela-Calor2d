#!/bin/bash

# Array of SLURM scripts
SCRIPTS=("jobs/benchmark_sec.slurm" "jobs/benchmark_mpi.slurm" "jobs/benchmark_omp.slurm")

# Submit each script
for script in "${SCRIPTS[@]}"; do
    if [[ -f "$script" ]]; then
        echo "Submitting $script..."
        sbatch "$script"
    else
        echo "Warning: $script not found!"
    fi
done
