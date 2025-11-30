#!/bin/bash
#SBATCH --job-name=heat_openmp
#SBATCH --output=slurm-out/heat_openmp-%j.out
#SBATCH --error=slurm-out/heat_openmp-%j.err
#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1                  # 1 process
#SBATCH --cpus-per-task=8           # number of OpenMP threads
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=2G

# load modules or set environment as required by your cluster
# module load gcc/12.2 openmpi/4.1.2  # example

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "Running on $SLURM_JOB_NODELIST"
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"

# run the program from a shared filesystem where out/ exists
srun ./out/version1  # replace with program args if any
