#!/bin/bash
#SBATCH --time=000:10:00
#SBATCH --account=def-akhadra
#SBATCH --mem=4G
#SBATCH --mail-user=niklas.brake@mail.mcgill.ca
#SBATCH --mail-type=FAIL
#SBATCH --cpus-per-task=2
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $1)
DIR=$1
folder1="$DIR/presynaptic_network/spikeTimesLong.csv"
folder2="$DIR/presynaptic_network/correlations.csv"
# ./functions/compute_correlation_matrix_parallel.exe $folder1 $folder2
./functions/compute_tiling_correlation.exe $folder1 $folder2