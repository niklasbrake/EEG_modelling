#!/bin/bash
#SBATCH --time=000:30:00
#SBATCH --account=def-akhadra
#SBATCH --mem=5000
#SBATCH --mail-user=niklas.brake@mail.mcgill.ca
#SBATCH --mail-type=FAIL,END
#SBATCH --cpus-per-task=12
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" folderList.txt)
folder1="$DIR/spikeFileLong.csv"
folder2="$DIR/correlationFile.csv"
./functions/compute_correlation_matrix_parallel.exe $folder1 $folder2