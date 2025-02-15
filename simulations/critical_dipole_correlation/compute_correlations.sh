#!/bin/bash
#SBATCH --time=002:00:00
#SBATCH --account=def-akhadra
#SBATCH --mem=4G
#SBATCH --mail-user=niklas.brake@mail.mcgill.ca
#SBATCH --mail-type=FAIL
#SBATCH --cpus-per-task=8
#SBATCH --output=correlation.log
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

folder=$1
subfolders=($(ls -d $folder/*))

for i in {0..9}; do
    b=$((10*(${SLURM_ARRAY_TASK_ID}-1)+i))
    tempFolder=${subfolders[b]}

    spikingFile=$tempFolder'/presynaptic_network/spikeTimesLong.csv'
    correlationFile=$tempFolder'/presynaptic_network/correlations.csv'

    /lustre04/scratch/nbrake/code/simulation_code/compute_tiling_correlation.exe $spikingFile $correlationFile
done