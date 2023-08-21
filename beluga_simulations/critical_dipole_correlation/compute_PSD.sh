#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --gpus-per-node=1
#SBATCH --account=def-akhadra
#SBATCH --mem=48G
#SBATCH --mail-user=niklas.brake@mail.mcgill.ca
#SBATCH --mail-type=FAIL
#SBATCH --output=analysis.log

folder=$1

module load matlab/2020a
matlab -nodisplay -r "compute_PSD('$folder')"