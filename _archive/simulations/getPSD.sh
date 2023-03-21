#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --gpus-per-node=1
#SBATCH --account=def-akhadra
#SBATCH --mem=4G
#SBATCH --mail-user=niklas.brake@mail.mcgill.ca
#SBATCH --mail-type=ALL

module load matlab/2020a
matlab -nodisplay -r "temp_sim_analysis"