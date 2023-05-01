#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --gpus-per-node=1
#SBATCH --account=def-akhadra
#SBATCH --mem=16G
#SBATCH --mail-user=niklas.brake@mail.mcgill.ca
#SBATCH --mail-type=ALL

module load matlab/2020a

matlab -nodisplay -r "analyze_results"