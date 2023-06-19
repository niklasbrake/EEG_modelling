#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --gpus-per-node=1
#SBATCH --account=def-akhadra
#SBATCH --mem=32G
#SBATCH --mail-user=niklas.brake@mail.mcgill.ca
#SBATCH --mail-type=FAIL,END
#SBATCH --output=analysis.log

module load matlab/2020a
matlab -nodisplay -r "analyze_simulations"