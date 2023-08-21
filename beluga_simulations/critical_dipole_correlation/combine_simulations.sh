#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --mem=4GB
#SBATCH --account=def-akhadra
#SBATCH --mail-user=niklas.brake@mail.mcgill.ca
#SBATCH --mail-type=FAIL,END
#SBATCH --output=analysis.log

folder=$1

module load matlab/2020a
matlab -nodisplay -r "combine_simulations('$folder')"