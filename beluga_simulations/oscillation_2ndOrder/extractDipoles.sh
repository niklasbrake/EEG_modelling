#!/bin/bash
#SBATCH --time=00:05:00
#SBATCH --account=def-akhadra
#SBATCH --mem=4G
#SBATCH --mail-user=niklas.brake@mail.mcgill.ca
#SBATCH --mail-type=FAIL,END

module load matlab/2020a
matlab -nodisplay -r "extractDipoles"

