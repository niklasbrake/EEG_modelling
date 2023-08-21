#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --account=def-akhadra
#SBATCH --mem=4G
#SBATCH --mail-user=niklas.brake@mail.mcgill.ca
#SBATCH --mail-type=FAIL,END
#SBATCH --output=analysis.log

folder=$1
S=$2

module load matlab/2020a
matlab -nodisplay -r "analyze_simulations('$folder',$S)"
