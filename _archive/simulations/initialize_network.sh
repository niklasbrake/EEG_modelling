#!/bin/bash
#SBATCH --time=001:30:00
#SBATCH --account=def-akhadra
#SBATCH --mem=4000
#SBATCH --mail-user=niklas.brake@mail.mcgill.ca
#SBATCH --mail-type=FAIL

# DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $1)
DIR=$1
module load matlab/2020a
echo $DIR
matlab -nodisplay -r "initialize_network('$DIR')"
