#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --account=def-akhadra
#SBATCH --mem=4G
#SBATCH --mail-user=niklas.brake@mail.mcgill.ca
#SBATCH --mail-type=FAIL
#SBATCH --output=initialization.log

folder=$1
branchNo=$2
module load matlab/2020a
matlab -nodisplay -r "initialize_network('$folder',$branchNo,${SLURM_ARRAY_TASK_ID})"

